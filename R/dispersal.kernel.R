#' Dispersal kernels for dispersal data
#'
#' Fits several pre-defined distributions to dispersal or movement data, computing several estimators: AIC, AICc, BIC,
#' Chi-squared, and Kolgomorov-Smirnov. It also estimates the parameter(s) value(s) and SE of each distribution,
#' as well as its Mean, Variance, Skewness, and Kurtosis.
#'
#' @param data A numeric vector of distances.
#' @param distribution A character string naming the distributions to fit. By default “all” distributions are selected,
#' but they may be selected manually by name. See Details.
#' @param order.by A character string giving the estimator by which the distributions are to be ordered in the output table.
#'  The default is “AICc”, but one can also choose “AIC” or “BIC”.
#' @param extreme.values Bolean (TRUE/FALSE) to whether “data” is to be fitted with extreme value distributions. Default is FALSE.
#' @details This function fits one or more dispersal kernels with 1-2 parameters, by estimating the
#' distribution of kernel parameters (θ∈R) maximizing the likelihood function
#' @return Returns a list with all maximum likelihood calculations for each selected distribution.
#' Furthermore, two summary tables may be accessed
#' \tabular{ll}{
#'   \code{distribution.selection} \tab A table listing the fitted distributions. \cr
#'   \code{distribution.parameters} \tab A table with the estimated parameters and moments for each fitted distribution. \cr
#' }
#'
#' @seealso \code{\link{plot.dispfit}}
#' @keywords dispersal kernel distribution
#' @import msm numDeriv stats
#' @export
#' @examples
#' ## simulate data from exponential distribution
#' sim <- rexp(100, rate = 0.01)
#'
#' ## run dispersal.kernel function
#' fit <- dispersal.kernel(sim)
#'
#' ## display table
#' fit$distribution.selection

dispersal.kernel <- function (data, distribution = "all", order.by = "AICc", extreme.values = FALSE) {
  # require("msm")
  # require("numDeriv")
  kernel.fit <- list()
  values <- setNames(data.frame(matrix(ncol = 19, nrow = 1)),  c("AIC", "AICc", "BIC", "Chi-squared value", "Chi-squared significance",
                                                                 "Kolmogorov-Smirnov value", "K-S significance",
                                                                 "Parameter 1", "Parameter 1 SE",
                                                                 "Parameter 2", "Parameter 2 SE",
                                                                 "Mean", "Mean SE",
                                                                 "Variance", "Variance SE",
                                                                 "Skewness", "Skewness SE",
                                                                 "Kurtosis", "Kurtosis SE"))

  ## Saving data to final object - to be used in kernel plot function
  kernel.fit$data <- data
  ## Used in Chi-Squared calculations
  chi.res.hist <- hist(data, plot = F, breaks = "FD")
  ## Used in Kolmogorov-Smirnov calculations
  ks.res.hist <- hist(data, plot = F, breaks = length(data))
  # res.hist <- hist(data, plot = F, breaks = length(data))
  # cumulative.data <- c(res.hist$counts[1])
  # for (i in 1+seq_along(res.hist$counts)) {
  #   cumulative.data[i] <- cumulative.data[i-1] + res.hist$counts[i]
  # }
  # cumulative.data <- cumulative.data/sum(res.hist$counts)
  # cumulative.data <- cumulative.data[!is.na(cumulative.data)]

  if (isFALSE(extreme.values)) {
    if ("all" %in% distribution | "rayleigh" %in% distribution) {
      message(paste("Fitting Rayleigh distribution..."))
      rayleigh.values <- rayleigh.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$rayleigh <- rayleigh.values$opt
      values <- rbind(values, setNames(rayleigh.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Rayleigh"
    }
    if ("all" %in% distribution | "exponential" %in% distribution) {
      message(paste("Fitting Exponential distribution..."))
      exponential.values <- exponential.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$exponential <- exponential.values$opt
      values <- rbind(values, setNames(exponential.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Exponential"
    }
    if ("all" %in% distribution | "general normal" %in% distribution) {
      message(paste("Fitting Generalized Normal distribution..."))
      generalnormal.values <- generalnormal.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$generalnormal <- generalnormal.values$opt
      values <- rbind(values, setNames(generalnormal.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Generalized Normal"
    }
    ### 4 ### 2Dt / BIVARIATE STUDENT'S T
    if ("all" %in% distribution | "2Dt" %in% distribution) {
      message(paste("Fitting 2Dt distribution..."))
      twodt.values <- twodt.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$twodt <- twodt.values$opt
      values <- rbind(values, setNames(twodt.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "2Dt"
    }
    ### 5 ### POWER-LAW / GEOMETRIC
    if ("all" %in% distribution | "geometric" %in% distribution) {
      message(paste("Fitting Geometric distribution..."))
      geometric.values <- geometric.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$geometric <- geometric.values$opt
      values <- rbind(values, setNames(geometric.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Geometric"
    }
    ### 6 ### INVERSE POWER-LAW - UNDEFINED

    ### 7 ### LOGISTIC ###
    if ("all" %in% distribution | "logistic" %in% distribution) {
      message(paste("Fitting Logistic distribution..."))
      logistic.values <- logistic.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$logistic <- logistic.values$opt
      values <- rbind(values, setNames(logistic.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Logistic"
    }
    ### 8 ### LOGNORMAL ### DONE!
    if ("all" %in% distribution | "lognormal" %in% distribution) {
      message(paste("Fitting Log-Normal distribution..."))
      lognorm.values <- lognorm.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$lognorm <- lognorm.values$opt
      values <- rbind(values, setNames(lognorm.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Log-Normal"
    }
    ### 9 ### GAUSSIAN MIXTURE

    ### 10 ### WALD (INVERSE GAUSSIAN) ### DONE!
    if ("all" %in% distribution | "wald" %in% distribution) {
      message(paste("Fitting Wald distribution..."))
      wald.values <- wald.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$wald <- wald.values$opt
      values <- rbind(values, setNames(wald.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Wald"
    }
    ### 11 ### WEIBULL ### DONE!
    if ("all" %in% distribution | "weibull" %in% distribution) {
      message(paste("Fitting Weibull distribution..."))
      weibull.values <- weibull.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$weibull <- weibull.values$opt
      values <- rbind(values, setNames(weibull.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Weibull"
    }
    ### 12 ### GAMMA ### DONE!
    if ("all" %in% distribution | "gamma" %in% distribution) {
      message(paste("Fitting Gamma distribution..."))
      gamma.values <- gamma.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$gamma <- gamma.values$opt
      values <- rbind(values, setNames(gamma.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Gamma"
      # end
    }
    ### 13 ### LOG-SECH / HYPERBOLIC SECANT ### DONE!
    if ("all" %in% distribution | "log-sech" %in% distribution) {
      message(paste("Fitting Log-sech distribution..."))
      logsech.values <- logsech.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$logsech <- logsech.values$opt
      values <- rbind(values, setNames(logsech.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Log-sech"
      # end
    }
    ### 14 ### CAUCHY
    if ("all" %in% distribution | "cauchy" %in% distribution) {
      message(paste("Fitting Cauchy distribution..."))
      cauchy.values <- cauchy.function(kernel.fit$data, chi.res.hist, ks.res.hist)

      kernel.fit$cauchy <- cauchy.values$opt
      values <- rbind(values, setNames(cauchy.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Cauchy"
      # end
    }
  }  else if (isTRUE(extreme.values)) {
    ### EXTREME DISTRIBUTIONS FROM GARCIA AND BORDA DE AGUA 2017
    ### 15 ### GUMBEL
    if ("all" %in% distribution | "gumbel" %in% distribution) {
      message(paste("Fitting Gumbel distribution..."))
      # function
      log.dist.gumbel <- function (par, r) {
        a <- par[1] # location, u
        b <- par[2] # scale > 0
        fgumbel <- (1/b) * exp(-((r - a)/b + exp(-(r - a)/b))) ## pdf from Forbes 2011
        -sum(log(fgumbel)) ## Negative Log Likelihood
      }
      dist.gumbel <- function (r, a, b) {
        fgumbel <- (1/b) * exp(-((r - a)/b + exp(-(r - a)/b))) ## pdf from Forbes 2011
      }
      # initial values estimation
      m <- mean(data)
      n <- length(data)
      b1 <- m - sum(data * exp(-data)) / sum(exp(-data))
      b2 <- m - sum(data * exp(-data / b1)) / sum(exp(-data / b1))
      a <- -b2 * log((1/n) * sum(exp(-data / b2)))
      # optimization procedure
      dist.gumbel.opt <- optim (par = c(a, b2), ## valor inicial para o "a"
                                fn = log.dist.gumbel, ## função a minimizar
                                r = data, ## dados
                                method = "Nelder-Mead"#, ## método quando para estimar apenas um parametro (par)
                                #lower = c(0.001, 0.001), ## limite inferior do parametro
                                #upper = c(Inf, Inf)
      )
      kernel.fit$gumbel <- dist.gumbel.opt
      # output values
      # AIC
      aic.gumbel <- 2*length(dist.gumbel.opt$par) + 2 * dist.gumbel.opt$value
      # AICc
      aicc.gumbel <- aic.gumbel + (2 * length(dist.gumbel.opt$par)^2 + 2 * length(dist.gumbel.opt$par))/(length(data) - length(dist.gumbel.opt$par) - 1 )
      # BIC
      bic.gumbel <-  2 * dist.gumbel.opt$value + length(dist.gumbel.opt$par)*log(length(data))
      # Chi-squared
      chi.expected.values.gumbel <- dist.gumbel(chi.res.hist$mids, dist.gumbel.opt$par[1], dist.gumbel.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
      chi.squared.statistic.gumbel <- sum((chi.res.hist$counts - chi.expected.values.gumbel)^2 / chi.expected.values.gumbel)
      chi.squared.pvalue.gumbel <- 1-pchisq(chi.squared.statistic.gumbel, length(chi.res.hist$counts)-3)
      # Kolmogorov-Smirnov
      ks.expected.values.gumbel <- dist.gumbel(ks.res.hist$mids, dist.gumbel.opt$par[1], dist.gumbel.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
      simul.gumbel <- c()
      for (i in seq_along(ks.res.hist$mids)) {
        simul.gumbel <- c(simul.gumbel, rep(ks.res.hist$mids[i], round(ks.expected.values.gumbel[i], 0)))
      }
      ks.gumbel <- ks.test(data, simul.gumbel)
      g.max.gumbel <- as.numeric(ks.gumbel$statistic)
      KS.gumbel <- as.numeric(ks.gumbel$p.value)

      # cumulative.expected.values.gumbel <- c(expected.values.gumbel[1])
      # for (i in 1+seq_along(expected.values.gumbel)) {
      #   cumulative.expected.values.gumbel[i] <- cumulative.expected.values.gumbel[i-1] + expected.values.gumbel[i]
      # }
      # cumulative.expected.values.gumbel <- cumulative.expected.values.gumbel/sum(expected.values.gumbel)
      # cumulative.expected.values.gumbel <- cumulative.expected.values.gumbel[!is.na(cumulative.expected.values.gumbel)]
      # g.max.gumbel <- max(abs(cumulative.data - cumulative.expected.values.gumbel))
      # if (g.max.gumbel < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
      #   KS.gumbel <- "Accept"
      # } else {KS.gumbel <- "Reject"}
      # parameter estimate
      par.1.gumbel <- dist.gumbel.opt$par[1]
      par.2.gumbel <- dist.gumbel.opt$par[2]
      # parameter estimate standard error
      par.1.se.gumbel <- sqrt(diag(solve(numDeriv::hessian(log.dist.gumbel, x=dist.gumbel.opt$par, r=data))))[1]
      par.2.se.gumbel <- sqrt(diag(solve(numDeriv::hessian(log.dist.gumbel, x=dist.gumbel.opt$par, r=data))))[2]
      # mean dispersal distance
      mean.gumbel <- dist.gumbel.opt$par[1] + dist.gumbel.opt$par[2]*0.5772
      mean.stderr.gumbel <- msm::deltamethod(~ x1 + x2*0.5772, mean = dist.gumbel.opt$par, cov = solve(numDeriv::hessian(log.dist.gumbel, x=dist.gumbel.opt$par, r=data)) )
      # variance
      variance.gumbel <- (pi^2/6) * dist.gumbel.opt$par[2]^2
      variance.stderr.gumbel <- msm::deltamethod(~ (pi^2/6) * x2^2, mean = dist.gumbel.opt$par, cov = solve(numDeriv::hessian(log.dist.gumbel, x=dist.gumbel.opt$par, r=data)) )
      # skewness
      skewness.gumbel <- 1.14
      skewness.stderr.gumbel <- NA
      # kurtosis
      kurtosis.gumbel <- 2.4
      kurtosis.stderr.gumbel <- NA
      # output
      gumbel.values <- data.frame(aic.gumbel, aicc.gumbel, bic.gumbel,
                                  chi.squared.statistic.gumbel, chi.squared.pvalue.gumbel,g.max.gumbel, KS.gumbel,
                                  par.1.gumbel, par.1.se.gumbel, par.2.gumbel, par.2.se.gumbel,
                                  mean.gumbel, mean.stderr.gumbel, variance.gumbel, variance.stderr.gumbel,
                                  skewness.gumbel, skewness.stderr.gumbel, kurtosis.gumbel, kurtosis.stderr.gumbel)
      values <- rbind(values, setNames(gumbel.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Gumbel"
      # end
    }
    ### 16 ### FRECHET
    if ("all" %in% distribution | "frechet" %in% distribution) {
      message(paste("Fitting Frechet distribution..."))
      # function
      log.dist.frechet <- function (par, r) {
        a <- par[1] # shape, a
        b <- par[2] # scale, s
        ffrechet <- (a/b) * ((r/b)^(-1-a)) * exp(-(r/b)^-a)
        -sum(log(ffrechet)) ## Negative Log Likelihood
      }
      dist.frechet <- function (r, a, b) {
        ffrechet <- (a/b) * ((r/b)^(-1-a)) * exp(-(r/b)^-a)
      }
      # initial values estimation
      # optimization procedure
      dist.frechet.opt <- optim (par = c(1, 1), ## valor inicial para o "a"
                                 fn = log.dist.frechet, ## função a minimizar
                                 r = data, ## dados
                                 method = "Nelder-Mead" #, ## método quando para estimar apenas um parametro (par)
                                 #lower = c(0.0001, 0.0001), ## limite inferior do parametro
                                 #upper = c(Inf, Inf)
      )
      kernel.fit$frechet <- dist.frechet.opt
      # output values
      # AIC
      aic.frechet <- 2*length(dist.frechet.opt$par) + 2 * dist.frechet.opt$value
      # AICc
      aicc.frechet <- aic.frechet + (2 * length(dist.frechet.opt$par)^2 + 2 * length(dist.frechet.opt$par))/(length(data) - length(dist.frechet.opt$par) - 1 )
      # BIC
      bic.frechet <-  2 * dist.frechet.opt$value + length(dist.frechet.opt$par)*log(length(data))
      # Chi-squared
      chi.expected.values.frechet <- dist.frechet(chi.res.hist$mids, dist.frechet.opt$par[1], dist.frechet.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
      chi.squared.statistic.frechet <- sum((chi.res.hist$counts - chi.expected.values.frechet)^2 / chi.expected.values.frechet)
      chi.squared.pvalue.frechet <- 1-pchisq(chi.squared.statistic.frechet, length(chi.res.hist$counts)-3)
      # Kolmogorov-Smirnov
      ks.expected.values.frechet <- dist.frechet(ks.res.hist$mids, dist.frechet.opt$par[1], dist.frechet.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
      simul.frechet <- c()
      for (i in seq_along(ks.res.hist$mids)) {
        simul.frechet <- c(simul.frechet, rep(ks.res.hist$mids[i], round(ks.expected.values.frechet[i], 0)))
      }
      ks.frechet <- ks.test(data, simul.frechet)
      g.max.frechet <- as.numeric(ks.frechet$statistic)
      KS.frechet <- as.numeric(ks.frechet$p.value)

      # cumulative.expected.values.frechet <- c(expected.values.frechet[1])
      # for (i in 1+seq_along(expected.values.frechet)) {
      #   cumulative.expected.values.frechet[i] <- cumulative.expected.values.frechet[i-1] + expected.values.frechet[i]
      # }
      # cumulative.expected.values.frechet <- cumulative.expected.values.frechet/sum(expected.values.frechet)
      # cumulative.expected.values.frechet <- cumulative.expected.values.frechet[!is.na(cumulative.expected.values.frechet)]
      # g.max.frechet <- max(abs(cumulative.data - cumulative.expected.values.frechet))
      # if (g.max.frechet < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
      #   KS.frechet <- "Accept"
      # } else {KS.frechet <- "Reject"}
      # parameter estimate
      par.1.frechet <- dist.frechet.opt$par[1]
      par.2.frechet <- dist.frechet.opt$par[2]
      # parameter estimate standard error
      par.1.se.frechet <- sqrt(diag(solve(numDeriv::hessian(log.dist.frechet, x=dist.frechet.opt$par, r=data))))[1]
      par.2.se.frechet <- sqrt(diag(solve(numDeriv::hessian(log.dist.frechet, x=dist.frechet.opt$par, r=data))))[2]
      # mean dispersal distance
      if (dist.frechet.opt$par[1] > 1) {
        mean.frechet <- dist.frechet.opt$par[2]*gamma(1-1/dist.frechet.opt$par[1])
      } else {
        mean.frechet <- "Infinite Value"
      }
      if (dist.frechet.opt$par[1] > 1) {
        mean.stderr.frechet <-  msm::deltamethod(~ x2*gamma(1-1/x1), mean = dist.frechet.opt$par, cov = solve(numDeriv::hessian(log.dist.frechet, x=dist.frechet.opt$par, r=data)) )
      } else {
        mean.stderr.frechet <- "Infinite Value"
      }
      # variance
      if (dist.frechet.opt$par[1] > 2) {
        variance.frechet <- dist.frechet.opt$par[2]^2 * (gamma(1-2/dist.frechet.opt$par[1]) - (gamma(1-1/dist.frechet.opt$par[1]))^2)
      } else {
        variance.frechet <- "Infinite Value"
      }
      if (dist.frechet.opt$par[1] > 2) {
        variance.stderr.frechet <- msm::deltamethod(~ x2^2 * (gamma(1-2/x1) - (gamma(1-1/x1))^2),
                                                    mean = dist.frechet.opt$par,
                                                    cov = solve(numDeriv::hessian(log.dist.frechet, x=dist.frechet.opt$par, r=data)) )
      } else {
        variance.stderr.frechet <- "Infinite Value"
      }
      # skewness
      if (dist.frechet.opt$par[1] > 3) {
        skewness.frechet <- (gamma(1-3/dist.frechet.opt$par[1]) - 3*gamma(1-2/dist.frechet.opt$par[1]) * gamma(1-1/dist.frechet.opt$par[1]) + 2*gamma(1-1/dist.frechet.opt$par[1])^3) / (sqrt((gamma(1-2/dist.frechet.opt$par[1]) - gamma(1-1/dist.frechet.opt$par[1])^2)^3))
      } else {
        skewness.frechet <- "Infinite Value"
      }
      if (dist.frechet.opt$par[1] > 3) {
        skewness.stderr.frechet <- msm::deltamethod(~ (gamma(1-3/x1) - 3*gamma(1-2/x1) * gamma(1-1/x1) + 2*gamma(1-1/x1)^3) / (sqrt((gamma(1-2/x1) - gamma(1-1/x1)^2)^3)),
                                                    mean = dist.frechet.opt$par,
                                                    cov = solve(numDeriv::hessian(log.dist.frechet, x=dist.frechet.opt$par, r=data)) )
      } else {
        skewness.stderr.frechet <- "Infinite Value"
      }
      # kurtosis
      if (dist.frechet.opt$par[1] > 4) {
        kurtosis.frechet <- -6 + ((gamma(1-4/dist.frechet.opt$par[1]) - 4*gamma(1-3/dist.frechet.opt$par[1]) * gamma(1-1/dist.frechet.opt$par[1]) + 3*gamma(1-2/dist.frechet.opt$par[1])^2) / (gamma(1-2/dist.frechet.opt$par[1]) - gamma(1-1/dist.frechet.opt$par[1])^2)^2)
      } else {
        kurtosis.frechet <- "Infinite Value"
      }
      if (dist.frechet.opt$par[1] > 4) {
        kurtosis.stderr.frechet <- msm::deltamethod(~ -6 + ((gamma(1-4/x1) - 4*gamma(1-3/x1) * gamma(1-1/x1) + 3*gamma(1-2/x1)^2) / (gamma(1-2/x1) - gamma(1-1/x1)^2)^2),
                                                    mean = dist.frechet.opt$par,
                                                    cov = solve(numDeriv::hessian(log.dist.frechet, x=dist.frechet.opt$par, r=data)) )
      } else {
        kurtosis.stderr.frechet <- "Infinite Value"
      }
      # output
      frechet.values <- data.frame(aic.frechet, aicc.frechet, bic.frechet,
                                   chi.squared.statistic.frechet, chi.squared.pvalue.frechet,g.max.frechet, KS.frechet,
                                   par.1.frechet, par.1.se.frechet, par.2.frechet, par.2.se.frechet,
                                   mean.frechet, mean.stderr.frechet, variance.frechet, variance.stderr.frechet,
                                   skewness.frechet, skewness.stderr.frechet, kurtosis.frechet, kurtosis.stderr.frechet)
      values <- rbind(values, setNames(frechet.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Frechet"
      # end
    }
    ### 11 ### WEIBULL (EXTREME VALUES) ### DONE!
    if ("all" %in% distribution | "weibull" %in% distribution) {
      message(paste("Fitting Weibull distribution..."))
      # function
      log.dist.weibull <- function (r, par) {
        a <- par[1] ## scale
        b <- par[2] ## shape
        fw <- (b/a) * (r/a)^(b-1) * exp(-(r^b/a^b)) ## function from Austerlitz 2004
        -sum(log(fw)) ##
      }
      dist.weibull <- function (r, a, b) {
        fw <- (b/a) * (r/a)^(b-1) * exp(-(r^b/a^b)) ## function from Austerlitz 2004
      }
      # initial values estimation
      m <- mean(log(data))
      v <- var(log(data))
      shape <- 1.2/sqrt(v)
      scale <- exp(m + 0.572/shape)
      # optimization procedure
      dist.weibull.opt <- optim (par = c(scale, shape), ## valor inicial para o "a"
                                 fn = log.dist.weibull, ## função a minimizar
                                 r = data, ## dados
                                 method = "Nelder-Mead",
                                 # lower = c(0.00001, 0.00001),
                                 control = list(maxit = 10000))
      kernel.fit$weibull <- dist.weibull.opt
      # output values
      # AIC
      aic.weibull <- 2*length(dist.weibull.opt$par) + 2 * dist.weibull.opt$value
      # AICc
      aicc.weibull <- aic.weibull + (2 * length(dist.weibull.opt$par)^2 + 2 * length(dist.weibull.opt$par))/(length(data) - length(dist.weibull.opt$par) - 1 )
      # BIC
      bic.weibull <-  2 * dist.weibull.opt$value + length(dist.weibull.opt$par)*log(length(data))
      # Chi-squared
      chi.expected.values.weibull <- dist.weibull(chi.res.hist$mids, dist.weibull.opt$par[1], dist.weibull.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
      chi.squared.statistic.weibull <- sum((chi.res.hist$counts - chi.expected.values.weibull)^2 / chi.expected.values.weibull)
      chi.squared.pvalue.weibull <- 1-pchisq(chi.squared.statistic.weibull, length(chi.res.hist$counts)-3)
      # Kolmogorov-Smirnov
      ks.expected.values.weibull <- dist.weibull(ks.res.hist$mids, dist.weibull.opt$par[1], dist.weibull.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
      simul.weibull <- c()
      for (i in seq_along(ks.res.hist$mids)) {
        simul.weibull <- c(simul.weibull, rep(ks.res.hist$mids[i], round(ks.expected.values.weibull[i], 0)))
      }
      ks.weibull <- ks.test(data, simul.weibull)
      g.max.weibull <- as.numeric(ks.weibull$statistic)
      KS.weibull <- as.numeric(ks.weibull$p.value)

      # cumulative.expected.values.weibull <- c(expected.values.weibull[1])
      # for (i in 1+seq_along(expected.values.weibull)) {
      #   cumulative.expected.values.weibull[i] <- cumulative.expected.values.weibull[i-1] + expected.values.weibull[i]
      # }
      # cumulative.expected.values.weibull <- cumulative.expected.values.weibull/sum(expected.values.weibull)
      # cumulative.expected.values.weibull <- cumulative.expected.values.weibull[!is.na(cumulative.expected.values.weibull)]
      # g.max.weibull <- max(abs(cumulative.data - cumulative.expected.values.weibull))
      # if (g.max.weibull < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
      #   KS.weibull <- "Accept"
      # } else {KS.weibull <- "Reject"}
      # parameter estimate
      par.1.weibull <- dist.weibull.opt$par[1]
      par.2.weibull <- dist.weibull.opt$par[2]
      # parameter estimate standard error
      par.1.se.weibull <- sqrt(diag(solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data))))[1]
      par.2.se.weibull <- sqrt(diag(solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data))))[2]
      # mean dispersal distance ## from Austerlitz 2004
      mean.weibull <- dist.weibull.opt$par[1] * (gamma(1 + 1/dist.weibull.opt$par[2]))
      mean.stderr.weibull <- msm::deltamethod(~ x1 * (gamma(1 + 1/x2)), mean = dist.weibull.opt$par, cov = solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data)) )
      # variance
      variance.weibull <- dist.weibull.opt$par[1]^2*(gamma(1+2/dist.weibull.opt$par[2])-gamma(1+1/dist.weibull.opt$par[2])^2)
      variance.stderr.weibull <- msm::deltamethod(~ x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2), mean = dist.weibull.opt$par,
                                                  cov = solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data)) )
      # skewness
      skewness.weibull <- (dist.weibull.opt$par[1]^3*gamma(1+3/dist.weibull.opt$par[2])-3*mean.weibull*variance.weibull-mean.weibull^3)/variance.weibull^(3/2)
      skewness.stderr.weibull <- msm::deltamethod(~ (x1^3*gamma(1+3/x2)-3*(x1 * (gamma(1 + 1/x2)))*(x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))-(x1 * (gamma(1 + 1/x2)))^3)/(x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))^(3/2),
                                                  mean = dist.weibull.opt$par,
                                                  cov = solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data)) )
      # kurtosis
      kurtosis.weibull <- (dist.weibull.opt$par[1]^4*gamma(1+4/dist.weibull.opt$par[2])-4*skewness.weibull*variance.weibull^(3/2)*mean.weibull-6*variance.weibull*mean.weibull^2-mean.weibull^4)/variance.weibull^2-3
      kurtosis.stderr.weibull <- msm::deltamethod(~ (x1^4*gamma(1+4/x2) - 4 * ((x1^3*gamma(1+3/x2)-3*(x1 * (gamma(1 + 1/x2)))*(x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))-(x1 * (gamma(1 + 1/x2)))^3)/(x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))^(3/2)) * (x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))^(3/2) * (x1 * (gamma(1 + 1/x2))) - 6 * (x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2)) * (x1 * (gamma(1 + 1/x2)))^2 - (x1 * (gamma(1 + 1/x2)))^4) / (x1^2*(gamma(1+2/x2)-gamma(1+1/x2)^2))^2-3,
                                                  mean = dist.weibull.opt$par,
                                                  cov = solve(numDeriv::hessian(log.dist.weibull, x=dist.weibull.opt$par, r=data)) )
      # output
      weibull.values <- data.frame(aic.weibull, aicc.weibull, bic.weibull,
                                   chi.squared.statistic.weibull, chi.squared.pvalue.weibull,g.max.weibull, KS.weibull,
                                   par.1.weibull, par.1.se.weibull, par.2.weibull, par.2.se.weibull,
                                   mean.weibull, mean.stderr.weibull, variance.weibull, variance.stderr.weibull,
                                   skewness.weibull, skewness.stderr.weibull, kurtosis.weibull, kurtosis.stderr.weibull)
      values <- rbind(values, setNames(weibull.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Weibull"
      # end
    }
    # ### 15 ### GENERALIZED EXTREME VALUE (GEV)
    #   # function
    # log.dist.gumbel <- function (par, r) {
    #   a <- par[1]
    #   b <- par[2]
    #   fgumbel <- exp(-exp((a/b - r/b)) + a/b - r/b)/(2 * pi * b * r) ## pdf from Forbes 2011
    #   -sum(log(fgumbel)) ## Negative Log Likelihood
    # }
    # (1/c) * ((1 + (a * ((r-b)/c)))^(-1/a))^(a+1) * exp(-((1 + (a * ((r-b)/c)))^(-1/a)))
  } else {message("Please select TRUE/FALSE in extreme.values")}

  ### FINAL OUTPUTS ###
  kernel.fit$values <- values[-1,]

  for (i in seq_along(kernel.fit$values$AIC)) {
    kernel.fit$values$delta.AIC[i] <- kernel.fit$values$AIC[i]-min(kernel.fit$values$AIC)
  }
  for (i in seq_along(kernel.fit$values$AICc)) {
    kernel.fit$values$delta.AICc[i] <- kernel.fit$values$AICc[i]-min(kernel.fit$values$AICc)
  }
  for (i in seq_along(kernel.fit$values$BIC)) {
    kernel.fit$values$delta.BIC[i] <- kernel.fit$values$BIC[i]-min(kernel.fit$values$BIC)
  }
  for (i in seq_along(kernel.fit$values$delta.AICc)) {
    kernel.fit$values$wi[i] <- exp(-.5*kernel.fit$values$delta.AICc[i]) / sum(exp(-.5*kernel.fit$values$delta.AICc))
  }
  colnames(kernel.fit$values)[colnames(kernel.fit$values)=="delta.AIC"] <- "Delta AIC"
  colnames(kernel.fit$values)[colnames(kernel.fit$values)=="delta.AICc"] <- "Delta AICc"
  colnames(kernel.fit$values)[colnames(kernel.fit$values)=="delta.BIC"] <- "Delta BIC"
  colnames(kernel.fit$values)[colnames(kernel.fit$values)=="wi"] <- "wi"
  kernel.fit$values <- kernel.fit$values[,c(1, 20, 2, 21, 3, 22, 23, 4:19)]

  if ("AIC" %in% order.by) {
    kernel.fit$values <- kernel.fit$values[order(kernel.fit$values$AIC),]
  }
  if ("AICc" %in% order.by) {
    kernel.fit$values <- kernel.fit$values[order(kernel.fit$values$AICc),]
  }
  if ("wi" %in% order.by) {
    kernel.fit$values <- kernel.fit$values[order(kernel.fit$values$wi, decreasing = T),]
  }
  if ("BIC" %in% order.by) {
    kernel.fit$values <- kernel.fit$values[order(kernel.fit$values$BIC),]
  }

  # kernel.fit$values <- data.frame(lapply(kernel.fit$values[,c(1,15,2,16,3:14,17)],
  #                                               function(y) if(is.numeric(y)) round(y, 3) else y))

  kernel.fit$distribution.selection <- kernel.fit$values[, 1:11]
  kernel.fit$distribution.parameters <- kernel.fit$values[, 12:23]
  class(kernel.fit) <- "dispfit"
  kernel.fit
}


