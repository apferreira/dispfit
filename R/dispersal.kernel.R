#' Dispersal kernels for movement data
#'
#' Fits several pre-defined distributions to dispersal or movement data, computing several estimators: AIC, AICc, BIC, Chi-squared, and Kolgomorov-Smirnov. It also estimates the parameter(s) value(s) and SE of each distribution, as well as its Mean, Variance, Skewness, and Kurtosis.
#' @param data A numeric vector of distances.
#' @param distribution A character string naming the distributions to fit. By default “all” distributions are selected, but they may be selected manually by name. See Details.
#' @param order.by A character string giving the estimator by which the distributions are to be ordered in the output table. The default is “AICc”, but one can also choose “AIC” or “BIC”.
#' @param extreme.values Bolean (TRUE/FALSE) to whether “data” is to be fitted with extreme value distributions. Default is FALSE.
#' @return
#' @keywords dispersal kernel distribution
#' @import msm
#' @import numDeriv
#' @export
#' @examples
#' simulated.data.exponential <- rexp(100, rate = 0.01)
#' dispersal.kernel(simulated.data.exponential, distribution = "all", order.by = "AICc", extreme.values = FALSE)

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
  kernel.fit$data <- data
  res.hist <- hist(data, plot = F)
  cumulative.data <- c(res.hist$counts[1])
  for (i in 1+seq_along(res.hist$counts)) {
    cumulative.data[i] <- cumulative.data[i-1] + res.hist$counts[i]
  }
  cumulative.data <- cumulative.data/sum(res.hist$counts)
  cumulative.data <- cumulative.data[!is.na(cumulative.data)]

  if (isFALSE(extreme.values)) {
    ### 1 ### GAUSSIAN / RAYLEIGH ### DONE!
    if ("all" %in% distribution | "rayleigh" %in% distribution) {
      print(paste("Fitting Rayleigh distribution..."))
      # function
      log.dist.rayleigh <- function (r, a) {
        fg <-(1/(pi*a^2)) * exp(-r^2/a^2) ## Rayleigh, as defined in Nathan 2012
        -sum(log(fg)) ## Negative Log Likelihood
      }
      dist.rayleigh <- function (r, a) {
        fg <- 2*pi*r*(1/(pi*a^2)) * exp(-r^2/a^2)
      }
      # initial values estimation
      n <- length(data)
      sd0 <- sqrt((n - 1)/n) * sd(data) ## parameter estimate to use as initial value
      # optimization procedure
      dist.rayleigh.opt <- optim (par = sd0, ## initial value
                                  fn = log.dist.rayleigh, ## function to minimize
                                  r = data, ## dados
                                  method = "Brent", ## one-parameter method
                                  lower = 0.000001, ## lower-bound for parameter
                                  upper = 100000)
      kernel.fit$rayleigh <- dist.rayleigh.opt

      # output values
      # AIC
      aic.rayleigh <- 2 * length(dist.rayleigh.opt$par) + 2 * dist.rayleigh.opt$value
      # AICc
      aicc.rayleigh <- aic.rayleigh + (2 * length(dist.rayleigh.opt$par)^2 + 2 * length(dist.rayleigh.opt$par))/(length(data) - length(dist.rayleigh.opt$par) - 1 )
      # BIC
      bic.rayleigh <-  2 * dist.rayleigh.opt$value + length(dist.rayleigh.opt$par)*log(length(data))
      # Chi-squared
      expected.values.rayleigh <- dist.rayleigh(res.hist$mids,dist.rayleigh.opt$par)*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.rayleigh <- sum((res.hist$counts - expected.values.rayleigh)^2 / expected.values.rayleigh)
      chi.squared.pvalue.rayleigh <- 1 - pchisq(chi.squared.statistic.rayleigh, length(res.hist$counts)-2)
      # Kolmogorov-Smirnov
      cumulative.expected.values.rayleigh <- c(expected.values.rayleigh[1])
      for (i in 1+seq_along(expected.values.rayleigh)) {
        cumulative.expected.values.rayleigh[i] <- cumulative.expected.values.rayleigh[i-1] + expected.values.rayleigh[i]
      }
      cumulative.expected.values.rayleigh <- cumulative.expected.values.rayleigh/sum(expected.values.rayleigh)
      cumulative.expected.values.rayleigh <- cumulative.expected.values.rayleigh[!is.na(cumulative.expected.values.rayleigh)]
      g.max.rayleigh <- max(abs(cumulative.data - cumulative.expected.values.rayleigh))
      if (g.max.rayleigh < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.rayleigh <- "Accept"
      } else {KS.rayleigh <- "Reject"}
      # parameter estimate
      par.1.rayleigh <- dist.rayleigh.opt$par
      par.2.rayleigh <- NA
      # parameter estimate standard error
      par.1.se.rayleigh <- sqrt(diag(solve(numDeriv::hessian(log.dist.rayleigh, x=dist.rayleigh.opt$par, r=data))))
      # par.1.se.rayleigh <- sd(replicate(1000, optim (par = sd0, ## initial value
      #                       fn = log.dist.rayleigh, ## function to minimize
      #                       r = sample(data, replace = T), ## dados
      #                       method = "Brent", ## one-parameter method
      #                       lower = 0.000001, ## lower-bound for parameter
      #                       upper = 100000)$par))
      par.2.se.rayleigh <- NA
      # mean
      mean.rayleigh <- dist.rayleigh.opt$par*sqrt(pi)/2
      mean.stderr.rayleigh <- msm::deltamethod(~ x1 * sqrt(pi) / 2, mean = dist.rayleigh.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.rayleigh.opt$par, r=data)))
      # variance
      variance.rayleigh <- ((4-pi)*(dist.rayleigh.opt$par^2))/4
      variance.stderr.rayleigh <- msm::deltamethod(~ ((4-pi)*(x1^2))/4, mean = dist.rayleigh.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.rayleigh.opt$par, r=data)))
      # standard deviation
      stdev.rayleigh <- sqrt(((4-pi)*(dist.rayleigh.opt$par^2))/4)
      stdev.stderr.rayleigh <- msm::deltamethod(~ sqrt(((4-pi)*(x1^2))/4), mean = dist.rayleigh.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.rayleigh.opt$par, r=data)))
      # skewness
      skewness.rayleigh <- ((3*sqrt(pi)*dist.rayleigh.opt$par^3)/4)/variance.rayleigh^(3/2)
      skewness.stderr.rayleigh <- msm::deltamethod(~ (( (3*sqrt(pi)*x1^3) /4)) / (( (4-pi) * (x1^2) )/4)^(3/2), mean = dist.rayleigh.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.rayleigh.opt$par, r=data)))
      # kurtosis
      kurtosis.rayleigh <- (2*dist.rayleigh.opt$par^4)/variance.rayleigh^2
      kurtosis.stderr.rayleigh <- msm::deltamethod(~ (2*x1^4) / (( (4-pi) * (x1^2) )/4)^(2), mean = dist.rayleigh.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.rayleigh.opt$par, r=data)))
      # output
      rayleigh.values <- data.frame(aic.rayleigh, aicc.rayleigh, bic.rayleigh,
                                    chi.squared.statistic.rayleigh, chi.squared.pvalue.rayleigh,g.max.rayleigh, KS.rayleigh,
                                    par.1.rayleigh, par.1.se.rayleigh, par.2.rayleigh, par.2.se.rayleigh,
                                    mean.rayleigh, mean.stderr.rayleigh, variance.rayleigh, variance.stderr.rayleigh,
                                    skewness.rayleigh, skewness.stderr.rayleigh, kurtosis.rayleigh, kurtosis.stderr.rayleigh)
      values <- rbind(values, setNames(rayleigh.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Rayleigh"
      # end
    }
    if ("all" %in% distribution | "exponential" %in% distribution) {
      print(paste("Fitting Exponential distribution..."))
      ### 2 ### EXPONENTIAL ### DONE!
      # function
      log.dist.exponential <- function (r, a) {
        fexponential <- (1 / (2 * pi * a * r )) * exp(-r/a) # corrected function, adapted from Nathan 2012
        -sum(log(fexponential))
      }
      dist.exponential <- function (r, a) {
        fexponential <-  2*pi*r*(1 / (2 * pi * a * r )) * exp(-r/a) # corrected function, adapted from Nathan 2012
      }
      # initial values estimation
      rate <- 1/mean(data)
      # optimization procedure
      dist.exponential.opt <- optim (par = rate,
                                     fn = log.dist.exponential,
                                     r = data,
                                     method = "Brent",
                                     lower = 0.000001,
                                     upper = 100000)
      kernel.fit$exponential <- dist.exponential.opt
      # output values
      # AIC
      aic.exponential <- 2 + 2 * dist.exponential.opt$value
      # AICc
      aicc.exponential <- aic.exponential + (2 * length(dist.exponential.opt$par)^2 + 2 * length(dist.exponential.opt$par))/(length(data) - length(dist.exponential.opt$par) - 1 )
      # BIC
      bic.exponential <-  2 * dist.exponential.opt$value + length(dist.exponential.opt$par)*log(length(data))
      # Chi-squared
      expected.values.exponential <- dist.exponential(res.hist$mids, dist.exponential.opt$par)*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.exponential <- sum((res.hist$counts - expected.values.exponential)^2 / expected.values.exponential)
      chi.squared.pvalue.exponential <- 1-pchisq(chi.squared.statistic.exponential, length(res.hist$counts)-2)
      # Kolmogorov-Smirnov
      cumulative.expected.values.exponential <- c(expected.values.exponential[1])
      for (i in 1+seq_along(expected.values.exponential)) {
        cumulative.expected.values.exponential[i] <- cumulative.expected.values.exponential[i-1] + expected.values.exponential[i]
      }
      cumulative.expected.values.exponential <- cumulative.expected.values.exponential/sum(expected.values.exponential)
      cumulative.expected.values.exponential <- cumulative.expected.values.exponential[!is.na(cumulative.expected.values.exponential)]
      g.max.exponential <- max(abs(cumulative.data - cumulative.expected.values.exponential))
      if (g.max.exponential < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.exponential <- "Accept"
      } else {KS.exponential <- "Reject"}
      # parameter estimate
      par.1.exponential <- dist.exponential.opt$par
      par.2.exponential <- NA
      # parameter estimate standard error
      par.1.se.exponential <- sqrt(diag(solve(numDeriv::hessian(log.dist.exponential, x=dist.exponential.opt$par, r=data))))
      par.2.se.exponential <- NA
      # mean dispersal distance
      mean.exponential <- dist.exponential.opt$par
      mean.stderr.exponential <- sqrt(diag(solve(numDeriv::hessian(log.dist.exponential, x=dist.exponential.opt$par, r=data))))
      # variance
      variance.exponential <- dist.exponential.opt$par^2
      variance.stderr.exponential <- msm::deltamethod(~ x1^2, mean = dist.exponential.opt$par, cov = solve(numDeriv::hessian(log.dist.exponential, x=dist.exponential.opt$par, r=data)))
      # standard deviation
      stdev.exponential <- dist.exponential.opt$par
      stdev.stderr.exponential <- sqrt(diag(solve(numDeriv::hessian(log.dist.exponential, x=dist.exponential.opt$par, r=data))))
      # skewness
      skewness.exponential <- 6
      skewness.stderr.exponential <- NA
      # kurtosis
      kurtosis.exponential <- 24
      kurtosis.stderr.exponential <- NA
      # output
      exponential.values <- data.frame(aic.exponential, aicc.exponential, bic.exponential,
                                       chi.squared.statistic.exponential, chi.squared.pvalue.exponential, g.max.exponential, KS.exponential,
                                       par.1.exponential, par.1.se.exponential, par.2.exponential, par.2.se.exponential,
                                       mean.exponential, mean.stderr.exponential, variance.exponential, variance.stderr.exponential,
                                       skewness.exponential, skewness.stderr.exponential, kurtosis.exponential, kurtosis.stderr.exponential)
      values <- rbind(values, setNames(exponential.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Exponential"
      # end
    }
    ## 3 ### EXPONENTIAL POWER / GENERALIZED NORMAL / GENERALIZED ERROR
    if ("all" %in% distribution | "general normal" %in% distribution) {
      print(paste("Fitting Generalized Normal distribution..."))
      # function
      log.dist.generalnormal <- function (r, par) {
        a <- par[1] ## scale
        b <- par[2] ## shape
        fgeneralnormal <- (b / (2 * pi * (a^2) * gamma(2 / b))) * exp(-(r^b / a^b))
        -sum(log(fgeneralnormal)) ##
      }
      dist.generalnormal <- function (r, a, b) {
        fgeneralnormal <- 2*pi*r*(b / (2 * pi * (a^2) * gamma(2 / b))) * exp(-(r / a) ^ b)
      }
      # log.dist.generalnormal <- function (r, par) {
      #   a <- par[1] # scale
      #   b <- par[2] # shape
      #   fgeneralnormal <- (1/(2*pi*r))*(1 / (2 * a * gamma(1 / b))) * exp(-(r / a) ^ b)
      #   -sum(log(fgeneralnormal)) ##
      # }
      # initial values estimation ##
      m <- mean(data)
      n <- length(data)
      s1 <- m
      s2 <- sum((data-m)^2)/n
      s3 <- (sum((data-m)^4))/(n*s2^2)
      # scale <- sqrt(s2 * gamma(2/shape)/gamma(4/shape))
      shape <- (gamma(1/s1)*gamma(5/s1))/(gamma(3/s1)^2) ## works
      # scale <- sqrt(s2/(gamma(3/shape)/gamma(1/shape))) ## works
      scale <- sqrt(s2*gamma(1/shape)/gamma(3/shape)) ## works
      # optimization procedure
      dist.generalnormal.opt <- optim (par = c(scale, shape), ## valor inicial para o "a"
                                       fn = log.dist.generalnormal, ## função a minimizar
                                       r = data,
                                       method = "Nelder-Mead",
                                       # lower = c(0.000001, 0.00001),
                                       #upper = c(Inf, Inf),
                                       control = list(maxit = 10000))
      if (dist.generalnormal.opt$par[1] <= 0 | dist.generalnormal.opt$par[2] <= 0) {
        while (TRUE) {
          dist.generalnormal.opt <- optim (par = c(scale, shape), ## valor inicial para o "a"
                                           fn = log.dist.generalnormal, ## função a minimizar
                                           r = data,
                                           method = "SANN",
                                           # lower = c(0, 0),
                                           control = list(maxit = 10000))
          try.generalnormal <- try(
            dist.generalnormal.opt.try <- optim (par = c(dist.generalnormal.opt$par[1], dist.generalnormal.opt$par[2]), ## valor inicial para o "a"
                                                 fn = log.dist.generalnormal, ## função a minimizar
                                                 r = data,
                                                 method = "L-BFGS-B",
                                                 lower = c(0.000001, 0.000001),
                                                 upper = c(Inf, Inf),
                                                 control = list(maxit = 10000)),
            silent=T)
          if (class(try.generalnormal) != "try-error") {
            dist.generalnormal.opt.try
            break
          }
        }
        dist.generalnormal.opt <- dist.generalnormal.opt.try
      }
      kernel.fit$generalizednormal <- dist.generalnormal.opt
      #   # output values
      # AIC
      aic.generalnormal <- 2*length(dist.generalnormal.opt$par) + 2 * dist.generalnormal.opt$value ## AIC
      # AICc
      aicc.generalnormal <- aic.generalnormal + (2 * length(dist.generalnormal.opt$par)^2 + 2 * length(dist.generalnormal.opt$par))/(length(data) - length(dist.generalnormal.opt$par) - 1 )
      # BIC
      bic.generalnormal <-  2 * dist.generalnormal.opt$value + length(dist.generalnormal.opt$par)*log(length(data))
      # Chi-squared
      expected.values.generalnormal <- dist.generalnormal(res.hist$mids, dist.generalnormal.opt$par[1], dist.generalnormal.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.generalnormal <- sum((res.hist$counts - expected.values.generalnormal)^2 / expected.values.generalnormal)
      chi.squared.pvalue.generalnormal <- 1-pchisq(chi.squared.statistic.generalnormal, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.generalnormal <- c(expected.values.generalnormal[1])
      for (i in 1+seq_along(expected.values.generalnormal)) {
        cumulative.expected.values.generalnormal[i] <- cumulative.expected.values.generalnormal[i-1] + expected.values.generalnormal[i]
      }
      cumulative.expected.values.generalnormal <- cumulative.expected.values.generalnormal/sum(expected.values.generalnormal)
      cumulative.expected.values.generalnormal <- cumulative.expected.values.generalnormal[!is.na(cumulative.expected.values.generalnormal)]
      g.max.generalnormal <- max(abs(cumulative.data - cumulative.expected.values.generalnormal))
      if (g.max.generalnormal < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.generalnormal <- "Accept"
      } else {KS.generalnormal <- "Reject"}
      # mean dispersal distance
      # mean.generalnormal <- dist.generalnormal.opt$par[1] * (gamma(3/dist.generalnormal.opt$par[2])/gamma(2/dist.generalnormal.opt$par[2]))
      # parameter estimate
      par.1.generalnormal <- dist.generalnormal.opt$par[1]
      par.2.generalnormal <- dist.generalnormal.opt$par[2]
      # parameter estimate standard error
      par.1.se.generalnormal <- sqrt(diag(solve(numDeriv::hessian(log.dist.generalnormal, x=dist.generalnormal.opt$par, r=data))))[1]
      par.2.se.generalnormal <- sqrt(diag(solve(numDeriv::hessian(log.dist.generalnormal, x=dist.generalnormal.opt$par, r=data))))[2]
      # mean dispersal distance
      mean.generalnormal <- dist.generalnormal.opt$par[1] * (gamma(3/dist.generalnormal.opt$par[2])/gamma(2/dist.generalnormal.opt$par[2]))
      mean.stderr.generalnormal <- msm::deltamethod(~ x1 * (gamma(3/x2)/gamma(2/x2)),
                                                    mean = dist.generalnormal.opt$par,
                                                    cov = solve(numDeriv::hessian(log.dist.generalnormal,
                                                                                  x=dist.generalnormal.opt$par,
                                                                                  r=data)) )
      # variance
      variance.generalnormal <- (dist.generalnormal.opt$par[1]^2 *gamma(4/dist.generalnormal.opt$par[2])) / gamma(2/dist.generalnormal.opt$par[2]) - mean.generalnormal^2
      variance.stderr.generalnormal <- msm::deltamethod(~ (x1^2 *gamma(4/x2)) / gamma(2/x2) - (x1 * (gamma(3/x2)/gamma(2/x2)))^2,
                                                        mean = dist.generalnormal.opt$par,
                                                        cov = solve(numDeriv::hessian(log.dist.generalnormal,
                                                                                      x=dist.generalnormal.opt$par,
                                                                                      r=data)) )
      # standard deviation
      stdev.generalnormal <- sqrt((dist.generalnormal.opt$par[1]^2 *gamma(4/dist.generalnormal.opt$par[2])) / gamma(2/dist.generalnormal.opt$par[2]) - mean.generalnormal^2)
      stdev.stderr.generalnormal <- msm::deltamethod(~ sqrt((x1^2 *gamma(4/x2)) / gamma(2/x2) - (x1 * (gamma(3/x2)/gamma(2/x2)))^2),
                                                     mean = dist.generalnormal.opt$par,
                                                     cov = solve(numDeriv::hessian(log.dist.generalnormal,
                                                                                   x=dist.generalnormal.opt$par,
                                                                                   r=data)) )
      # skewness
      skewness.generalnormal <- ((dist.generalnormal.opt$par[1]^3 * gamma(5/dist.generalnormal.opt$par[2])) / gamma(2/dist.generalnormal.opt$par[2])) / (variance.generalnormal^(3/2))
      skewness.stderr.generalnormal <- msm::deltamethod(~ ((x1^3 * gamma(5/x2)) / gamma(2/x2)) / ((x1^2 *gamma(4/x2)) / gamma(2/x2) - (x1 * (gamma(3/x2)/gamma(2/x2)))^2)^(3/2),
                                                        mean = dist.generalnormal.opt$par,
                                                        cov = solve(numDeriv::hessian(log.dist.generalnormal,
                                                                                      x=dist.generalnormal.opt$par,
                                                                                      r=data)) )
      # kurtosis
      kurtosis.generalnormal <- ((dist.generalnormal.opt$par[1]^4 * gamma(6/dist.generalnormal.opt$par[2])) / gamma(2/dist.generalnormal.opt$par[2])) / (variance.generalnormal^2)
      # (gamma(6/dist.generalnormal.opt$par[2]) * gamma(2/dist.generalnormal.opt$par[2])) / (gamma(4/dist.generalnormal.opt$par[2])^2)
      # ((dist.generalnormal.opt$par[1]^4 * gamma(6/dist.generalnormal.opt$par[2])) / gamma(2/dist.generalnormal.opt$par[2])) / ((dist.generalnormal.opt$par[1]^2 *gamma(4/dist.generalnormal.opt$par[2])) / gamma(2/dist.generalnormal.opt$par[2]))^2
      kurtosis.stderr.generalnormal <- msm::deltamethod(~ ((x1^4 * gamma(6/x2)) / gamma(2/x2)) / ((x1^2 *gamma(4/x2)) / gamma(2/x2) - (x1 * (gamma(3/x2)/gamma(2/x2)))^2)^2,
                                                        mean = dist.generalnormal.opt$par,
                                                        cov = solve(numDeriv::hessian(log.dist.generalnormal,
                                                                                      x=dist.generalnormal.opt$par,
                                                                                      r=data)) )
      # output
      generalnormal.values <- data.frame(aic.generalnormal, aicc.generalnormal, bic.generalnormal,
                                         chi.squared.statistic.generalnormal, chi.squared.pvalue.generalnormal,g.max.generalnormal, KS.generalnormal,
                                         par.1.generalnormal, par.1.se.generalnormal, par.2.generalnormal, par.2.se.generalnormal,
                                         mean.generalnormal, mean.stderr.generalnormal, variance.generalnormal, variance.stderr.generalnormal,
                                         skewness.generalnormal, skewness.stderr.generalnormal, kurtosis.generalnormal, kurtosis.stderr.generalnormal)
      values <- rbind(values, setNames(generalnormal.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Generalized Normal"
      # end
    }
    ### 4 ### 2Dt / BIVARIATE STUDENT'S T
    if ("all" %in% distribution | "2Dt" %in% distribution) {
      print(paste("Fitting 2Dt distribution..."))
      # function
      log.dist.2dt <- function (r, par) {
        a <- par[1] ## scale parameter
        b <- par[2] ## shape parameter
        f2dt <- ((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))
        -sum(log(f2dt)) ##
      }
      dist.2dt <- function (r, a, b) {
        f2dt <- 2*pi*r*((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))
      }
      # initial values estimation
      while (TRUE) {
        SANN.2dt.opt <- optim (par = c(1, 1.000001), ## valor inicial para o "a"
                               fn = log.dist.2dt, ## função a minimizar
                               r = data,
                               method = "SANN",
                               # lower = c(0, 0),
                               control = list(maxit = 10000))
        try.2dt <- try(
          dist.2dt.opt.try <- optim (par = c(SANN.2dt.opt$par[1], SANN.2dt.opt$par[2]), ## valor inicial para o "a"
                                     fn = log.dist.2dt, ## função a minimizar
                                     r = data,
                                     method = "L-BFGS-B",
                                     lower = c(0.000001, 1.000001),
                                     upper = c(Inf, Inf),
                                     control = list(maxit = 10000)),
          silent=T)

        if (class(try.2dt) != "try-error") {
          dist.2dt.opt.try
          break
        }
      }
      # optimization procedure
      dist.2dt.opt <- dist.2dt.opt.try
      # dist.2dt.opt <- optim (par = c(0.001, 1.001), ## valor inicial para o "a"
      #                        fn = log.dist.2dt, ## função a minimizar
      #                        r = data, ## dados
      #                        control = list(maxit = 10000), ## limite superior do parametro
      #                        hessian = T)
      kernel.fit$twodt <- dist.2dt.opt
      # output values
      # AIC
      aic.2dt <- 2*length(dist.2dt.opt$par) + 2 * dist.2dt.opt$value
      # AICc
      aicc.2dt <- aic.2dt + (2 * length(dist.2dt.opt$par)^2 + 2 * length(dist.2dt.opt$par))/(length(data) - length(dist.2dt.opt$par) - 1 )
      # BIC
      bic.2dt <-  2 * dist.2dt.opt$value + length(dist.2dt.opt$par)*log(length(data))
      # Chi-squared
      expected.values.2dt <- dist.2dt(res.hist$mids, dist.2dt.opt$par[1], dist.2dt.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.2dt <- sum((res.hist$counts - expected.values.2dt)^2 / expected.values.2dt)
      chi.squared.pvalue.2dt <- 1-pchisq(chi.squared.statistic.2dt, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.2dt <- c(expected.values.2dt[1])
      for (i in 1+seq_along(expected.values.2dt)) {
        cumulative.expected.values.2dt[i] <- cumulative.expected.values.2dt[i-1] + expected.values.2dt[i]
      }
      cumulative.expected.values.2dt <- cumulative.expected.values.2dt/sum(expected.values.2dt)
      cumulative.expected.values.2dt <- cumulative.expected.values.2dt[!is.na(cumulative.expected.values.2dt)]
      g.max.2dt <- max(abs(cumulative.data - cumulative.expected.values.2dt))
      if (g.max.2dt < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.2dt <- "Accept"
      } else {KS.2dt <- "Reject"}
      # parameter estimate
      par.1.2dt <- dist.2dt.opt$par[1]
      par.2.2dt <- dist.2dt.opt$par[2]
      # parameter estimate standard error
      par.1.se.2dt <- sqrt(diag(solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data))))[1]
      par.2.se.2dt <- sqrt(diag(solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data))))[2]
      # mean dispersal distance
      if (dist.2dt.opt$par[2] >= 3/2) {
        mean.2dt <- dist.2dt.opt$par[1] * (sqrt(pi)/2) * (exp(lgamma(dist.2dt.opt$par[2]-(3/2))-lgamma(dist.2dt.opt$par[2]-1)))
      } else {
        mean.2dt <- "Infinite Value"
      }
      if (dist.2dt.opt$par[2] >= 1.5) {
        mean.stderr.2dt <- msm::deltamethod(~ x1 * (sqrt(pi)/2)*(exp(lgamma(x2-(3/2))-lgamma(x2-1))), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
      } else {
        mean.stderr.2dt <- "Infinite Value"
      }
      # variance
      if (dist.2dt.opt$par[2] > 2) {
        variance.2dt <- dist.2dt.opt$par[1] * sqrt(1/(dist.2dt.opt$par[2]-2))
      } else {
        variance.2dt <- "Infinite Value"
      }
      if (dist.2dt.opt$par[2] > 2) {
        variance.stderr.2dt <- msm::deltamethod(~ x1 * sqrt(1/(x2-2)), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
      } else {
        variance.stderr.2dt <- "Infinite Value"
      }
      # standard deviation
      if (dist.2dt.opt$par[2] > 2) {
        stdev.2dt <- sqrt(dist.2dt.opt$par[1] * sqrt(1/(dist.2dt.opt$par[2]-2)))
      } else {
        stdev.2dt <- "Infinite Value"
      }
      if (dist.2dt.opt$par[2] > 2) {
        stdev.stderr.2dt <- msm::deltamethod(~ sqrt(x1 * sqrt(1/(x2-2))), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
      } else {
        stdev.stderr.2dt <- "Infinite Value"
      }
      # skewness
      if (2*(dist.2dt.opt$par[2]-1) > 2.5) {
        skewness.2dt <- (dist.2dt.opt$par[1]*(gamma(2.5) * exp(lgamma(dist.2dt.opt$par[2]-2.5) - lgamma(dist.2dt.opt$par[2]-1)))^(1/3))
      } else {
        skewness.2dt <- "Infinite Value"
      }
      if (2*(dist.2dt.opt$par[2]-1) > 2.5) {
        skewness.stderr.2dt <- msm::deltamethod(~ x1 * (gamma(2.5) * exp(lgamma(x2-2.5) - lgamma(x2-1)))^(1/3), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
      } else {
        skewness.stderr.2dt <- "Infinite Value"
      }
      # kurtosis
      if (2*(dist.2dt.opt$par[2]-1) > 3) {
        kurtosis.2dt <- 2^(1/4) * dist.2dt.opt$par[1] * (1/((dist.2dt.opt$par[2] - 3) * (dist.2dt.opt$par[2] - 2)))^(1/4)
      } else {
        kurtosis.2dt <- "Infinite Value"
      }
      if (2*(dist.2dt.opt$par[2]-1) > 3) {
        kurtosis.stderr.2dt <- msm::deltamethod(~ 2^(1/4) * x1 * (1/((x2 - 3) * (x2 - 2)))^(1/4),
                                                mean = dist.2dt.opt$par,
                                                cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
      } else {
        kurtosis.stderr.2dt <- "Infinite Value"
      }
      # output
      twodt.values <- data.frame(aic.2dt, aicc.2dt, bic.2dt,
                                 chi.squared.statistic.2dt, chi.squared.pvalue.2dt,g.max.2dt, KS.2dt,
                                 par.1.2dt, par.1.se.2dt, par.2.2dt, par.2.se.2dt,
                                 mean.2dt, mean.stderr.2dt, variance.2dt, variance.stderr.2dt,
                                 skewness.2dt, skewness.stderr.2dt, kurtosis.2dt, kurtosis.stderr.2dt)
      values <- rbind(values, setNames(twodt.values, names(values)))
      row.names(values)[length(values$AIC)] <- "2Dt"
      # end
    }
    ### 5 ### POWER-LAW / GEOMETRIC
    if ("all" %in% distribution | "geometric" %in% distribution) {
      print(paste("Fitting Geometric distribution..."))
      # function
      log.dist.geometric <- function (r, par) {
        a <- par[1]
        b <- par[2]
        fgeometric <- (((b - 2) * (b - 1)) / (2 * pi * (a^2))) * ((1 + (r / a)) ^ -b)
        -sum(log(fgeometric)) ##
      }
      dist.geometric <- function (r, a, b) {
        fgeometric <- 2*pi*r*(((b - 2) * (b - 1)) / (2 * pi * (a^2))) * ((1 + (r / a)) ^ -b)
      }
      # initial values estimation
      while (TRUE) {
        SANN.geometric.opt <- optim (par = c(1, 2.000001), ## valor inicial para o "a"
                                     fn = log.dist.geometric, ## função a minimizar
                                     r = data,
                                     method = "SANN",
                                     # lower = c(0, 0),
                                     control = list(maxit = 10000))
        try.geometric <- try(
          dist.geometric.opt.try <- optim (par = c(SANN.geometric.opt$par[1], SANN.geometric.opt$par[2]), ## valor inicial para o "a"
                                           fn = log.dist.geometric, ## função a minimizar
                                           r = data,
                                           method = "L-BFGS-B",
                                           lower = c(0.000001, 2.000001),
                                           upper = c(Inf, Inf),
                                           control = list(maxit = 10000)),
          silent=T)

        if (class(try.geometric) != "try-error") {
          dist.geometric.opt.try
          break
        }
      }
      # optimization procedure
      dist.geometric.opt <- dist.geometric.opt.try
      # dist.geometric.opt <- optim (par = c(0.00001, 3.00001), ## valor inicial para o "a"
      #                       fn = log.dist.geometric, ## função a minimizar
      #                       r = data, ## dados
      #                       lower = c(0.00001, 3.00001),
      #                       control = list(maxit = 10000), ## limite superior do parametro
      #                       hessian = T)
      kernel.fit$geometric <- dist.geometric.opt
      # output values
      # AIC
      aic.geometric <- 2*length(dist.geometric.opt$par) + 2 * dist.geometric.opt$value
      # AICc
      aicc.geometric <- aic.geometric + (2 * length(dist.geometric.opt$par)^2 + 2 * length(dist.geometric.opt$par))/(length(data) - length(dist.geometric.opt$par) - 1 )
      # BIC
      bic.geometric <-  2 * dist.geometric.opt$value + length(dist.geometric.opt$par)*log(length(data))
      # Chi-squared
      expected.values.geometric <- dist.geometric(res.hist$mids, dist.geometric.opt$par[1], dist.geometric.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.geometric <- sum((res.hist$counts - expected.values.geometric)^2 / expected.values.geometric)
      chi.squared.pvalue.geometric <- 1-pchisq(chi.squared.statistic.geometric, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.geometric <- c(expected.values.geometric[1])
      for (i in 1+seq_along(expected.values.geometric)) {
        cumulative.expected.values.geometric[i] <- cumulative.expected.values.geometric[i-1] + expected.values.geometric[i]
      }
      cumulative.expected.values.geometric <- cumulative.expected.values.geometric/sum(expected.values.geometric)
      cumulative.expected.values.geometric <- cumulative.expected.values.geometric[!is.na(cumulative.expected.values.geometric)]
      g.max.geometric <- max(abs(cumulative.data - cumulative.expected.values.geometric))
      if (g.max.geometric < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.geometric <- "Accept"
      } else {KS.geometric <- "Reject"}
      # parameter estimate
      par.1.geometric <- dist.geometric.opt$par[1]
      par.2.geometric <- dist.geometric.opt$par[2]
      # parameter estimate standard error
      par.1.se.geometric <- sqrt(diag(solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data))))[1]
      par.2.se.geometric <- sqrt(diag(solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data))))[2]
      # mean dispersal distance
      if (dist.geometric.opt$par[2] >= 3) {
        mean.geometric <- (2 * dist.geometric.opt$par[1]) / (dist.geometric.opt$par[2]-3)
      } else {
        mean.geometric <- "Infinite Value"
      }
      if (dist.geometric.opt$par[2] >= 3) {
        mean.stderr.geometric <- msm::deltamethod(~ (2 * x1) / (x2-3),
                                                  mean = dist.geometric.opt$par,
                                                  cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
      } else {
        mean.stderr.geometric <- "Infinite Value"
      }
      # variance
      if (dist.geometric.opt$par[2] >= 4) {
        variance.geometric <- sqrt(6) * dist.geometric.opt$par[1] * sqrt(1/((dist.geometric.opt$par[2]-4) * (dist.geometric.opt$par[2]-3)))
      } else {
        variance.geometric <- "Infinite Value"
      }
      if (dist.geometric.opt$par[2] >= 4) {
        variance.stderr.geometric <- msm::deltamethod(~ sqrt(6) * x1 * sqrt(1/((x2-4) * (x2-3))),
                                                      mean = dist.geometric.opt$par,
                                                      cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
      } else {
        variance.stderr.geometric <- "Infinite Value"
      }
      # standard deviation
      if (dist.geometric.opt$par[2] >= 4) {
        stdev.geometric <- sqrt(sqrt(6) * dist.geometric.opt$par[1] * sqrt(1/((dist.geometric.opt$par[2]-4) * (dist.geometric.opt$par[2]-3))))
      } else {
        stdev.geometric <- "Infinite Value"
      }
      if (dist.geometric.opt$par[2] >= 4) {
        stdev.stderr.geometric <- msm::deltamethod(~ sqrt(sqrt(6) * x1 * sqrt(1/((x2-4) * (x2-3)))),
                                                   mean = dist.geometric.opt$par,
                                                   cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
      } else {
        stdev.stderr.geometric <- "Infinite Value"
      }
      # skewness
      if (dist.geometric.opt$par[2] >= 5) {
        skewness.geometric <- 2 * 3^(1/3) * dist.geometric.opt$par[1] * (1/((dist.geometric.opt$par[2] - 5) * (dist.geometric.opt$par[2] - 4) * (dist.geometric.opt$par[2] - 3)))^(1/3)
      } else {
        skewness.geometric <- "Infinite Value"
      }
      if (dist.geometric.opt$par[2] >= 5) {
        skewness.stderr.geometric <- msm::deltamethod(~ 2 * 3^(1/3) * x1 * (1/((x2 - 5) * (x2 - 4) * (x2 - 3)))^(1/3),
                                                      mean = dist.geometric.opt$par,
                                                      cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
      } else {
        skewness.stderr.geometric <- "Infinite Value"
      }
      # kurtosis
      if (dist.geometric.opt$par[2] >= 6) {
        kurtosis.geometric <- 2^(3/4) * 15^(1/4) * dist.geometric.opt$par[1] * (1/((dist.geometric.opt$par[2] - 6) * (dist.geometric.opt$par[2] - 5) * (dist.geometric.opt$par[2] - 4) * (dist.geometric.opt$par[2] - 3)))^(1/4)
      } else {
        kurtosis.geometric <- "Infinite Value"
      }
      if (dist.geometric.opt$par[2] >= 6) {
        kurtosis.stderr.geometric <- msm::deltamethod(~ 2^(3/4) * 15^(1/4) * x2 * (1/((x2 - 6) * (x2 - 5) * (x2 - 4) * (x2 - 3)))^(1/4),
                                                      mean = dist.geometric.opt$par,
                                                      cov = solve(numDeriv::hessian(log.dist.geometric, x=dist.geometric.opt$par, r=data)) )
      } else {
        kurtosis.geometric <- "Infinite Value"
      }
      # output
      geometric.values <- data.frame(aic.geometric, aicc.geometric, bic.geometric,
                                     chi.squared.statistic.geometric, chi.squared.pvalue.geometric,g.max.geometric, KS.geometric,
                                     par.1.geometric, par.1.se.geometric, par.2.geometric, par.2.se.geometric,
                                     mean.geometric, mean.stderr.geometric, variance.geometric, variance.stderr.geometric,
                                     skewness.geometric, skewness.stderr.geometric, kurtosis.geometric, kurtosis.stderr.geometric)
      values <- rbind(values, setNames(geometric.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Geometric"
      # end
    }
    ### 6 ### INVERSE POWER-LAW - UNDEFINED

    ### 7 ### LOGISTIC ###
    if ("all" %in% distribution | "logistic" %in% distribution) {
      print(paste("Fitting Logistic distribution..."))
      # function
      log.dist.logistic <- function (par, r) {
        a <- par[1] ## location, mean
        b <- par[2] ## scale
        flogistic <- (b / (2 * pi * (a^2) * gamma(2/b) * gamma(1-(2/b)) )) * ((1 + ((r^b) / (a^b)))^(-1))
        -sum(log(flogistic)) ## Negative Log Likelihood
      }
      dist.logistic <- function (r, a, b) {
        flogistic <- 2*pi*r*(b / (2 * pi * (a^2) * gamma(2/b) * gamma(1-(2/b)) )) * ((1 + ((r^b) / (a^b)))^(-1))
      }
      # initial values estimation
      # n <- length(data)
      # location <- mean(data)
      # v <- (n - 1)/n * var(data)
      # scale <- sqrt(3 * v)/pi

      while (TRUE) {
        SANN.logistic.opt <- optim (par = c(1, 3), ## valor inicial para o "a"
                                    fn = log.dist.logistic, ## função a minimizar
                                    r = data,
                                    method = "SANN",
                                    # lower = c(0, 0),
                                    control = list(maxit = 10000))
        try.logistic <- try(
          dist.logistic.opt.try <- optim (par = c(SANN.logistic.opt$par[1], SANN.logistic.opt$par[2]), ## valor inicial para o "a"
                                          fn = log.dist.logistic, ## função a minimizar
                                          r = data,
                                          method = "L-BFGS-B",
                                          lower = c(0.000001, 2.000001),
                                          upper = c(Inf, Inf),
                                          control = list(maxit = 10000)),
          silent=T)

        if (class(try.logistic) != "try-error") {
          dist.logistic.opt.try
          break
        }
      }
      # optimization procedure
      dist.logistic.opt <- dist.logistic.opt.try
      # dist.logistic.opt <- optim (par = c(location, scale), ## valor inicial para o "a"
      #                          fn = log.dist.logistic, ## função a minimizar
      #                          r = data, ## dados
      #                          method = "L-BFGS-B", ## método quando para estimar apenas um parametro (par)
      #                          lower = c(0.00001, 2.00001), ## limite inferior do parametro
      #                          upper = c(Inf, Inf), ## limite superior do parametro
      #                          hessian = T)
      kernel.fit$logistic <- dist.logistic.opt
      # output values
      # AIC
      aic.logistic <- 2 * length(dist.logistic.opt$par) + 2 * dist.logistic.opt$value
      # AICc
      aicc.logistic <- aic.logistic + (2 * length(dist.logistic.opt$par)^2 + 2 * length(dist.logistic.opt$par))/(length(data) - length(dist.logistic.opt$par) - 1 )
      # BIC
      bic.logistic <-  2 * dist.logistic.opt$value + length(dist.logistic.opt$par)*log(length(data))
      # Chi-squared
      expected.values.logistic <- dist.logistic(res.hist$mids, dist.logistic.opt$par[1], dist.logistic.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.logistic <- sum((res.hist$counts - expected.values.logistic)^2 / expected.values.logistic)
      chi.squared.pvalue.logistic <- 1-pchisq(chi.squared.statistic.logistic, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.logistic <- c(expected.values.logistic[1])
      for (i in 1+seq_along(expected.values.logistic)) {
        cumulative.expected.values.logistic[i] <- cumulative.expected.values.logistic[i-1] + expected.values.logistic[i]
      }
      cumulative.expected.values.logistic <- cumulative.expected.values.logistic/sum(expected.values.logistic)
      cumulative.expected.values.logistic <- cumulative.expected.values.logistic[!is.na(cumulative.expected.values.logistic)]
      g.max.logistic <- max(abs(cumulative.data - cumulative.expected.values.logistic))
      if (g.max.logistic < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.logistic <- "Accept"
      } else {KS.logistic <- "Reject"}
      # parameter estimate
      par.1.logistic <- dist.logistic.opt$par[1]
      par.2.logistic <- dist.logistic.opt$par[2]
      # parameter estimate standard error
      par.1.se.logistic <- sqrt(diag(solve(numDeriv::hessian(log.dist.logistic, x=dist.logistic.opt$par, r=data))))[1]
      par.2.se.logistic <- sqrt(diag(solve(numDeriv::hessian(log.dist.logistic, x=dist.logistic.opt$par, r=data))))[2]
      # mean dispersal distance
      mean.logistic <- dist.logistic.opt$par[1] * (( gamma(3/dist.logistic.opt$par[2]) * gamma(1-(3/dist.logistic.opt$par[2])) ) /
                                                     ( gamma(2/dist.logistic.opt$par[2]) * gamma(1-(2/dist.logistic.opt$par[2])) ) )

      mean.logistic <- dist.logistic.opt$par[1]^(1/dist.logistic.opt$par[2]) * (( gamma(3/dist.logistic.opt$par[2]) * gamma(1-(3/dist.logistic.opt$par[2])) ) /
                                                                                  ( gamma(2/dist.logistic.opt$par[2]) * gamma(1-(2/dist.logistic.opt$par[2])) ) )

      mean.stderr.logistic <- msm::deltamethod(~ x1 * ((gamma(3/x2) * gamma(1-(3/x2))) / (gamma(2/x2) * gamma(1-(2/x2))) ), mean = dist.logistic.opt$par, cov = solve(numDeriv::hessian(log.dist.logistic, x=dist.logistic.opt$par, r=data)) )
      # variance
      variance.logistic <- "in progress"
      # 1, 1/b, 1 + 1/b, -a^(-b) r^b
      # x <- 1000
      # par.2.logistic*x*hypergeo::hypergeo_buhring(1, 1/par.2.logistic, 1 + 1/par.2.logistic, -par.1.logistic^(-par.2.logistic)* x^par.2.logistic)
      #
      # ((x^2)/2)*(par.1.logistic^-par.2.logistic)*hypergeo::hypergeo_buhring(1, 2/par.2.logistic, (par.2.logistic+2)/par.2.logistic, -par.1.logistic^(-par.2.logistic) * x^par.2.logistic)
      #
      #
      #
      #
      # buhring_eqn11
      # hypergeo
      variance.stderr.logistic <-"in progress"
      # skewness
      skewness.logistic <- "in progress"
      skewness.stderr.logistic <- "in progress"
      # kurtosis
      kurtosis.logistic <- "in progress"
      kurtosis.stderr.logistic <- "in progress"
      # output
      logistic.values <- data.frame(aic.logistic, aicc.logistic, bic.logistic,
                                    chi.squared.statistic.logistic, chi.squared.pvalue.logistic,g.max.logistic, KS.logistic,
                                    par.1.logistic, par.1.se.logistic, par.2.logistic, par.2.se.logistic,
                                    mean.logistic, mean.stderr.logistic, variance.logistic, variance.stderr.logistic,
                                    skewness.logistic, skewness.stderr.logistic, kurtosis.logistic, kurtosis.stderr.logistic)
      values <- rbind(values, setNames(logistic.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Logistic"
      # end
    }
    ### 8 ### LOGNORMAL ### DONE!
    if ("all" %in% distribution | "lognormal" %in% distribution) {
      print(paste("Fitting Log-Normal distribution..."))
      # function
      log.dist.lognorm <- function (par, r) {
        a <- par[1]
        b <- par[2]
        flognorm <- (1 / (((2 * pi) ^ (3/2)) * (b * (r ^ 2)))) * exp(-(log(r / a)^2) / (2 * (b ^ 2)))
        -sum(log(flognorm)) ## Negative Log Likelihood da lognormiana
      }
      dist.lognorm <- function (r, a, b) {
        flognorm <- 2*pi*r * (1 / (((2 * pi) ^ (3/2)) * (b * (r ^ 2)))) * exp(-(log(r / a)^2) / (2 * (b ^ 2)))
      }
      # initial values estimation
      n <- length(data)
      lx <- log(data)
      sd0 <- sqrt((n - 1)/n) * sd(lx)
      ml <- mean(lx)
      # optimization procedure
      dist.lognorm.opt <- optim (par = c(ml, sd0), ## valor inicial para o "a"
                                 fn = log.dist.lognorm, ## função a minimizar
                                 r = data,
                                 method = "Nelder-Mead")
      kernel.fit$lognormal <- dist.lognorm.opt
      # output values
      # AIC
      aic.lognorm <- 2 * length(dist.lognorm.opt$par) + 2 * dist.lognorm.opt$value
      # AICc
      aicc.lognorm <- aic.lognorm + (2 * length(dist.lognorm.opt$par)^2 + 2 * length(dist.lognorm.opt$par))/(length(data) - length(dist.lognorm.opt$par) - 1 )
      # BIC
      bic.lognorm <-  2 * dist.lognorm.opt$value + length(dist.lognorm.opt$par)*log(length(data))
      # Chi-squared
      expected.values.lognorm <- dist.lognorm(res.hist$mids, dist.lognorm.opt$par[1], dist.lognorm.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.lognorm <- sum((res.hist$counts - expected.values.lognorm)^2 / expected.values.lognorm)
      chi.squared.pvalue.lognorm <- 1-pchisq(chi.squared.statistic.lognorm, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.lognorm <- c(expected.values.lognorm[1])
      for (i in 1+seq_along(expected.values.lognorm)) {
        cumulative.expected.values.lognorm[i] <- cumulative.expected.values.lognorm[i-1] + expected.values.lognorm[i]
      }
      cumulative.expected.values.lognorm <- cumulative.expected.values.lognorm/sum(expected.values.lognorm)
      cumulative.expected.values.lognorm <- cumulative.expected.values.lognorm[!is.na(cumulative.expected.values.lognorm)]
      g.max.lognorm <- max(abs(cumulative.data - cumulative.expected.values.lognorm))
      if (g.max.lognorm < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.lognorm <- "Accept"
      } else {KS.lognorm <- "Reject"}
      # parameter estimate
      par.1.lognorm <- dist.lognorm.opt$par[1]
      par.2.lognorm <- dist.lognorm.opt$par[2]
      # parameter estimate standard error
      par.1.se.lognorm <- sqrt(diag(solve(numDeriv::hessian(log.dist.lognorm, x=dist.lognorm.opt$par, r=data))))[1]
      par.2.se.lognorm <- sqrt(diag(solve(numDeriv::hessian(log.dist.lognorm, x=dist.lognorm.opt$par, r=data))))[2]
      # mean dispersal distance
      mean.lognorm <- dist.lognorm.opt$par[1] * exp((dist.lognorm.opt$par[2]^2)/2)
      mean.stderr.lognorm <- msm::deltamethod(~ x1 * exp((x2^2)/2), mean = dist.lognorm.opt$par, cov = solve(numDeriv::hessian(log.dist.lognorm, x=dist.lognorm.opt$par, r=data)) )
      # variance
      variance.lognorm <- (exp(dist.lognorm.opt$par[2]^2)-1) * dist.lognorm.opt$par[1]^2 * exp(dist.lognorm.opt$par[2]^2)
      variance.stderr.lognorm <- msm::deltamethod(~ (exp(x2^2)-1) * x1^2 * exp(x2^2), mean = dist.lognorm.opt$par, cov = solve(numDeriv::hessian(log.dist.lognorm, x=dist.lognorm.opt$par, r=data)))
      # skewness
      skewness.lognorm <- (exp(dist.lognorm.opt$par[2]^2)+2) * sqrt(exp(dist.lognorm.opt$par[2]^2)-1)
      skewness.stderr.lognorm <- msm::deltamethod(~ (exp(x2^2)+2) * sqrt(exp(x2^2)-1), mean = dist.lognorm.opt$par, cov = solve(numDeriv::hessian(log.dist.lognorm, x=dist.lognorm.opt$par, r=data)) )
      # kurtosis
      kurtosis.lognorm <- exp(4*dist.lognorm.opt$par[2]^2) + 2*exp(3*dist.lognorm.opt$par[2]^2) + 3*exp(2*dist.lognorm.opt$par[2]^2) - 6
      kurtosis.stderr.lognorm <- msm::deltamethod(~ exp(4*x2^2) + 2*exp(3*x2^2) + 3*exp(2*x2^2) - 6, mean = dist.lognorm.opt$par, cov = solve(numDeriv::hessian(log.dist.lognorm, x=dist.lognorm.opt$par, r=data)) )
      # output
      lognorm.values <- data.frame(aic.lognorm, aicc.lognorm, bic.lognorm,
                                   chi.squared.statistic.lognorm, chi.squared.pvalue.lognorm,g.max.lognorm, KS.lognorm,
                                   par.1.lognorm, par.1.se.lognorm, par.2.lognorm, par.2.se.lognorm,
                                   mean.lognorm, mean.stderr.lognorm, variance.lognorm, variance.stderr.lognorm,
                                   skewness.lognorm, skewness.stderr.lognorm, kurtosis.lognorm, kurtosis.stderr.lognorm)
      values <- rbind(values, setNames(lognorm.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Log-Normal"
      # end
    }
    ### 9 ### GAUSSIAN MIXTURE
    # log.dist.gaussianmix <- function (r, par) {
    #   a1 <- par[1]
    #   a2 <- par[2]
    #   p <- par[3]
    #   fgm <- ((p/(pi*a1^2)) * exp(-r^2/a1^2)) + ((1 - p/(pi*a2^2)) * exp(-r^2/a2^2))
    #   -sum(log(fgm)) ## Negative Log Likelihood da rayleighmixiana
    # }
    #
    # dist.gaussianmix.opt <- optim (par = c(0.001, 0.001, 0.001), ## valor inicial para o "a"
    #                             fn = log.dist.gaussianmix, ## função a minimizar
    #                             r = data, ## dados
    #                             method = "L-BFGS-B", ## método quando para estimar apenas um parametro (par)
    #                             lower = c(0, 0, 0), ## limite inferior do parametro
    #                             upper = c(Inf, Inf, 1)) ## limite superior do parametro
    #
    # aic.gaussianmix <- 2 - 2 * dist.rayleighmix.opt$value ## AIC
    #
    # dist.gaussianmix.opt$par*sqrt(pi)/2 ## mean dispersal distance

    # function
    # initial values estimation
    # optimization procedure

    # output values
    # AIC
    # AICc
    # mean dispersal distance
    # mean dispersal distance lower limit with 95% CI
    # mean dispersal distance upper limit with 95% CI

    ### 10 ### WALD (INVERSE GAUSSIAN) ### DONE!
    if ("all" %in% distribution | "wald" %in% distribution) {
      print(paste("Fitting Wald distribution..."))
      # function
      log.dist.wald <- function (r, par) {
        a <- par[1] ## location parameter, mean
        b <- par[2] ## scale parameter
        fwald <- (sqrt(b)/sqrt(8 * (pi^3) * (r^5))) * exp(-(b * ((r - a)^2))/(2 * (a^2) * r))
        -sum(log(fwald)) ##
      }
      dist.wald <- function (r, a, b) {
        fwald <- 2*pi*r * (sqrt(b)/sqrt(8 * (pi^3) * (r^5))) * exp(-(b * ((r - a)^2))/(2 * (a^2) * r))
      }
      # initial values estimation
      m <- mean(data)
      scale <- length(data)/(sum(1/data)-mean(1/data))
      # optimization procedure
      dist.wald.opt <- optim (par = c(m, scale), ##
                              fn = log.dist.wald, ##
                              r = data, ##
                              method = "Nelder-Mead",
                              # lower = c(0.00001, 0.00001), ## parameters minimum values
                              control = list(maxit = 10000))
      kernel.fit$wald <- dist.wald.opt
      # output values
      # AIC
      aic.wald <- 2*length(dist.wald.opt$par) + 2 * dist.wald.opt$value
      # AICc
      aicc.wald <- aic.wald + (2 * length(dist.wald.opt$par)^2 + 2 * length(dist.wald.opt$par))/(length(data) - length(dist.wald.opt$par) - 1 )
      # BIC
      bic.wald <-  2 * dist.wald.opt$value + length(dist.wald.opt$par)*log(length(data))
      # Chi-squared
      expected.values.wald <- dist.wald(res.hist$mids, dist.wald.opt$par[1], dist.wald.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.wald <- sum((res.hist$counts - expected.values.wald)^2 / expected.values.wald)
      chi.squared.pvalue.wald <- 1-pchisq(chi.squared.statistic.wald, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.wald <- c(expected.values.wald[1])
      for (i in 1+seq_along(expected.values.wald)) {
        cumulative.expected.values.wald[i] <- cumulative.expected.values.wald[i-1] + expected.values.wald[i]
      }
      cumulative.expected.values.wald <- cumulative.expected.values.wald/sum(expected.values.wald)
      cumulative.expected.values.wald <- cumulative.expected.values.wald[!is.na(cumulative.expected.values.wald)]
      g.max.wald <- max(abs(cumulative.data - cumulative.expected.values.wald))
      if (g.max.wald < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.wald <- "Accept"
      } else {KS.wald <- "Reject"}
      # parameter estimate
      par.1.wald <- dist.wald.opt$par[1]
      par.2.wald <- dist.wald.opt$par[2]
      # parameter estimate standard error
      par.1.se.wald <- sqrt(diag(solve(numDeriv::hessian(log.dist.wald, x=dist.wald.opt$par, r=data))))[1]
      par.2.se.wald <- sqrt(diag(solve(numDeriv::hessian(log.dist.wald, x=dist.wald.opt$par, r=data))))[2]
      # mean dispersal distance
      mean.wald <- dist.wald.opt$par[1]
      mean.stderr.wald <- par.1.se.wald
      # variance
      variance.wald <- dist.wald.opt$par[1]^3/dist.wald.opt$par[2]
      variance.stderr.wald <- msm::deltamethod(~ x1^3/x2, mean = dist.wald.opt$par, cov = solve(numDeriv::hessian(log.dist.wald, x=dist.wald.opt$par, r=data)) )
      # skewness
      skewness.wald <- 3 * sqrt((dist.wald.opt$par[1]*dist.wald.opt$par[2]))
      skewness.stderr.wald <- msm::deltamethod(~ 3 * sqrt((x1*x2)), mean = dist.wald.opt$par, cov = solve(numDeriv::hessian(log.dist.wald, x=dist.wald.opt$par, r=data)) )
      # kurtosis
      kurtosis.wald <- (15*dist.wald.opt$par[1])/dist.wald.opt$par[2]
      kurtosis.stderr.wald <- msm::deltamethod(~ (15*x1)/x2, mean = dist.wald.opt$par, cov = solve(numDeriv::hessian(log.dist.wald, x=dist.wald.opt$par, r=data)) )
      # output
      wald.values <- data.frame(aic.wald, aicc.wald, bic.wald,
                                chi.squared.statistic.wald, chi.squared.pvalue.wald,g.max.wald, KS.wald,
                                par.1.wald, par.1.se.wald, par.2.wald, par.2.se.wald,
                                mean.wald, mean.stderr.wald, variance.wald, variance.stderr.wald,
                                skewness.wald, skewness.stderr.wald, kurtosis.wald, kurtosis.stderr.wald)
      values <- rbind(values, setNames(wald.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Wald"
      # end
    }
    ### 11 ### WEIBULL ### DONE!
    if ("all" %in% distribution | "weibull" %in% distribution) {
      print(paste("Fitting Weibull distribution..."))
      # function
      log.dist.weibull <- function (r, par) {
        a <- par[1] ## scale
        b <- par[2] ## shape
        fw <- (b/(2*pi*a^b)) * (r^(b-2)) * exp(-(r^b/a^b)) ## function from Austerlitz 2004
        -sum(log(fw)) ##
      }
      dist.weibull <- function (r, a, b) {
        fw <- 2*pi*r * (b/(2*pi*a^b)) * (r^(b-2)) * exp(-(r^b/a^b)) ## function from Austerlitz 2004
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
      expected.values.weibull <- dist.weibull(res.hist$mids, dist.weibull.opt$par[1], dist.weibull.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.weibull <- sum((res.hist$counts - expected.values.weibull)^2 / expected.values.weibull)
      chi.squared.pvalue.weibull <- 1-pchisq(chi.squared.statistic.weibull, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.weibull <- c(expected.values.weibull[1])
      for (i in 1+seq_along(expected.values.weibull)) {
        cumulative.expected.values.weibull[i] <- cumulative.expected.values.weibull[i-1] + expected.values.weibull[i]
      }
      cumulative.expected.values.weibull <- cumulative.expected.values.weibull/sum(expected.values.weibull)
      cumulative.expected.values.weibull <- cumulative.expected.values.weibull[!is.na(cumulative.expected.values.weibull)]
      g.max.weibull <- max(abs(cumulative.data - cumulative.expected.values.weibull))
      if (g.max.weibull < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.weibull <- "Accept"
      } else {KS.weibull <- "Reject"}
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
    ### 12 ### GAMMA ### DONE!
    if ("all" %in% distribution | "gamma" %in% distribution) {
      print(paste("Fitting Gamma distribution..."))
      # function
      log.dist.gamma <- function (par, r) {
        a <- par[1] ## scale
        b <- par[2] ## shape
        fgamma <- (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
        -sum(log(fgamma)) ## Negative Log Likelihood
      }
      dist.gamma <- function (r, a, b) {
        fgamma <- 2*pi*r * (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
      }
      # initial values estimation
      shape <- (mean(data)/sd(data))
      scale <- sd(data)^2/mean(data)
      # optimization procedure
      dist.gamma.opt <- optim (par = c(scale, shape), ##
                               fn = log.dist.gamma, ## função a minimizar
                               r = data,
                               method = "L-BFGS-B",
                               lower = c(0.00001, 0.00001)
      )
      kernel.fit$gamma <- dist.gamma.opt
      # output values
      # AIC
      aic.gamma <- 2 * length(dist.gamma.opt$par) + 2 * dist.gamma.opt$value
      # AICc
      aicc.gamma <- aic.gamma + (2 * length(dist.gamma.opt$par)^2 + 2 * length(dist.gamma.opt$par))/(length(data) - length(dist.gamma.opt$par) - 1 )
      # BIC
      bic.gamma <-  2 * dist.gamma.opt$value + length(dist.gamma.opt$par)*log(length(data))
      # Chi-squared
      expected.values.gamma <- dist.gamma(res.hist$mids, dist.gamma.opt$par[1],  dist.gamma.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.gamma <- sum((res.hist$counts - expected.values.gamma)^2 / expected.values.gamma)
      chi.squared.pvalue.gamma <- 1-pchisq(chi.squared.statistic.gamma, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.gamma <- c(expected.values.gamma[1])
      for (i in 1+seq_along(expected.values.gamma)) {
        cumulative.expected.values.gamma[i] <- cumulative.expected.values.gamma[i-1] + expected.values.gamma[i]
      }
      cumulative.expected.values.gamma <- cumulative.expected.values.gamma/sum(expected.values.gamma)
      cumulative.expected.values.gamma <- cumulative.expected.values.gamma[!is.na(cumulative.expected.values.gamma)]
      g.max.gamma <- max(abs(cumulative.data - cumulative.expected.values.gamma))
      if (g.max.gamma < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.gamma <- "Accept"
      } else {KS.gamma <- "Reject"}
      # parameter estimate
      par.1.gamma <- dist.gamma.opt$par[1]
      par.2.gamma <- dist.gamma.opt$par[2]
      # parameter estimate standard error
      par.1.se.gamma <- sqrt(diag(solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data))))[1]
      par.2.se.gamma <- sqrt(diag(solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data))))[2]
      # mean dispersal distance
      mean.gamma <- dist.gamma.opt$par[1] * dist.gamma.opt$par[2]
      mean.stderr.gamma <- msm::deltamethod(~ x1 * x2,
                                            mean = dist.gamma.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data)) )
      # variance
      variance.gamma <- dist.gamma.opt$par[1]^2 * dist.gamma.opt$par[2]
      variance.stderr.gamma <- msm::deltamethod(~ x2*x1^2,
                                                mean = dist.gamma.opt$par,
                                                cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data)) )
      # skewness
      skewness.gamma <- 2/sqrt(dist.gamma.opt$par[2])
      skewness.stderr.gamma <- msm::deltamethod(~ 2/sqrt(x2),
                                                mean = dist.gamma.opt$par,
                                                cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data)) )
      # kurtosis
      kurtosis.gamma <- 6/dist.gamma.opt$par[2]
      kurtosis.stderr.gamma <- msm::deltamethod(~ 6/x2,
                                                mean = dist.gamma.opt$par,
                                                cov = solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data)) )
      # output
      gamma.values <- data.frame(aic.gamma, aicc.gamma, bic.gamma,
                                 chi.squared.statistic.gamma, chi.squared.pvalue.gamma,g.max.gamma, KS.gamma,
                                 par.1.gamma, par.1.se.gamma, par.2.gamma, par.2.se.gamma,
                                 mean.gamma, mean.stderr.gamma, variance.gamma, variance.stderr.gamma,
                                 skewness.gamma, skewness.stderr.gamma, kurtosis.gamma, kurtosis.stderr.gamma)
      values <- rbind(values, setNames(gamma.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Gamma"
      # end
    }
    ### 13 ### LOG-SECH / HYPERBOLIC SECANT ### DONE!
    if ("all" %in% distribution | "log-sech" %in% distribution) {
      print(paste("Fitting Log-sech distribution..."))
      # function
      log.dist.logsech <- function (par, r) {
        a <- par[1] ## location, median
        b <- par[2] ## scale
        flogsech <- (1 / ((pi^2) * b * (r^2))) / (((r / a)^(1 / b)) + ((r / a) ^ -(1 / b)))
        -sum(log(flogsech)) ##
      }
      dist.logsech <- function (r, a, b) {
        flogsech <- 2*pi*r * (1 / ((pi^2) * b * (r^2))) / (((r / a)^(1 / b)) + ((r / a) ^ -(1 / b)))
      }
      # initial values estimation # from Halley and Inchausti 2002
      location <- median(log(data))
      scale <- IQR(log(data))/(2*-log(tan(pi/8)))
      # optimization procedure
      dist.logsech.opt <- optim (par = c(location, scale), ## valor inicial para o "a"
                                 fn = log.dist.logsech, ## função a minimizar
                                 r = data, ## dados
                                 control = list(maxit = 10000),
                                 method = "Nelder-Mead")
      #lower = c(0.00001, 0.00001),

      kernel.fit$logsech <- dist.logsech.opt
      # output values
      # AIC
      aic.logsech <- 2*length(dist.logsech.opt$par) + 2 * dist.logsech.opt$value
      # AICc
      aicc.logsech <- aic.logsech + (2 * length(dist.logsech.opt$par)^2 + 2 * length(dist.logsech.opt$par))/(length(data) - length(dist.logsech.opt$par) - 1 )
      # BIC
      bic.logsech <-  2 * dist.logsech.opt$value + length(dist.logsech.opt$par)*log(length(data))
      # Chi-squared
      expected.values.logsech <- dist.logsech(res.hist$mids, dist.logsech.opt$par[1], dist.logsech.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.logsech <- sum((res.hist$counts - expected.values.logsech)^2 / expected.values.logsech)
      chi.squared.pvalue.logsech <- 1-pchisq(chi.squared.statistic.logsech, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.logsech <- c(expected.values.logsech[1])
      for (i in 1+seq_along(expected.values.logsech)) {
        cumulative.expected.values.logsech[i] <- cumulative.expected.values.logsech[i-1] + expected.values.logsech[i]
      }
      cumulative.expected.values.logsech <- cumulative.expected.values.logsech/sum(expected.values.logsech)
      cumulative.expected.values.logsech <- cumulative.expected.values.logsech[!is.na(cumulative.expected.values.logsech)]
      g.max.logsech <- max(abs(cumulative.data - cumulative.expected.values.logsech))
      if (g.max.logsech < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.logsech <- "Accept"
      } else {KS.logsech <- "Reject"}
      # parameter estimate
      par.1.logsech <- dist.logsech.opt$par[1]
      par.2.logsech <- dist.logsech.opt$par[2]
      # parameter estimate standard error
      par.1.se.logsech <- sqrt(diag(solve(numDeriv::hessian(log.dist.logsech, x=dist.logsech.opt$par, r=data))))[1]
      par.2.se.logsech <- sqrt(diag(solve(numDeriv::hessian(log.dist.logsech, x=dist.logsech.opt$par, r=data))))[2]
      # mean dispersal distance
      if (dist.logsech.opt$par[2] < 1) {
        mean.logsech <- 2*pi* (integrate(function(r) (r^2 / ((pi^2) * dist.logsech.opt$par[2] * (r^2))) / (((r / dist.logsech.opt$par[1])^(1 / dist.logsech.opt$par[2])) + ((r / dist.logsech.opt$par[1]) ^ -(1 / dist.logsech.opt$par[2]))),
                                         lower = 0,
                                         upper = Inf))$value
      } else {
        mean.logsech <- "Infinite Value"
      }
      if (dist.logsech.opt$par[2] < 1) {
        mean.stderr.logsech <- "in progress"
      } else {
        mean.stderr.logsech <- "Infinite Value"
      }
      # variance
      variance.logsech <- "in progress"
      variance.stderr.logsech <-"in progress"
      # skewness
      skewness.logsech <- "in progress"
      skewness.stderr.logsech <- "in progress"
      # kurtosis
      kurtosis.logsech <- "in progress"
      kurtosis.stderr.logsech <- "in progress"
      # output
      logsech.values <- data.frame(aic.logsech, aicc.logsech, bic.logsech,
                                   chi.squared.statistic.logsech, chi.squared.pvalue.logsech,g.max.logsech, KS.logsech,
                                   par.1.logsech, par.1.se.logsech, par.2.logsech, par.2.se.logsech,
                                   mean.logsech, mean.stderr.logsech, variance.logsech, variance.stderr.logsech,
                                   skewness.logsech, skewness.stderr.logsech, kurtosis.logsech, kurtosis.stderr.logsech)
      values <- rbind(values, setNames(logsech.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Log-sech"
      # end
    }
    ### 14 ### CAUCHY
    if ("all" %in% distribution | "cauchy" %in% distribution) {
      print(paste("Fitting Cauchy distribution..."))
      # function
      log.dist.cauchy <- function (par, r) {
        a <- par[1]
        b <- par[2]
        fcauchy <- (1/(2*pi*r))*(1 / (pi*b)) / (1 + ((r-a)/b)^2)
        -sum(log(fcauchy)) ##
      }
      dist.cauchy <- function (r, a, b) {
        fcauchy <- 2*pi*r * (1/(2*pi*r))*(1 / (pi*b)) / (1 + ((r-a)/b)^2)
      }
      # initial values estimation # from Halley and Inchausti 2002
      location <- median(data)
      scale <- IQR(data)/2
      # optimization procedure
      dist.cauchy.opt <- optim (par = c(location, scale), ## valor inicial para o "a"
                                fn = log.dist.cauchy, ## função a minimizar
                                r = data, ## dados
                                #control = list(maxit = 10000), ## limite superior do parametro
                                method = "Nelder-Mead"
                                # lower = c(0.00001, 0.00001)
      )
      kernel.fit$cauchy <- dist.cauchy.opt
      # output values
      # AIC
      aic.cauchy <- 2*length(dist.cauchy.opt$par) + 2 * dist.cauchy.opt$value
      # AICc
      aicc.cauchy <- aic.cauchy + (2 * length(dist.cauchy.opt$par)^2 + 2 * length(dist.cauchy.opt$par))/(length(data) - length(dist.cauchy.opt$par) - 1 )
      # BIC
      bic.cauchy <-  2 * dist.cauchy.opt$value + length(dist.cauchy.opt$par)*log(length(data))
      # Chi-squared
      expected.values.cauchy <- dist.cauchy(res.hist$mids, dist.cauchy.opt$par[1], dist.cauchy.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.cauchy <- sum((res.hist$counts - expected.values.cauchy)^2 / expected.values.cauchy)
      chi.squared.pvalue.cauchy <- 1-pchisq(chi.squared.statistic.cauchy, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.cauchy <- c(expected.values.cauchy[1])
      for (i in 1+seq_along(expected.values.cauchy)) {
        cumulative.expected.values.cauchy[i] <- cumulative.expected.values.cauchy[i-1] + expected.values.cauchy[i]
      }
      cumulative.expected.values.cauchy <- cumulative.expected.values.cauchy/sum(expected.values.cauchy)
      cumulative.expected.values.cauchy <- cumulative.expected.values.cauchy[!is.na(cumulative.expected.values.cauchy)]
      g.max.cauchy <- max(abs(cumulative.data - cumulative.expected.values.cauchy))
      if (g.max.cauchy < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.cauchy <- "Accept"
      } else {KS.cauchy <- "Reject"}
      # parameter estimate
      par.1.cauchy <- dist.cauchy.opt$par[1]
      par.2.cauchy <- dist.cauchy.opt$par[2]
      # parameter estimate standard error
      par.1.se.cauchy <- sqrt(diag(solve(numDeriv::hessian(log.dist.cauchy, x=dist.cauchy.opt$par, r=data))))[1]
      par.2.se.cauchy <- sqrt(diag(solve(numDeriv::hessian(log.dist.cauchy, x=dist.cauchy.opt$par, r=data))))[2]
      # mean dispersal distance # the mean is undefined, so we used the median in the case of Cauchy
      mean.cauchy <- par.1.cauchy
      mean.stderr.cauchy <- par.1.se.cauchy
      # variance
      variance.cauchy <- "undefined"
      variance.stderr.cauchy <- "undefined"
      # skewness
      skewness.cauchy <- "undefined"
      skewness.stderr.cauchy <- "undefined"
      # kurtosis
      kurtosis.cauchy <- "undefined"
      kurtosis.stderr.cauchy <- "undefined"
      # output
      cauchy.values <- data.frame(aic.cauchy, aicc.cauchy, bic.cauchy,
                                  chi.squared.statistic.cauchy, chi.squared.pvalue.cauchy,g.max.cauchy, KS.cauchy,
                                  par.1.cauchy, par.1.se.cauchy, par.2.cauchy, par.2.se.cauchy,
                                  mean.cauchy, mean.stderr.cauchy, variance.cauchy, variance.stderr.cauchy,
                                  skewness.cauchy, skewness.stderr.cauchy, kurtosis.cauchy, kurtosis.stderr.cauchy)
      values <- rbind(values, setNames(cauchy.values, names(values)))
      row.names(values)[length(values$AIC)] <- "Cauchy"
      # end
    }
  }  else if (isTRUE(extreme.values)) {
    ### EXTREME DISTRIBUTIONS FROM GARCIA AND BORDA DE AGUA 2017
    ### 15 ### GUMBEL
    if ("all" %in% distribution | "gumbel" %in% distribution) {
      print(paste("Fitting Gumbel distribution..."))
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
      expected.values.gumbel <- dist.gumbel(res.hist$mids, dist.gumbel.opt$par[1], dist.gumbel.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.gumbel <- sum((res.hist$counts - expected.values.gumbel)^2 / expected.values.gumbel)
      chi.squared.pvalue.gumbel <- 1-pchisq(chi.squared.statistic.gumbel, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.gumbel <- c(expected.values.gumbel[1])
      for (i in 1+seq_along(expected.values.gumbel)) {
        cumulative.expected.values.gumbel[i] <- cumulative.expected.values.gumbel[i-1] + expected.values.gumbel[i]
      }
      cumulative.expected.values.gumbel <- cumulative.expected.values.gumbel/sum(expected.values.gumbel)
      cumulative.expected.values.gumbel <- cumulative.expected.values.gumbel[!is.na(cumulative.expected.values.gumbel)]
      g.max.gumbel <- max(abs(cumulative.data - cumulative.expected.values.gumbel))
      if (g.max.gumbel < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.gumbel <- "Accept"
      } else {KS.gumbel <- "Reject"}
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
      print(paste("Fitting Frechet distribution..."))
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
      expected.values.frechet <- dist.frechet(res.hist$mids, dist.frechet.opt$par[1], dist.frechet.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.frechet <- sum((res.hist$counts - expected.values.frechet)^2 / expected.values.frechet)
      chi.squared.pvalue.frechet <- 1-pchisq(chi.squared.statistic.frechet, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.frechet <- c(expected.values.frechet[1])
      for (i in 1+seq_along(expected.values.frechet)) {
        cumulative.expected.values.frechet[i] <- cumulative.expected.values.frechet[i-1] + expected.values.frechet[i]
      }
      cumulative.expected.values.frechet <- cumulative.expected.values.frechet/sum(expected.values.frechet)
      cumulative.expected.values.frechet <- cumulative.expected.values.frechet[!is.na(cumulative.expected.values.frechet)]
      g.max.frechet <- max(abs(cumulative.data - cumulative.expected.values.frechet))
      if (g.max.frechet < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.frechet <- "Accept"
      } else {KS.frechet <- "Reject"}
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
      print(paste("Fitting Weibull distribution..."))
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
      expected.values.weibull <- dist.weibull(res.hist$mids, dist.weibull.opt$par[1], dist.weibull.opt$par[2])*length(data)*(res.hist$breaks[2] - res.hist$breaks[1])
      chi.squared.statistic.weibull <- sum((res.hist$counts - expected.values.weibull)^2 / expected.values.weibull)
      chi.squared.pvalue.weibull <- 1-pchisq(chi.squared.statistic.weibull, length(res.hist$counts)-3)
      # Kolmogorov-Smirnov
      cumulative.expected.values.weibull <- c(expected.values.weibull[1])
      for (i in 1+seq_along(expected.values.weibull)) {
        cumulative.expected.values.weibull[i] <- cumulative.expected.values.weibull[i-1] + expected.values.weibull[i]
      }
      cumulative.expected.values.weibull <- cumulative.expected.values.weibull/sum(expected.values.weibull)
      cumulative.expected.values.weibull <- cumulative.expected.values.weibull[!is.na(cumulative.expected.values.weibull)]
      g.max.weibull <- max(abs(cumulative.data - cumulative.expected.values.weibull))
      if (g.max.weibull < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
        KS.weibull <- "Accept"
      } else {KS.weibull <- "Reject"}
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
  } else {print("Please select TRUE/FALSE in extreme.values")}

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

  return(kernel.fit)
}
