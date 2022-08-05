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
#' @details This function fits one or more dispersal kernels with 1-2 parameters, by estimating the
#' distribution of kernel parameters (θ∈R) maximizing the likelihood function.
#' @return Returns a list with all maximum likelihood calculations for each selected distribution.
#' Furthermore, two summary tables may be accessed.
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

dispersal.kernel <- function (data, distribution = "all", order.by = "AICc", confidence.level = 0.95) {
  # require("msm")
  # require("numDeriv")
  if (confidence.level > 1 || confidence.level < 0)
  {stop('value for confidence.level must be between 0 and 1')}

  kernel.fit <- list()
  values <- setNames(data.frame(matrix(ncol = 21, nrow = 1)),  c("AIC", "AICc", "BIC", "Chi-squared value", "Chi-squared significance",
                                                                 "Kolmogorov-Smirnov value", "K-S significance",
                                                                 "Parameter 1", "Parameter 1 lower CI", "Parameter 1 upper CI",
                                                                 "Parameter 2", "Parameter 2 lower CI", "Parameter 2 upper CI",
                                                                 "Mean", "Mean SE",
                                                                 "Standard Deviation", "Standard Deviation SE",
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
  possible.distributions <- c("all", "rayleigh", "exponential", "general normal", "2Dt", "geometric", "lognormal", "wald", "weibull", "gamma")
  if(any(!(distribution %in% possible.distributions)))
  	stop("Possible values for 'distribution' are any of: ", paste(possible.distributions, collapse=", "))

    if ("all" %in% distribution | "rayleigh" %in% distribution) {
      message(paste("Fitting Rayleigh distribution..."))
      rayleigh.values <- rayleigh.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$rayleigh <- rayleigh.values$opt
      values <- rbind(values, setNames(rayleigh.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Rayleigh"
    }
    if ("all" %in% distribution | "exponential" %in% distribution) {
      message(paste("Fitting Exponential distribution..."))
      exponential.values <- exponential.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$exponential <- exponential.values$opt
      values <- rbind(values, setNames(exponential.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Exponential"
    }
    if ("all" %in% distribution | "general normal" %in% distribution) {
      message(paste("Fitting Generalized Normal distribution..."))
      generalnormal.values <- generalnormal.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$generalnormal <- generalnormal.values$opt
      values <- rbind(values, setNames(generalnormal.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Generalized Normal"
    }
    ### 4 ### 2Dt / BIVARIATE STUDENT'S T
    if ("all" %in% distribution | "2Dt" %in% distribution) {
      message(paste("Fitting 2Dt distribution..."))
      twodt.values <- twodt.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$twodt <- twodt.values$opt
      values <- rbind(values, setNames(twodt.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "2Dt"
    }
    ### 5 ### POWER-LAW / GEOMETRIC
    if ("all" %in% distribution | "geometric" %in% distribution) {
      message(paste("Fitting Geometric distribution..."))
      geometric.values <- geometric.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$geometric <- geometric.values$opt
      values <- rbind(values, setNames(geometric.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Geometric"
    }
    ### 6 ### INVERSE POWER-LAW - UNDEFINED

    ### 7 ### LOGISTIC ###
    # if ("all" %in% distribution | "logistic" %in% distribution) {
    #   message(paste("Fitting Logistic distribution..."))
    #   logistic.values <- logistic.function(kernel.fit$data, chi.res.hist, ks.res.hist)
    #
    #   kernel.fit$logistic <- logistic.values$opt
    #   values <- rbind(values, setNames(logistic.values$res, names(values)))
    #   row.names(values)[length(values$AIC)] <- "Logistic"
    # }
    ### 8 ### LOGNORMAL ### DONE!
    if ("all" %in% distribution | "lognormal" %in% distribution) {
      message(paste("Fitting Log-Normal distribution..."))
      lognorm.values <- lognorm.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$lognorm <- lognorm.values$opt
      values <- rbind(values, setNames(lognorm.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Log-Normal"
    }
    ### 9 ### GAUSSIAN MIXTURE

    ### 10 ### WALD (INVERSE GAUSSIAN) ### DONE!
    if ("all" %in% distribution | "wald" %in% distribution) {
      message(paste("Fitting Wald distribution..."))
      wald.values <- wald.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$wald <- wald.values$opt
      values <- rbind(values, setNames(wald.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Wald"
    }
    ### 11 ### WEIBULL ### DONE!
    if ("all" %in% distribution | "weibull" %in% distribution) {
      message(paste("Fitting Weibull distribution..."))
      weibull.values <- weibull.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$weibull <- weibull.values$opt
      values <- rbind(values, setNames(weibull.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Weibull"
    }
    ### 12 ### GAMMA ### DONE!
    if ("all" %in% distribution | "gamma" %in% distribution) {
      message(paste("Fitting Gamma distribution..."))
      gamma.values <- gamma.function(kernel.fit$data, chi.res.hist, ks.res.hist, confidence.level)

      kernel.fit$gamma <- gamma.values$opt
      values <- rbind(values, setNames(gamma.values$res, names(values)))
      row.names(values)[length(values$AIC)] <- "Gamma"
      # end
    }
    ### 13 ### LOG-SECH / HYPERBOLIC SECANT ### DONE!
    # if ("all" %in% distribution | "log-sech" %in% distribution) {
    #   message(paste("Fitting Log-sech distribution..."))
    #   logsech.values <- logsech.function(kernel.fit$data, chi.res.hist, ks.res.hist)
    #
    #   kernel.fit$logsech <- logsech.values$opt
    #   values <- rbind(values, setNames(logsech.values$res, names(values)))
    #   row.names(values)[length(values$AIC)] <- "Log-sech"
    #   # end
    # }
    ### 14 ### CAUCHY
    # if ("all" %in% distribution | "cauchy" %in% distribution) {
    #   message(paste("Fitting Cauchy distribution..."))
    #   cauchy.values <- cauchy.function(kernel.fit$data, chi.res.hist, ks.res.hist)
    #
    #   kernel.fit$cauchy <- cauchy.values$opt
    #   values <- rbind(values, setNames(cauchy.values$res, names(values)))
    #   row.names(values)[length(values$AIC)] <- "Cauchy"
    #   # end
    # }

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
  kernel.fit$values <- kernel.fit$values[,c(1, 22, 2, 23, 3, 24, 25, 4:14, 16, 18, 20)]

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
  kernel.fit$distribution.parameters <- kernel.fit$values[, 12:21]
  class(kernel.fit) <- "dispfit"
  kernel.fit
}


