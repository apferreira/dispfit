#' Predicted data from dispersal kernel fit for dispersal data
#'
#' Plots the distributions previously fitted by \code{\link{dispersal.kernel}} against the data density plot.
#'
#' @param data Output object from the \code{\link{dispersal.kernel}} function.
#' @param fit.criteria Either a numeric value referring to the number of distributions to plot (ordered according to the order.by argument from dispersal.kernel); a character string choosing to plot “all” distributions, to plot the distributions with an “AIC”, “AICc”, or “BIC” difference to the top distribution bellow the number set by criteria.dif, or to plot specific distributions. By default, only the top distribution is plotted.
#' @param criteria.dif A numeric value used when “AIC”, “AICc”, or "BIC" are also selected in fit.criteria. Its value refers to the estimator difference to the top distribution. Only distributions with their estimator difference to the top model falling below this value will be plotted. The default is 2.
#' @param envelopes A boolean (TRUE/FALSE) value choosing whether to calculate the confidence envelopes associated with each distribution. Default is TRUE.
#' @param level A numeric value varying from 0 to 1, defining the confidence interval to be used when calculating confidence envelopes.
#' @seealso \code{\link{dispersal.kernel}}
#' @keywords kernel predict
#' @import ggplot2 reshape2 graphics
#' @export
#' @examples
#'
#' sim <- rexp(100, rate = 0.01)
#'
#' test <- dispersal.kernel(sim)
#'
#' predict(test)

predict.dispfit <- function(data, fit.criteria = NULL, criteria.dif = 2,
                            envelopes = TRUE,
                            # level = 0.95,
                            # se.fit = FALSE,
                            #  interval = c("none", "confidence", "prediction"),
                            #  type = c("response", "terms"),
                            #  terms = NULL, na.action = na.pass,
                            #  pred.var = res.var/weights, weights = 1,
                            ...)
{
  if (is.null(fit.criteria)) {
    best.fit <- row.names(data$values[1,])
  } else if ("AIC" %in% fit.criteria) {
    best.fit <- row.names(data$values[which(data$values$`Delta AIC` < criteria.dif),])
  } else if ("AICc" %in% fit.criteria) {
    best.fit <- row.names(data$values[which(data$values$`Delta AICc` < criteria.dif),])
  } else if ("all" %in% fit.criteria) {
    best.fit <- row.names(data$values)
  } else if (is.character(fit.criteria)) {
    best.fit <- fit.criteria
  } else if (is.numeric(fit.criteria)) {
    best.fit <- row.names(data$values[c(1:fit.criteria),])
  }

  n <- 100
  x <- seq(1, floor(max(data$data)), length.out = n)
  # n <- length(x)

  pred.disp <- list()

  if ("Rayleigh" %in% best.fit) {
    dist.rayleigh <- function (r, a) {
      fg <- 2*pi*r*(1/(pi*a^2)) * exp(-r^2/a^2)
      return(fg)
    }
    pred.disp$rayleigh <- data.frame(distance = x, rayleigh = dist.rayleigh(x, data$values["Rayleigh","Parameter 1"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["Rayleigh","Parameter 1 lower CI"]) && !is.na(data$values["Rayleigh","Parameter 1 upper CI"]) &&
         !data$values["Rayleigh","Parameter 1 upper CI"] == Inf) {
      seq.rayleigh <- seq(data$values["Rayleigh","Parameter 1 lower CI"], data$values["Rayleigh","Parameter 1 upper CI"], length.out = n)
      df.seq.rayleigh <- data.frame(1:n)
      for (i in 1:n) {
        df.seq.rayleigh[i] <- dist.rayleigh(x, seq.rayleigh[i])
      }
      pred.disp$rayleigh$lwr <- apply(df.seq.rayleigh, 1, FUN = min)
      pred.disp$rayleigh$upr <- apply(df.seq.rayleigh, 1, FUN = max)
      }
      else {
          warning("Parameter CI for Rayleigh are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$rayleigh$lwr <- NA
      pred.disp$rayleigh$upr <- NA
        }
    }
  }
  if ("Exponential" %in% best.fit) {
    dist.exponential <- function (r, a) {
      fexponential <-  2*pi*r*(1 / (2 * pi * a ^ 2 )) * exp(-r/a) # corrected function, adapted from Nathan 2012
      return(fexponential)
    }
    pred.disp$exponential <- data.frame(distance = x, exponential = dist.exponential(x, data$values["Exponential","Parameter 1"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["Exponential","Parameter 1 lower CI"]) && !is.na(data$values["Exponential","Parameter 1 upper CI"]) &&
         !data$values["Exponential","Parameter 1 upper CI"] == Inf) {
      seq.exponential <- seq(data$values["Exponential","Parameter 1 lower CI"], data$values["Exponential","Parameter 1 upper CI"], length.out = n)
      df.seq.exponential <- data.frame(1:n)
      for (i in 1:n) {
        df.seq.exponential[i] <- dist.exponential(x, seq.exponential[i])
      }
      pred.disp$exponential$lwr <- apply(df.seq.exponential, 1, FUN = min)
      pred.disp$exponential$upr <- apply(df.seq.exponential, 1, FUN = max)
      }
      else {
        warning("Parameter CI for Exponential are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$exponential$lwr <- NA
      pred.disp$exponential$upr <- NA
      }
    }
  }
  if ("Generalized Normal" %in% best.fit) {
    dist.generalnormal <- function (r, a, b) {
      fgeneralnormal <- 2*pi*r*(b / (2 * pi * (a^2) * gamma(2 / b))) * exp(-(r / a) ^ b)
      return(fgeneralnormal)
    }
    pred.disp$generalnormal <- data.frame(distance = x, generalnormal = dist.generalnormal(x, data$values["Generalized Normal","Parameter 1"], data$values["Generalized Normal","Parameter 2"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["Generalized Normal","Parameter 1 lower CI"]) && !is.na(data$values["Generalized Normal","Parameter 1 upper CI"]) &&
         !data$values["Generalized Normal","Parameter 1 upper CI"] == Inf &&
         !is.na(data$values["Generalized Normal","Parameter 2 lower CI"]) && !is.na(data$values["Generalized Normal","Parameter 2 upper CI"]) &&
         !data$values["Generalized Normal","Parameter 2 upper CI"] == Inf) {
      seq.generalnormal.par.1 <- seq(data$values["Generalized Normal","Parameter 1 lower CI"], data$values["Generalized Normal","Parameter 1 upper CI"], length.out = n)
      seq.generalnormal.par.2 <- seq(data$values["Generalized Normal","Parameter 2 lower CI"], data$values["Generalized Normal","Parameter 2 upper CI"], length.out = n)
      df.seq.generalnormal <- data.frame(1:n)

      list.seq.generalnormal <- list()
      c = 1
      for (i in seq_along(1:n)) {
        for (j in seq_along(1:n)) {
          list.seq.generalnormal[[c]] <- dist.generalnormal(x, seq.generalnormal.par.1[i], seq.generalnormal.par.2[j])
          c = c + 1
        }
      }
      df.seq.generalnormal <- data.frame()
      df.seq.generalnormal <- do.call(cbind, list.seq.generalnormal)

      pred.disp$generalnormal$lwr <- apply(df.seq.generalnormal, 1, FUN = min)
      pred.disp$generalnormal$upr <- apply(df.seq.generalnormal, 1, FUN = max)
      }
      else {
        warning("Parameter(s) CI for Generalized Normal are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$generalnormal$lwr <- NA
      pred.disp$generalnormal$upr <- NA
      }
    }
    }
  if ("2Dt" %in% best.fit) {
    dist.2dt <- function (r, a, b) {
      f2dt <- 2*pi*r*((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))
      return(f2dt)
    }
    pred.disp$twodt <- data.frame(distance = x, twodt = dist.2dt(x, data$values["2Dt","Parameter 1"], data$values["2Dt","Parameter 2"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["2Dt","Parameter 1 lower CI"]) && !is.na(data$values["2Dt","Parameter 1 upper CI"]) &&
         !data$values["2Dt","Parameter 1 upper CI"] == Inf &&
         !is.na(data$values["2Dt","Parameter 2 lower CI"]) && !is.na(data$values["2Dt","Parameter 2 upper CI"]) &&
         !data$values["2Dt","Parameter 2 upper CI"] == Inf) {
      seq.2dt.par.1 <- seq(data$values["2Dt","Parameter 1 lower CI"], data$values["2Dt","Parameter 1 upper CI"], length.out = n)
      seq.2dt.par.2 <- seq(data$values["2Dt","Parameter 2 lower CI"], data$values["2Dt","Parameter 2 upper CI"], length.out = n)

      list.seq.2dt <- list()
      c = 1
      for (i in seq_along(1:n)) {
        for (j in seq_along(1:n)) {
          list.seq.2dt[[c]] <- dist.2dt(x, seq.2dt.par.1[i], seq.2dt.par.2[j])
          c = c + 1
        }
      }
      df.seq.2dt <- data.frame()
      df.seq.2dt <- do.call(cbind, list.seq.2dt)

      pred.disp$twodt$lwr <- apply(df.seq.2dt, 1, FUN = min)
      pred.disp$twodt$upr <- apply(df.seq.2dt, 1, FUN = max)
      }
      else {
        warning("Parameter(s) CI for 2Dt are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$twodt$lwr <- NA
      pred.disp$twodt$upr <- NA
      }
    }
  }
  if ("Geometric" %in% best.fit) {
    dist.geometric <- function (r, a, b) {
      fgeometric <- 2*pi*r*(((b - 2) * (b - 1)) / (2 * pi * (a^2))) * ((1 + (r / a)) ^ -b)
      return(fgeometric)
    }
    pred.disp$geometric <- data.frame(distance = x, geometric = dist.geometric(x, data$values["Geometric","Parameter 1"], data$values["Geometric","Parameter 2"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["Geometric","Parameter 1 lower CI"]) && !is.na(data$values["Geometric","Parameter 1 upper CI"]) &&
         !data$values["Geometric","Parameter 1 upper CI"] == Inf &&
         !is.na(data$values["Geometric","Parameter 2 lower CI"]) && !is.na(data$values["Geometric","Parameter 2 upper CI"]) &&
         !data$values["Geometric","Parameter 2 upper CI"] == Inf) {
      seq.geometric.par.1 <- seq(data$values["Geometric","Parameter 1 lower CI"], data$values["Geometric","Parameter 1 upper CI"], length.out = n)
      seq.geometric.par.2 <- seq(data$values["Geometric","Parameter 2 lower CI"], data$values["Geometric","Parameter 2 upper CI"], length.out = n)

      list.seq.geometric <- list()
      c = 1
      for (i in seq_along(1:n)) {
        for (j in seq_along(1:n)) {
          list.seq.geometric[[c]] <- dist.geometric(x, seq.geometric.par.1[i], seq.geometric.par.2[j])
          c = c + 1
        }
      }
      df.seq.geometric <- data.frame()
      df.seq.geometric <- do.call(cbind, list.seq.geometric)

      pred.disp$geometric$lwr <- apply(df.seq.geometric, 1, FUN = min)
      pred.disp$geometric$upr <- apply(df.seq.geometric, 1, FUN = max)
      }
      else {
        warning("Parameter(s) CI for Geometric are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$geometric$lwr <- NA
      pred.disp$geometric$upr <- NA
      }
    }
  }
  if ("Log-Normal" %in% best.fit) {
    dist.lognormal <- function (r, a, b) {
      flognorm <- 2*pi*r * (1 / (((2 * pi) ^ (3/2)) * (b * (r ^ 2)))) * exp(-(log(r / a)^2) / (2 * (b ^ 2)))
      return(flognorm)
    }
    pred.disp$lognormal <- data.frame(distance = x, lognormal = dist.lognormal(x, data$values["Log-Normal","Parameter 1"], data$values["Log-Normal","Parameter 2"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["Log-Normal","Parameter 1 lower CI"]) && !is.na(data$values["Log-Normal","Parameter 1 upper CI"]) &&
         !data$values["Log-Normal","Parameter 1 upper CI"] == Inf &&
         !is.na(data$values["Log-Normal","Parameter 2 lower CI"]) && !is.na(data$values["Log-Normal","Parameter 2 upper CI"]) &&
         !data$values["Log-Normal","Parameter 2 upper CI"] == Inf) {
      seq.lognormal.par.1 <- seq(data$values["Log-Normal","Parameter 1 lower CI"], data$values["Log-Normal","Parameter 1 upper CI"], length.out = n)
      seq.lognormal.par.2 <- seq(data$values["Log-Normal","Parameter 2 lower CI"], data$values["Log-Normal","Parameter 2 upper CI"], length.out = n)
      list.seq.lognorma <- list()
      c = 1
      for (i in seq_along(1:n)) {
        for (j in seq_along(1:n)) {
          list.seq.lognorma[[c]] <- dist.lognormal(x, seq.lognormal.par.1[i], seq.lognormal.par.2[j])
          c = c + 1
        }
      }
      df.seq.lognormal <- data.frame()
      df.seq.lognormal <- do.call(cbind, list.seq.lognorma)

      pred.disp$lognormal$lwr <- apply(df.seq.lognormal, 1, FUN = min)
      pred.disp$lognormal$upr <- apply(df.seq.lognormal, 1, FUN = max)
      }
      else {
        warning("Parameter(s) CI for Log-Normal are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$lognormal$lwr <- NA
      pred.disp$lognormal$upr <- NA
      }
    }
  }
  if ("Wald" %in% best.fit) {
    dist.wald <- function (r, a, b) {
      fwald <- 2*pi*r * (sqrt(b)/sqrt(8 * (pi^3) * (r^5))) * exp(-(b * ((r - a)^2))/(2 * (a^2) * r))
      return(fwald)
    }
    pred.disp$wald <- data.frame(distance = x, wald = dist.wald(x, data$values["Wald","Parameter 1"], data$values["Wald","Parameter 2"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["Wald","Parameter 1 lower CI"]) && !is.na(data$values["Wald","Parameter 1 upper CI"]) &&
         !data$values["Wald","Parameter 1 upper CI"] == Inf &&
         !is.na(data$values["Wald","Parameter 2 lower CI"]) && !is.na(data$values["Wald","Parameter 2 upper CI"]) &&
         !data$values["Wald","Parameter 2 upper CI"] == Inf) {
      seq.wald.par.1 <- seq(data$values["Wald","Parameter 1 lower CI"], data$values["Wald","Parameter 1 upper CI"], length.out = n)
      seq.wald.par.2 <- seq(data$values["Wald","Parameter 2 lower CI"], data$values["Wald","Parameter 2 upper CI"], length.out = n)
      list.seq.wald <- list()
      c = 1
      for (i in seq_along(1:n)) {
        for (j in seq_along(1:n)) {
          list.seq.wald[[c]] <- dist.wald(x, seq.wald.par.1[i], seq.wald.par.2[j])
          c = c + 1
        }
      }
      df.seq.wald <- data.frame()
      df.seq.wald <- do.call(cbind, list.seq.wald)

      pred.disp$wald$lwr <- apply(df.seq.wald, 1, FUN = min)
      pred.disp$wald$upr <- apply(df.seq.wald, 1, FUN = max)
      }
      else {
        warning("Parameter(s) CI for Wald are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$wald$lwr <- NA
      pred.disp$wald$upr <- NA
      }
    }
  }
  if ("Weibull" %in% best.fit) {
    dist.weibull <- function (r, a, b) {
      fw <- 2*pi*r * (b/(2*pi*a^b)) * (r^(b-2)) * exp(-(r^b/a^b)) ## function from Austerlitz 2004
      return(fw)
    }
    pred.disp$weibull <- data.frame(distance = x, weibull = dist.weibull(x, data$values["Weibull","Parameter 1"], data$values["Weibull","Parameter 2"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["Weibull","Parameter 1 lower CI"]) && !is.na(data$values["Weibull","Parameter 1 upper CI"]) &&
         !data$values["Weibull","Parameter 1 upper CI"] == Inf &&
         !is.na(data$values["Weibull","Parameter 2 lower CI"]) && !is.na(data$values["Weibull","Parameter 2 upper CI"]) &&
         !data$values["Weibull","Parameter 2 upper CI"] == Inf) {
      seq.weibull.par.1 <- seq(data$values["Weibull","Parameter 1 lower CI"], data$values["Weibull","Parameter 1 upper CI"], length.out = n)
      seq.weibull.par.2 <- seq(data$values["Weibull","Parameter 2 lower CI"], data$values["Weibull","Parameter 2 upper CI"], length.out = n)
      list.seq.weibull <- list()
      c = 1
      for (i in seq_along(1:n)) {
        for (j in seq_along(1:n)) {
          list.seq.weibull[[c]] <- dist.weibull(x, seq.weibull.par.1[i], seq.weibull.par.2[j])
          c = c + 1
        }
      }
      df.seq.weibull <- data.frame()
      df.seq.weibull <- do.call(cbind, list.seq.weibull)

      pred.disp$weibull$lwr <- apply(df.seq.weibull, 1, FUN = min)
      pred.disp$weibull$upr <- apply(df.seq.weibull, 1, FUN = max)
      }
      else {
        warning("Parameter(s) CI for Weibull are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$weibull$lwr <- NA
      pred.disp$weibull$upr <- NA
      }
    }
  }
  if ("Gamma" %in% best.fit) {
    dist.gamma <- function (r, a, b) {
      fgamma <- 2*pi*r * (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
      return(fgamma)
    }
    pred.disp$gamma <- data.frame(distance = x, gamma = dist.gamma(x, data$values["Gamma","Parameter 1"], data$values["Gamma","Parameter 2"]))
    if (isTRUE(envelopes)) {
      if(!is.na(data$values["Gamma","Parameter 1 lower CI"]) && !is.na(data$values["Gamma","Parameter 1 upper CI"]) &&
         !data$values["Gamma","Parameter 1 upper CI"] == Inf &&
         !is.na(data$values["Gamma","Parameter 2 lower CI"]) && !is.na(data$values["Gamma","Parameter 2 upper CI"]) &&
         !data$values["Gamma","Parameter 2 upper CI"] == Inf) {
      seq.gamma.par.1 <- seq(data$values["Gamma","Parameter 1 lower CI"], data$values["Gamma","Parameter 1 upper CI"], length.out = n)
      seq.gamma.par.2 <- seq(data$values["Gamma","Parameter 2 lower CI"], data$values["Gamma","Parameter 2 upper CI"], length.out = n)
      list.seq.gamma <- list()
      c = 1
      for (i in seq_along(1:n)) {
        for (j in seq_along(1:n)) {
          list.seq.gamma[[c]] <- dist.gamma(x, seq.gamma.par.1[i], seq.gamma.par.2[j])
          c = c + 1
        }
      }
      df.seq.gamma <- data.frame()
      df.seq.gamma <- do.call(cbind, list.seq.gamma)

      pred.disp$gamma$lwr <- apply(df.seq.gamma, 1, FUN = min)
      pred.disp$gamma$upr <- apply(df.seq.gamma, 1, FUN = max)
      }
      else {
        warning("Parameter(s) CI for Gamma are NA or infinite. Unable to estimate confidence envelopes.")
      pred.disp$gamma$lwr <- NA
      pred.disp$gamma$upr <- NA
      }
    }
  }
    pred.disp
}
