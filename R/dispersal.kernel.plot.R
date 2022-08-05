#' Plots from dispersal kernel fit for dispersal data
#'
#' Plots the distributions previously fitted by \code{\link{dispersal.kernel}} against the data density plot.
#'
#' @param data Output object from the \code{\link{dispersal.kernel}} function.
#' @param fit.criteria Either a numeric value referring to the number of distributions to plot (ordered according to the order.by argument from dispersal.kernel); a character string choosing to plot “all” distributions, to plot the distributions with an “AIC”, “AICc”, or “BIC” difference to the top distribution bellow the number set by criteria.dif, or to plot specific distributions. By default, only the top distribution is plotted.
#' @param criteria.dif A numeric value used when “AIC”, “AICc”, or "BIC" are also selected in fit.criteria. Its value refers to the estimator difference to the top distribution. Only distributions with their estimator difference to the top model falling below this value will be plotted. The default is 2.
#' @param envelopes A boolean (TRUE/FALSE) value choosing whether to plot the confidence envelopes associated with each distribution. Default is TRUE.
#' @param plot.data A boolean (TRUE/FALSE) value selecting to plot the density of the data.
#' @seealso \code{\link{dispersal.kernel}}
#' @keywords kernel plot
#' @import ggplot2 reshape2 graphics
#' @export
#' @examples
#'
#' sim <- rexp(100, rate = 0.01)
#'
#' test <- dispersal.kernel(sim)
#'
#' plot(test)

plot.dispfit <- function (data, fit.criteria = NULL, criteria.dif = 2, envelopes = TRUE, plot.data = TRUE) {

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
  } else if (isFALSE(fit.criteria)){
    return(ggplot2::ggplot() +
             ggplot2::theme_light() +
             ggplot2::stat_density(data = data.frame(x=data$data), ggplot2::aes(x=x), colour = "black", geom = "line", size = 1) +
             ggplot2::labs(x = "Distance (m)", y = "Probability"))
  }

  pred <- predict.dispfit(data = data, fit.criteria = fit.criteria, criteria.dif = criteria.dif, envelopes = envelopes)
  pred.basic <- lapply(pred, function(x) x[,1:2])

  pred.basic$melt <- reshape2::melt(pred.basic, id = "distance")
  pred.basic$melt$variable <- as.character(pred.basic$melt$variable)
  pred.basic$melt$variable[pred.basic$melt$variable == "rayleigh"] <- "Rayleigh"
  pred.basic$melt$variable[pred.basic$melt$variable == "exponential"] <- "Exponential"
  pred.basic$melt$variable[pred.basic$melt$variable == "generalnormal"] <- "Generalized Normal"
  pred.basic$melt$variable[pred.basic$melt$variable == "twodt"] <- "2Dt"
  pred.basic$melt$variable[pred.basic$melt$variable == "geometric"] <- "Geometric"
  pred.basic$melt$variable[pred.basic$melt$variable == "logistic"] <- "Logistic"
  pred.basic$melt$variable[pred.basic$melt$variable == "lognormal"] <- "Log-Normal"
  pred.basic$melt$variable[pred.basic$melt$variable == "wald"] <- "Wald"
  pred.basic$melt$variable[pred.basic$melt$variable == "weibull"] <- "Weibull"
  pred.basic$melt$variable[pred.basic$melt$variable == "gamma"] <- "Gamma"
  pred.basic$melt$variable[pred.basic$melt$variable == "logsech"] <- "Log-sech"

  if (isTRUE(envelopes)) {
    if ("Rayleigh" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Rayleigh"] <- pred$rayleigh$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Rayleigh"] <- pred$rayleigh$lwr
      }
    if ("Exponential" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Exponential"] <- pred$exponential$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Exponential"] <- pred$exponential$lwr
    }
    if ("Generalized Normal" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Generalized Normal"] <- pred$generalnormal$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Generalized Normal"] <- pred$generalnormal$lwr
    }
    if ("2Dt" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "2Dt"] <- pred$twodt$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "2Dt"] <- pred$twodt$lwr
      }
    if ("Geometric" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Geometric"] <- pred$geometric$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Geometric"] <- pred$geometric$lwr
      }
    if ("Logistic" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Logistic"] <- pred$logistic$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Logistic"] <- pred$logistic$lwr
    }
    if ("Log-Normal" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Log-Normal"] <- pred$lognormal$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Log-Normal"] <- pred$lognormal$lwr
    }
    if ("Wald" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Wald"] <- pred$wald$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Wald"] <- pred$wald$lwr
    }
    if ("Weibull" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Weibull"] <- pred$weibull$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Weibull"] <- pred$weibull$lwr
    }
    if ("Gamma" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Gamma"] <- pred$gamma$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Gamma"] <- pred$gamma$lwr
    }
    if ("Log-sech" %in% best.fit) {
      pred.basic$melt$upr[pred.basic$melt$variable == "Log-sech"] <- pred$logsech$upr
      pred.basic$melt$lwr[pred.basic$melt$variable == "Log-sech"] <- pred$logsech$lwr
    }

    ggplot2::ggplot() +
      ggplot2::theme_light() +
      ggplot2::geom_line(data = pred.basic$melt, ggplot2::aes(x = distance, y = value, colour = variable), lineend = "round") +
      ggplot2::geom_line(data = na.omit(pred.basic$melt), ggplot2::aes(x = distance, y=upr, col = variable), linetype = "dotted", show.legend = F) +
      ggplot2::geom_line(data = na.omit(pred.basic$melt), ggplot2::aes(x = distance, y=lwr, col = variable), linetype = "dotted", show.legend = F) +
      ggplot2::geom_ribbon(data = na.omit(pred.basic$melt), ggplot2::aes(x = distance, ymax = upr, ymin = lwr, fill = variable), alpha=0.15, show.legend = F) +
      { if(plot.data) ggplot2::stat_density(data = data.frame(x=data$data), ggplot2::aes(x=x), colour = "black", geom = "line", size = 1) } +
      #ggplot2::ylim(0, (1.2*max(pred.basic$melt$upr, na.rm=T))) +
      ggplot2::coord_cartesian(ylim=c(0, max(hist(data$data, breaks = length(data$data), plot = F)$density))) +
      ggplot2::labs(x = "Distance", y = "Density", colour = "Distribution")
  }
  else {
    ggplot2::ggplot() +
      ggplot2::theme_light() +
      ggplot2::geom_line(data = pred.basic$melt, ggplot2::aes(x = distance, y = value, colour = variable), lineend = "round") +
      { if(plot.data) ggplot2::stat_density(data = data.frame(x=data$data), ggplot2::aes(x=x), colour = "black", geom = "line", size = 1)} +
      ggplot2::labs(x = "Distance", y = "Density", colour = "Distribution")
  }
}
