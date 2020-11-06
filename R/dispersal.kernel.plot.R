#' Plots from dispersal kernel fit for dispersal data
#'
#' Plots the distributions previously fitted by dispersal.kernel against the data density plot.
#'
#' @param data Output object from the {\link{dispersal.kernel}} function.
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
  # require("ggplot2")
  # require("reshape2")

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

  all.sim <- list()
  x <- 1:floor(max(data$data))
  if ("Rayleigh" %in% best.fit) {
    dist.rayleigh <- function (r, a) {
      fg <- 2*pi*r*(1/(pi*a^2)) * exp(-r^2/a^2)
      return(fg)
    }
    all.sim$rayleigh <- data.frame(distance = x, rayleigh = dist.rayleigh(x, data$rayleigh$par))
  }
  if ("Exponential" %in% best.fit) {
    dist.exponential <- function (r, a) {
      fexponential <-  2*pi*r*(1 / (2 * pi * a ^ 2 )) * exp(-r/a) # corrected function, adapted from Nathan 2012
      return(fexponential)
    }
    all.sim$exponential <- data.frame(distance = x, exponential = dist.exponential(x, data$exponential$par))
  }
  if ("Generalized Normal" %in% best.fit) {
    dist.generalnormal <- function (r, a, b) {
      fgeneralnormal <- 2*pi*r*(b / (2 * pi * (a^2) * gamma(2 / b))) * exp(-(r / a) ^ b)
      return(fgeneralnormal)
    }
    all.sim$generalnormal <- data.frame(distance = x, generalnormal = dist.generalnormal(x, data$generalnormal$par[1], data$generalnormal$par[2]))
  }
  if ("2Dt" %in% best.fit) {
    dist.2dt <- function (r, a, b) {
      f2dt <- 2*pi*r*((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))
      return(f2dt)
    }
    all.sim$twodt <- data.frame(distance = x, twodt = dist.2dt(x, data$values["2Dt", "Parameter 1"], data$values["2Dt", "Parameter 2"]))
  }
  if ("Geometric" %in% best.fit) {
    dist.geometric <- function (r, a, b) {
      fgeometric <- 2*pi*r*(((b - 2) * (b - 1)) / (2 * pi * (a^2))) * ((1 + (r / a)) ^ -b)
      return(fgeometric)
    }
    all.sim$geometric <- data.frame(distance = x, geometric = dist.geometric(x, data$geometric$par[1], data$geometric$par[2]))
  }
  if ("Logistic" %in% best.fit) {
    dist.logistic <- function (r, a, b) {
      flogistic <- 2*pi*r*(b / (2 * pi * (a^2) * gamma(2/b) * gamma(1-(2/b)) )) * ((1 + ((r^b) / (a^b)))^(-1))
      return(flogistic)
    }
    all.sim$logistic <- data.frame(distance = x, logistic = dist.logistic(x, data$logistic$par[1], data$logistic$par[2]))
  }
  if ("Log-Normal" %in% best.fit) {
    dist.lognormal <- function (r, a, b) {
      flognorm <- 2*pi*r * (1 / (((2 * pi) ^ (3/2)) * (b * (r ^ 2)))) * exp(-(log(r / a)^2) / (2 * (b ^ 2)))
      return(flognorm)
    }
    all.sim$lognormal <- data.frame(distance = x, lognormal =  dist.lognormal(x, data$lognorm$par[1], data$lognorm$par[2]))
  }
  if ("Wald" %in% best.fit) {
    dist.wald <- function (r, a, b) {
      fwald <- 2*pi*r * (sqrt(b)/sqrt(8 * (pi^3) * (r^5))) * exp(-(b * ((r - a)^2))/(2 * (a^2) * r))
      return(fwald)
    }
    all.sim$wald <- data.frame(distance = x, wald = dist.wald(x, data$wald$par[1], data$wald$par[2]))
  }
  if ("Weibull" %in% best.fit) {
    dist.weibull <- function (r, a, b) {
      fw <- 2*pi*r * (b/(2*pi*a^b)) * (r^(b-2)) * exp(-(r^b/a^b)) ## function from Austerlitz 2004
      return(fw)
    }
    all.sim$weibull <- data.frame(distance = x, weibull = dist.weibull(x, data$weibull$par[1], data$weibull$par[2]))
  }
  if ("Gamma" %in% best.fit) {
    dist.gamma <- function (r, a, b) {
      fgamma <- 2*pi*r * (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
      return(fgamma)
    }
    all.sim$gamma <- data.frame(distance = x, gamma = dist.gamma(x, data$gamma$par[1], data$gamma$par[2]))
  }
  if ("Log-sech" %in% best.fit) {
    dist.logsech <- function (r, a, b) {
      flogsech <- 2*pi*r * (1 / ((pi^2) * b * (r^2))) / (((r / a)^(1 / b)) + ((r / a) ^ -(1 / b)))
      return(flogsech)
    }
    all.sim$logsech <- data.frame(distance = x, logsech = dist.logsech(x, data$logsech$par[1], data$logsech$par[2]))
  }
  # if ("Cauchy" %in% best.fit) {
  #   dist.cauchy <- function (r, a, b) {
  #     fcauchy <- 2*pi*r * (1/(2*pi*r))*(1 / (pi*b)) / (1 + ((r-a)/b)^2)
  #     return(fcauchy)
  #   }
  #   all.sim$cauchy <- data.frame(distance = 0:floor(max(data$data)), cauchy = dist.cauchy(0:floor(max(data$data)), data$cauchy$par[1], data$cauchy$par[2]))
  # }
  # if ("Gumbel" %in% best.fit) {
  #   dist.gumbel <- function (r, a, b) {
  #     fgumbel <- 2*pi*r * (1/(2*pi*r))*(1/b) * exp(-((r - a)/b + exp(-(r - a)/b))) ## pdf from Forbes 2011
  #     return(fgumbel)
  #   }
  #   all.sim$gumbel <- data.frame(distance = 0:floor(max(data$data)), gumbel = dist.gumbel(0:floor(max(data$data)), data$gumbel$par[1], data$gumbel$par[2]))
  # }
  # if ("Frechet" %in% best.fit) {
  #   dist.frechet <- function (r, a, b) {
  #     ffrechet <- 2*pi*r * (1/(2*pi*r)) * (a/b) * ((r/b)^(-1-a)) * exp(-(r/b)^-a)
  #     return(ffrechet)
  #   }
  #   all.sim$frechet <- data.frame(distance = 0:floor(max(data$data)), frechet = dist.frechet(0:floor(max(data$data)), data$frechet$par[1], data$frechet$par[2]))
  # }

  all.sim$melt <- reshape2::melt(all.sim, id = "distance")
  all.sim$melt$variable <- as.character(all.sim$melt$variable)
  all.sim$melt$variable[which(all.sim$melt$variable == "rayleigh")] <- "Rayleigh"
  all.sim$melt$variable[all.sim$melt$variable == "exponential"] <- "Exponential"
  all.sim$melt$variable[all.sim$melt$variable == "generalnormal"] <- "Generalized Normal"
  all.sim$melt$variable[all.sim$melt$variable == "twodt"] <- "2Dt"
  all.sim$melt$variable[all.sim$melt$variable == "geometric"] <- "Geometric"
  all.sim$melt$variable[all.sim$melt$variable == "logistic"] <- "Logistic"
  all.sim$melt$variable[all.sim$melt$variable == "lognormal"] <- "Log-Normal"
  all.sim$melt$variable[all.sim$melt$variable == "wald"] <- "Wald"
  all.sim$melt$variable[all.sim$melt$variable == "weibull"] <- "Weibull"
  all.sim$melt$variable[all.sim$melt$variable == "gamma"] <- "Gamma"
  all.sim$melt$variable[all.sim$melt$variable == "logsech"] <- "Log-sech"
  all.sim$melt$variable[all.sim$melt$variable == "cauchy"] <- "Cauchy"
  all.sim$melt$variable[all.sim$melt$variable == "gumbel"] <- "Gumbel"
  all.sim$melt$variable[all.sim$melt$variable == "frechet"] <- "Frechet"

  if (isTRUE(envelopes)) {
    n <- 4000
    if ("Rayleigh" %in% best.fit) {
      a <- rnorm(n,data$values["Rayleigh", "Parameter 1"],data$values["Rayleigh", "Parameter 1 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        y <- dist.rayleigh(x, as)
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Rayleigh Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Rayleigh"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Rayleigh"] <- dist[.025*m,]
    }
    if ("Exponential" %in% best.fit) {
      a <- rnorm(n,data$values["Exponential", "Parameter 1"],data$values["Exponential", "Parameter 1 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        y <- dist.exponential(x, as)
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Exponential Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Exponential"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Exponential"] <- dist[.025*m,]
    }
    if ("Generalized Normal" %in% best.fit) {
      a <- rnorm(n,data$values["Generalized Normal", "Parameter 1"],data$values["Generalized Normal", "Parameter 1 SE"])
      b <- rnorm(n,data$values["Generalized Normal", "Parameter 2"],data$values["Generalized Normal", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.generalnormal(x, as, bs)
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Generalized Normal Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Generalized Normal"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Generalized Normal"] <- dist[.025*m,]
    }
    if ("2Dt" %in% best.fit) {
      a <- rnorm(n,data$values["2Dt", "Parameter 1"],data$values["2Dt", "Parameter 1 SE"])
      b <- rnorm(n,data$values["2Dt", "Parameter 2"],data$values["2Dt", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.2dt(x, as, bs)
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the 2Dt Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "2Dt"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "2Dt"] <- dist[.025*m,]
    }
    if ("Geometric" %in% best.fit) {
      a <- rnorm(n,data$values["Geometric", "Parameter 1"],data$values["Geometric", "Parameter 1 SE"])
      b <- rnorm(n,data$values["Geometric", "Parameter 2"],data$values["Geometric", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.geometric(x, as, bs)
        dist[i,] <- y
      }
      dist <- apply(dist, 2, sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m)} # doesn't seem replicable, have to find another solution
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Geometric Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Geometric"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Geometric"] <- dist[.025*m,]
    }
    if ("Logistic" %in% best.fit) {
      a <- rnorm(n,data$values["Logistic", "Parameter 1"],data$values["Logistic", "Parameter 1 SE"])
      b <- rnorm(n,data$values["Logistic", "Parameter 2"],data$values["Logistic", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.logistic(x, as, bs)
        dist[i,] <- y
      }
      dist <- apply(dist, 2, sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Logistic Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Logistic"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Logistic"] <- dist[.025*m,]
    }
    if ("Log-Normal" %in% best.fit) {
      a <- rnorm(n,data$values["Log-Normal", "Parameter 1"],data$values["Log-Normal", "Parameter 1 SE"])
      b <- rnorm(n,data$values["Log-Normal", "Parameter 2"],data$values["Log-Normal", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.lognormal(x, as, bs)
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Log-Normal Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Log-Normal"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Log-Normal"] <- dist[.025*m,]
    }
    if ("Wald" %in% best.fit) {
      a <- rnorm(n,data$values["Wald", "Parameter 1"],data$values["Wald", "Parameter 1 SE"])
      b <- rnorm(n,data$values["Wald", "Parameter 2"],data$values["Wald", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.wald(x, as, bs)
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Wald Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Wald"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Wald"] <- dist[.025*m,]
    }
    if ("Weibull" %in% best.fit) {
      a <- rnorm(n,data$values["Weibull", "Parameter 1"],data$values["Weibull", "Parameter 1 SE"])
      b <- rnorm(n,data$values["Weibull", "Parameter 2"],data$values["Weibull", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.weibull(x, as, bs)
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Weibull Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Weibull"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Weibull"] <- dist[.025*m,]
    }
    if ("Gamma" %in% best.fit) {
      a <- rnorm(n,data$values["Gamma", "Parameter 1"],data$values["Gamma", "Parameter 1 SE"])
      b <- rnorm(n,data$values["Gamma", "Parameter 2"],data$values["Gamma", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.gamma(x, as, bs)
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Gamma Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Gamma"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Gamma"] <- dist[.025*m,]
    }
    if ("Log-sech" %in% best.fit) {
      a <- rnorm(n,data$values["Log-sech", "Parameter 1"],data$values["Log-sech", "Parameter 1 SE"])
      b <- rnorm(n,data$values["Log-sech", "Parameter 2"],data$values["Log-sech", "Parameter 2 SE"])
      dist <- matrix(0,nrow=n,ncol=length(x))
      for(i in 1:n){
        as <- a[i]
        bs <- b[i]
        y <- dist.logsech(x, as, bs)
        # y[is.nan(y)] <- NA
        dist[i,] <- y
      }
      dist <- apply(dist,2,sort)
      if(inherits(dist, "list")) {
        m <- min(sapply(dist, length))
        dist <- sapply(dist, '[', 1:m) }
      m <- length(dist)/length(x)
      if(m < n) message(n-m, ' random distributions were excluded from the Log-sech Distribution confidence envelopes')
      all.sim$melt$upp95[all.sim$melt$variable == "Log-sech"] <- dist[.975*m,]
      all.sim$melt$low95[all.sim$melt$variable == "Log-sech"] <- dist[.025*m,]
    }
    all.sim$melt$low95[all.sim$melt$low95<0] <- 0

    ggplot2::ggplot() +
      ggplot2::theme_light() +
      ggplot2::geom_line(data = all.sim$melt, ggplot2::aes(x = distance, y = value, colour = variable), lineend = "round") +
      ggplot2::geom_line(data = all.sim$melt, ggplot2::aes(x = distance, y=upp95, col = variable), linetype = "dotted", show.legend = F) +
      ggplot2::geom_line(data = all.sim$melt, ggplot2::aes(x = distance, y=low95, col = variable), linetype = "dotted", show.legend = F) +
      ggplot2::geom_ribbon(data = all.sim$melt, ggplot2::aes(x = distance, ymax = upp95, ymin = low95, fill = variable), alpha=0.15, show.legend = F) +
      { if(plot.data) ggplot2::stat_density(data = data.frame(x=data$data), ggplot2::aes(x=x), colour = "black", geom = "line", size = 1) } +
      ggplot2::ylim(0, max(all.sim$melt$upp95, na.rm=T)) +
      ggplot2::labs(x = "Distance (m)", y = "Probability", colour = "Distribution")
  }
  else {
   ggplot2::ggplot() +
      ggplot2::theme_light() +
      ggplot2::geom_line(data = all.sim$melt, ggplot2::aes(x = distance, y = value, colour = variable), lineend = "round") +
      { if(plot.data) ggplot2::stat_density(data = data.frame(x=data$data), ggplot2::aes(x=x), colour = "black", geom = "line", size = 1)} +
      ggplot2::labs(x = "Distance (m)", y = "Probability", colour = "Distribution")
  }
}
