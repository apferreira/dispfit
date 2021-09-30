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
                            envelopes = TRUE, level = 0.95,
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
  upr <- 1-(1-level)/2
  lwr <- (1-level)/2

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

  ## Calculate Confidence Envelopes
  n <- 4000
  if ("Rayleigh" %in% best.fit) {
    a <- msm::rtnorm(n, data$values["Rayleigh", "Parameter 1"], data$values["Rayleigh", "Parameter 1 SE"], 0, Inf)
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
    all.sim$rayleigh$upr <- dist[upr*m,]
    all.sim$rayleigh$lwr <- dist[lwr*m,]
  }
  if ("Exponential" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Exponential", "Parameter 1"],data$values["Exponential", "Parameter 1 SE"], 0, Inf)
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
    all.sim$exponential$upr <- dist[upr*m,]
    all.sim$exponential$lwr <- dist[lwr*m,]
  }
  if ("Generalized Normal" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Generalized Normal", "Parameter 1"],data$values["Generalized Normal", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["Generalized Normal", "Parameter 2"],data$values["Generalized Normal", "Parameter 2 SE"], 0, Inf)
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
    all.sim$generalnormal$upr <- dist[upr*m,]
    all.sim$generalnormal$lwr <- dist[lwr*m,]
  }
  if ("2Dt" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["2Dt", "Parameter 1"],data$values["2Dt", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["2Dt", "Parameter 2"],data$values["2Dt", "Parameter 2 SE"], 0, Inf)
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
    all.sim$twodt$upr <- dist[upr*m,]
    all.sim$twodt$lwr <- dist[lwr*m,]
  }
  if ("Geometric" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Geometric", "Parameter 1"],data$values["Geometric", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["Geometric", "Parameter 2"],data$values["Geometric", "Parameter 2 SE"], 0, Inf)
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
    all.sim$geometric$upr <- dist[upr*m,]
    all.sim$geometric$lwr <- dist[lwr*m,]
  }
  if ("Logistic" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Logistic", "Parameter 1"],data$values["Logistic", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["Logistic", "Parameter 2"],data$values["Logistic", "Parameter 2 SE"], 0, Inf)
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
    all.sim$logistic$upr <- dist[upr*m,]
    all.sim$logistic$lwr <- dist[lwr*m,]
  }
  if ("Log-Normal" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Log-Normal", "Parameter 1"],data$values["Log-Normal", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["Log-Normal", "Parameter 2"],data$values["Log-Normal", "Parameter 2 SE"], 0, Inf)
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
    all.sim$lognormal$upr <- dist[upr*m,]
    all.sim$lognormal$lwr <- dist[lwr*m,]
  }
  if ("Wald" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Wald", "Parameter 1"],data$values["Wald", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["Wald", "Parameter 2"],data$values["Wald", "Parameter 2 SE"], 0, Inf)
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
    all.sim$wald$upr <- dist[upr*m,]
    all.sim$wald$lwr <- dist[lwr*m,]
  }
  if ("Weibull" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Weibull", "Parameter 1"],data$values["Weibull", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["Weibull", "Parameter 2"],data$values["Weibull", "Parameter 2 SE"], 0, Inf)
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
    all.sim$weibull$upr <- dist[upr*m,]
    all.sim$weibull$lwr <- dist[lwr*m,]
  }
  if ("Gamma" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Gamma", "Parameter 1"],data$values["Gamma", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["Gamma", "Parameter 2"],data$values["Gamma", "Parameter 2 SE"], 0, Inf)
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
    all.sim$gamma$upr <- dist[upr*m,]
    all.sim$gamma$lwr <- dist[lwr*m,]
  }
  if ("Log-sech" %in% best.fit) {
    a <- msm::rtnorm(n,data$values["Log-sech", "Parameter 1"],data$values["Log-sech", "Parameter 1 SE"], 0, Inf)
    b <- msm::rtnorm(n,data$values["Log-sech", "Parameter 2"],data$values["Log-sech", "Parameter 2 SE"], 0, Inf)
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
    all.sim$logsech$upr <- dist[upr*m,]
    all.sim$logsech$lwr <- dist[lwr*m,]
  }
  # all.sim$melt$low95[all.sim$melt$low95<0] <- 0
  all.sim
}
