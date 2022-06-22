generalnormal.function <- function (data, chi.res.hist, ks.res.hist, confidence.level) {
  log.dist.generalnormal <- function (r, par) {
    a <- par[1] ## scale
    b <- par[2] ## shape
    fgeneralnormal <- (b / (2 * pi * (a^2) * gamma(2 / b))) * exp(-(r^b / a^b))
    if(all(is.nan(fgeneralnormal)) || all(fgeneralnormal < 0))
      return(NaN)
    else
      return(-sum(log(fgeneralnormal)))
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

  #   # output values
  # AIC
  aic.generalnormal <- 2*length(dist.generalnormal.opt$par) + 2 * dist.generalnormal.opt$value ## AIC
  # AICc
  aicc.generalnormal <- aic.generalnormal + (2 * length(dist.generalnormal.opt$par)^2 + 2 * length(dist.generalnormal.opt$par))/(length(data) - length(dist.generalnormal.opt$par) - 1 )
  # BIC
  bic.generalnormal <-  2 * dist.generalnormal.opt$value + length(dist.generalnormal.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.generalnormal <- dist.generalnormal(chi.res.hist$mids, dist.generalnormal.opt$par[1], dist.generalnormal.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.generalnormal <- sum((chi.res.hist$counts - chi.expected.values.generalnormal)^2 / chi.expected.values.generalnormal)
  chi.squared.pvalue.generalnormal <- 1-pchisq(chi.squared.statistic.generalnormal, length(chi.res.hist$counts)-3)
  # Kolmogorov-Smirnov
  ks.expected.values.generalnormal <- dist.generalnormal(ks.res.hist$mids, dist.generalnormal.opt$par[1], dist.generalnormal.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
  simul.generalnormal <- c()
  for (i in seq_along(ks.res.hist$mids)) {
    simul.generalnormal <- c(simul.generalnormal, rep(ks.res.hist$mids[i], round(ks.expected.values.generalnormal[i], 0)))
  }
  ks.generalnormal <- ks.test(data, simul.generalnormal)
  g.max.generalnormal <- as.numeric(ks.generalnormal$statistic)
  KS.generalnormal <- as.numeric(ks.generalnormal$p.value)

  # cumulative.expected.values.generalnormal <- c(expected.values.generalnormal[1])
  # for (i in 1+seq_along(expected.values.generalnormal)) {
  #   cumulative.expected.values.generalnormal[i] <- cumulative.expected.values.generalnormal[i-1] + expected.values.generalnormal[i]
  # }
  # cumulative.expected.values.generalnormal <- cumulative.expected.values.generalnormal/sum(expected.values.generalnormal)
  # cumulative.expected.values.generalnormal <- cumulative.expected.values.generalnormal[!is.na(cumulative.expected.values.generalnormal)]
  # g.max.generalnormal <- max(abs(cumulative.data - cumulative.expected.values.generalnormal))
  # if (g.max.generalnormal < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
  #   KS.generalnormal <- "Accept"
  # } else {KS.generalnormal <- "Reject"}
  # mean dispersal distance
  # mean.generalnormal <- dist.generalnormal.opt$par[1] * (gamma(3/dist.generalnormal.opt$par[2])/gamma(2/dist.generalnormal.opt$par[2]))
  # parameter estimate
  par.1.generalnormal <- dist.generalnormal.opt$par[1]
  par.2.generalnormal <- dist.generalnormal.opt$par[2]
  # parameter estimate standard error
  par.1.se.generalnormal <- sqrt(diag(solve(numDeriv::hessian(log.dist.generalnormal, x=dist.generalnormal.opt$par, r=data))))[1]
  par.2.se.generalnormal <- sqrt(diag(solve(numDeriv::hessian(log.dist.generalnormal, x=dist.generalnormal.opt$par, r=data))))[2]
  # parameter estimate confidence intervals
  log.dist.generalnormal.ci <- function (r, a, b) {
    # a <- par[1] ## scale
    # b <- par[2] ## shape
    fgeneralnormal <- (b / (2 * pi * (a^2) * gamma(2 / b))) * exp(-(r^b / a^b))
    -sum(log(fgeneralnormal)) ##
  }
  n.se <- 30
  len <- 1000
  par.1.ini <- par.1.generalnormal - n.se * par.1.se.generalnormal
  if (par.1.ini <= 0) {
    par.1.ini <- 0.01
  }
  par.1.fin <- par.1.generalnormal + n.se * par.1.se.generalnormal
  par.1.est <- seq(par.1.ini, par.1.fin, length.out = len)

  par.1.prof <- numeric(len)
  for (i in 1:len) {
    possibleError <- tryCatch(
      par.1.prof[i] <- optim(log.dist.generalnormal.ci, par = par.2.generalnormal, a = par.1.est[i],
                             r = data,
                             method="Nelder-Mead")$value,
      error = function(e) e)
    if(!inherits(possibleError, "error")){
      par.1.prof[i] <- optim(log.dist.generalnormal.ci, par = par.2.generalnormal, a = par.1.est[i],
                             r = data,
                             method="Nelder-Mead")$value
    }
  }

  if (length(which(par.1.prof == 0) > 0)) {
    par.1.prof <- par.1.prof[-which(par.1.prof == 0)]
  }

  prof.lower <- par.1.prof[1:which.min(par.1.prof)]
  prof.par.1.lower <- par.1.est[1:which.min(par.1.prof)]

  prof.upper <- par.1.prof[which.min(par.1.prof):length(par.1.prof)]
  prof.par.1.upper <- par.1.est[which.min(par.1.prof):length(par.1.prof)]

  par.1.generalnormal.CIlow <- approx(prof.lower, prof.par.1.lower, xout = dist.generalnormal.opt$value + qchisq(confidence.level, 1)/2)$y
  par.1.generalnormal.CIupp <- approx(prof.upper, prof.par.1.upper, xout = dist.generalnormal.opt$value + qchisq(confidence.level, 1)/2)$y

  par.2.ini <- par.2.generalnormal - n.se * par.2.se.generalnormal
  if (par.2.ini <= 0) {
    par.2.ini <- 0.01
  }
  par.2.fin <- par.2.generalnormal + n.se * par.2.se.generalnormal
  par.2.est <- seq(par.2.ini , par.2.fin, length.out = len)

  par.2.prof <- numeric(len)
  for (i in 1:len) {
    possibleError <- tryCatch(
      par.2.prof[i] <- optim(log.dist.generalnormal.ci, par = par.1.generalnormal, b = par.2.est[i],
                             r = data,
                             method = "Nelder-Mead")$value,
      error = function(e) e)
    if(!inherits(possibleError, "error")){
      par.2.prof[i] <- optim(log.dist.generalnormal.ci, par = par.1.generalnormal, b = par.2.est[i],
                             r = data,
                             method = "Nelder-Mead")$value
    }
  }

  if (length(which(par.2.prof == 0) > 0)) {
    par.2.prof <- par.2.prof[-which(par.2.prof == 0)]
  }

  prof.lower <- par.2.prof[1:which.min(par.2.prof)]
  prof.par.2.lower <- par.2.est[1:which.min(par.2.prof)]

  prof.upper <- par.2.prof[which.min(par.2.prof):length(par.2.prof)]
  prof.par.2.upper <- par.2.est[which.min(par.2.prof):length(par.2.prof)]

  par.2.generalnormal.CIlow <- approx(prof.lower, prof.par.2.lower, xout = dist.generalnormal.opt$value + qchisq(confidence.level, 1)/2)$y
  par.2.generalnormal.CIupp <- approx(prof.upper, prof.par.2.upper, xout = dist.generalnormal.opt$value + qchisq(confidence.level, 1)/2)$y
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
  res <- data.frame(aic.generalnormal, aicc.generalnormal, bic.generalnormal,
                    chi.squared.statistic.generalnormal, chi.squared.pvalue.generalnormal,g.max.generalnormal, KS.generalnormal,
                    par.1.generalnormal, par.1.generalnormal.CIlow, par.1.generalnormal.CIupp,
                    par.2.generalnormal, par.2.generalnormal.CIlow, par.2.generalnormal.CIupp,
                    mean.generalnormal, mean.stderr.generalnormal, stdev.generalnormal, stdev.stderr.generalnormal,
                    skewness.generalnormal, skewness.stderr.generalnormal, kurtosis.generalnormal, kurtosis.stderr.generalnormal)
  generalnormal.values <- list("opt" = dist.generalnormal.opt, "res" = res)
}
