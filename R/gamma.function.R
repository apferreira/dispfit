gamma.function <- function (data, chi.res.hist, ks.res.hist, confidence.level) {
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
  # output values
  # AIC
  aic.gamma <- 2 * length(dist.gamma.opt$par) + 2 * dist.gamma.opt$value
  # AICc
  aicc.gamma <- aic.gamma + (2 * length(dist.gamma.opt$par)^2 + 2 * length(dist.gamma.opt$par))/(length(data) - length(dist.gamma.opt$par) - 1 )
  # BIC
  bic.gamma <-  2 * dist.gamma.opt$value + length(dist.gamma.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.gamma <- dist.gamma(chi.res.hist$mids, dist.gamma.opt$par[1],  dist.gamma.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.gamma <- sum((chi.res.hist$counts - chi.expected.values.gamma)^2 / chi.expected.values.gamma)
  chi.squared.pvalue.gamma <- 1-pchisq(chi.squared.statistic.gamma, length(chi.res.hist$counts)-3)
  # Kolmogorov-Smirnov
  ks.expected.values.gamma <- dist.gamma(ks.res.hist$mids, dist.gamma.opt$par[1],  dist.gamma.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
  simul.gamma <- c()
  for (i in seq_along(ks.res.hist$mids)) {
    simul.gamma <- c(simul.gamma, rep(ks.res.hist$mids[i], round(ks.expected.values.gamma[i], 0)))
  }
  ks.gamma <- ks.test(data, simul.gamma)
  g.max.gamma <- as.numeric(ks.gamma$statistic)
  KS.gamma <- as.numeric(ks.gamma$p.value)

  # cumulative.expected.values.gamma <- c(expected.values.gamma[1])
  # for (i in 1+seq_along(expected.values.gamma)) {
  #   cumulative.expected.values.gamma[i] <- cumulative.expected.values.gamma[i-1] + expected.values.gamma[i]
  # }
  # cumulative.expected.values.gamma <- cumulative.expected.values.gamma/sum(expected.values.gamma)
  # cumulative.expected.values.gamma <- cumulative.expected.values.gamma[!is.na(cumulative.expected.values.gamma)]
  # g.max.gamma <- max(abs(cumulative.data - cumulative.expected.values.gamma))
  # if (g.max.gamma < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
  #   KS.gamma <- "Accept"
  # } else {KS.gamma <- "Reject"}
  # parameter estimate
  par.1.gamma <- dist.gamma.opt$par[1]
  par.2.gamma <- dist.gamma.opt$par[2]
  # parameter estimate standard error
  par.1.se.gamma <- sqrt(diag(solve(numDeriv    ::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data))))[1]
  par.2.se.gamma <- sqrt(diag(solve(numDeriv::hessian(log.dist.gamma, x=dist.gamma.opt$par, r=data))))[2]
  # parameter estimate confidence intervals
  log.dist.gamma.ci <- function (a, b, r) {
    # a <- par[1] ## scale
    # b <- par[2] ## shape
    fgamma <- (1 / (2 * pi * (a^2) * gamma(b))) * ((r/a)^(b-2)) * exp(-r/a)
    -sum(log(fgamma)) ## Negative Log Likelihood
  }
  n.se <- 30
  len <- 1000
  par.1.ini <- par.1.gamma - n.se * par.1.se.gamma
  if (par.1.ini <= 0) {
    par.1.ini <- 0.01
  }
  par.1.fin <- par.1.gamma + n.se * par.1.se.gamma
  par.1.est <- seq(par.1.ini, par.1.fin, length.out = len)

  par.1.prof = numeric(len)
  for (i in 1:len) {
    possibleError <- tryCatch(
      par.1.prof[i] <- optim(log.dist.gamma.ci, par = par.2.gamma, a = par.1.est[i],
                             r = data,
                             method = "Nelder-Mead")$value,
      error = function(e) e)
    if(!inherits(possibleError, "error")){
      par.1.prof[i] <- optim(log.dist.gamma.ci, par = par.2.gamma, a = par.1.est[i],
                             r = data,
                             method = "Nelder-Mead")$value
    }
  }

  if (length(which(par.1.prof == 0) > 0)) {
    par.1.prof <- par.1.prof[-which(par.1.prof == 0)]
  }

  prof.lower <- par.1.prof[1:which.min(par.1.prof)]
  prof.par.1.lower <- par.1.est[1:which.min(par.1.prof)]

  prof.upper <- par.1.prof[which.min(par.1.prof):length(par.1.prof)]
  prof.par.1.upper <- par.1.est[which.min(par.1.prof):length(par.1.prof)]

  par.1.gamma.CIlow <- approx(prof.lower, prof.par.1.lower, xout = dist.gamma.opt$value + qchisq(confidence.level, 1)/2)$y
  par.1.gamma.CIupp <- approx(prof.upper, prof.par.1.upper, xout = dist.gamma.opt$value + qchisq(confidence.level, 1)/2)$y

  par.2.ini <- par.2.gamma - n.se * par.2.se.gamma
  if (par.2.ini <= 0) {
    par.2.ini <- 0.01
  }
  par.2.fin <- par.2.gamma + n.se * par.2.se.gamma
  par.2.est <- seq(par.2.ini , par.2.fin, length.out = len)

  par.2.prof = numeric(len)
  for (i in 1:len) {
    possibleError <- tryCatch(
      par.2.prof[i] <- optim(log.dist.gamma.ci, par = par.1.gamma, b = par.2.est[i],
                             r = data,
                             method = "Nelder-Mead")$value,
      error = function(e) e)
    if(!inherits(possibleError, "error")){
      par.2.prof[i] <- optim(log.dist.gamma.ci, par = par.1.gamma, b = par.2.est[i],
                             r = data,
                             method = "Nelder-Mead")$value
    }
  }

  if (length(which(par.2.prof == 0) > 0)) {
    par.2.prof <- par.2.prof[-which(par.2.prof == 0)]
  }

  prof.lower = par.2.prof[1:which.min(par.2.prof)]
  prof.par.2.lower = par.2.est[1:which.min(par.2.prof)]

  prof.upper <- par.2.prof[which.min(par.2.prof):length(par.2.prof)]
  prof.par.2.upper <- par.2.est[which.min(par.2.prof):length(par.2.prof)]

  par.2.gamma.CIlow <- approx(prof.lower, prof.par.2.lower, xout = dist.gamma.opt$value + qchisq(confidence.level, 1)/2)$y
  par.2.gamma.CIupp <- approx(prof.upper, prof.par.2.upper, xout = dist.gamma.opt$value + qchisq(confidence.level, 1)/2)$y
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
  # standard deviation
  stdev.gamma <- sqrt(dist.gamma.opt$par[1]^2 * dist.gamma.opt$par[2])
  stdev.stderr.gamma <- msm::deltamethod(~ sqrt(x2*x1^2),
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
  res <- data.frame(aic.gamma, aicc.gamma, bic.gamma,
                    chi.squared.statistic.gamma, chi.squared.pvalue.gamma,g.max.gamma, KS.gamma,
                    par.1.gamma, par.1.gamma.CIlow, par.1.gamma.CIupp, par.2.gamma, par.2.gamma.CIlow, par.2.gamma.CIupp,
                    mean.gamma, mean.stderr.gamma, stdev.gamma, stdev.stderr.gamma,
                    skewness.gamma, skewness.stderr.gamma, kurtosis.gamma, kurtosis.stderr.gamma)
  gamma.values <- list("opt" = dist.gamma.opt, "res" = res)
}
