rayleigh.function <- function (data, chi.res.hist, ks.res.hist, confidence.level) {
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
  # output values
  # AIC
  aic.rayleigh <- 2 * length(dist.rayleigh.opt$par) + 2 * dist.rayleigh.opt$value
  # AICc
  aicc.rayleigh <- aic.rayleigh + (2 * length(dist.rayleigh.opt$par)^2 + 2 * length(dist.rayleigh.opt$par))/(length(data) - length(dist.rayleigh.opt$par) - 1 )
  # BIC
  bic.rayleigh <-  2 * dist.rayleigh.opt$value + length(dist.rayleigh.opt$par)*log(length(data))

  # Chi-squared - from Press 1992, pp. 621-622
  chi.expected.values.rayleigh <- dist.rayleigh(chi.res.hist$mids,dist.rayleigh.opt$par)*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.rayleigh <- sum((chi.res.hist$counts - chi.expected.values.rayleigh)^2 / chi.expected.values.rayleigh)
  chi.squared.pvalue.rayleigh <- 1 - pchisq(chi.squared.statistic.rayleigh, length(chi.res.hist$counts)-2)
  # Kolmogorov-Smirnov - from Sokal 1995, pp. 223-224
  ks.expected.values.rayleigh <- dist.rayleigh(ks.res.hist$mids,dist.rayleigh.opt$par)*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
  simul.rayleigh <- c()
  for (i in seq_along(ks.res.hist$mids)) {
    simul.rayleigh <- c(simul.rayleigh, rep(ks.res.hist$mids[i], round(ks.expected.values.rayleigh[i], 0)))
  }
  ks.rayleigh <- ks.test(data, simul.rayleigh)
  ks.d.rayleigh <- as.numeric(ks.rayleigh$statistic)
  ks.p.rayleigh <- as.numeric(ks.rayleigh$p.value)

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
  # parameter estimate confidence intervals
  n.se <- 30
  len <- 1000
  par.1.ini <- par.1.rayleigh - n.se * par.1.se.rayleigh
  if (par.1.ini <= 0) {
    par.1.ini <- 0.01
  }
  par.1.fin <- par.1.rayleigh + n.se * par.1.se.rayleigh
  par.1.est <- seq(par.1.ini , par.1.fin, length.out = len)

  par.1.prof = numeric(1000)
  for (i in 1:1000) {
    par.1.prof[i] = log.dist.rayleigh(a = par.1.est[i],
                                      r = data)
  }

  prof.lower <- par.1.prof[1:which.min(par.1.prof)]
  prof.par.1.lower <- par.1.est[1:which.min(par.1.prof)]

  prof.upper <- par.1.prof[which.min(par.1.prof):length(par.1.prof)]
  prof.par.1.upper <- par.1.est[which.min(par.1.prof):length(par.1.prof)]

  par.1.rayleigh.CIlow <- approx(prof.lower, prof.par.1.lower, xout = dist.rayleigh.opt$value + qchisq(confidence.level, 1)/2)$y
  par.1.rayleigh.CIupp <- approx(prof.upper, prof.par.1.upper, xout = dist.rayleigh.opt$value + qchisq(confidence.level, 1)/2)$y

  par.2.rayleigh.CIlow <- NA
  par.2.rayleigh.CIupp <- NA
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
  res <- data.frame(aic.rayleigh, aicc.rayleigh, bic.rayleigh,
                    chi.squared.statistic.rayleigh, chi.squared.pvalue.rayleigh, ks.d.rayleigh, ks.p.rayleigh,
                    par.1.rayleigh, par.1.rayleigh.CIlow, par.1.rayleigh.CIupp, par.2.rayleigh, par.2.rayleigh.CIlow, par.2.rayleigh.CIupp,
                    mean.rayleigh, mean.stderr.rayleigh, stdev.rayleigh, stdev.stderr.rayleigh,
                    skewness.rayleigh, skewness.stderr.rayleigh, kurtosis.rayleigh, kurtosis.stderr.rayleigh)
  rayleigh.values <- list("opt" = dist.rayleigh.opt, "res" = res)
}
