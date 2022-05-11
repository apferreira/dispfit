exponential.function <- function (data, chi.res.hist, ks.res.hist) {
  log.dist.exponential <- function (r, a) {
    fexponential <- (1 / (2 * pi * a ^ 2 )) * exp(-r/a) # corrected function, adapted from Nathan 2012
    -sum(log(fexponential))
  }
  dist.exponential <- function (r, a) {
    fexponential <-  2*pi*r*(1 / (2 * pi * a ^ 2 )) * exp(-r/a) # corrected function, adapted from Nathan 2012
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
  # output values
  # AIC
  aic.exponential <- 2 + 2 * dist.exponential.opt$value
  # AICc
  aicc.exponential <- aic.exponential + (2 * length(dist.exponential.opt$par)^2 + 2 * length(dist.exponential.opt$par))/(length(data) - length(dist.exponential.opt$par) - 1 )
  # BIC
  bic.exponential <-  2 * dist.exponential.opt$value + length(dist.exponential.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.exponential <- dist.exponential(chi.res.hist$mids, dist.exponential.opt$par)*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.exponential <- sum((chi.res.hist$counts - chi.expected.values.exponential)^2 / chi.expected.values.exponential)
  chi.squared.pvalue.exponential <- 1-pchisq(chi.squared.statistic.exponential, length(chi.res.hist$counts)-2)
  # Kolmogorov-Smirnov
  ks.expected.values.exponential <- dist.exponential(ks.res.hist$mids, dist.exponential.opt$par)*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
  simul.exponential <- c()
  for (i in seq_along(ks.res.hist$mids)) {
    simul.exponential <- c(simul.exponential, rep(ks.res.hist$mids[i], round(ks.expected.values.exponential[i], 0)))
  }
  ks.exponential <- ks.test(data, simul.exponential)
  g.max.exponential <- as.numeric(ks.exponential$statistic)
  KS.exponential <- as.numeric(ks.exponential$p.value)

  # cumulative.expected.values.exponential <- c(expected.values.exponential[1])
  # for (i in 1+seq_along(expected.values.exponential)) {
  #   cumulative.expected.values.exponential[i] <- cumulative.expected.values.exponential[i-1] + expected.values.exponential[i]
  # }
  # cumulative.expected.values.exponential <- cumulative.expected.values.exponential/sum(expected.values.exponential)
  # cumulative.expected.values.exponential <- cumulative.expected.values.exponential[!is.na(cumulative.expected.values.exponential)]
  # g.max.exponential <- max(abs(cumulative.data - cumulative.expected.values.exponential))
  # if (g.max.exponential < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
  #   KS.exponential <- "Accept"
  # } else {KS.exponential <- "Reject"}
  # parameter estimate
  par.1.exponential <- dist.exponential.opt$par
  par.2.exponential <- NA
  # parameter estimate standard error
  par.1.se.exponential <- sqrt(diag(solve(numDeriv::hessian(log.dist.exponential, x=dist.exponential.opt$par, r=data))))
  par.2.se.exponential <- NA
  # parameter estimate confidence intervals
  n.se <- 30
  len <- 1000
  par.1.ini <- par.1.exponential - n.se * par.1.se.exponential
  if (par.1.ini <= 0) {
    par.1.ini <- 0.01
  }
  par.1.fin <- par.1.exponential + n.se * par.1.se.exponential
  par.1.est <- seq(par.1.ini , par.1.fin, length = len)

  par.1.prof = numeric(len)
  for (i in 1:len) {
    par.1.prof[i] = log.dist.exponential(a = par.1.est[i],
                                         r = data)
  }

  prof.lower <- par.1.prof[1:which.min(par.1.prof)]
  prof.par.1.lower <- par.1.est[1:which.min(par.1.prof)]

  prof.upper <- par.1.prof[which.min(par.1.prof):length(par.1.prof)]
  prof.par.1.upper <- par.1.est[which.min(par.1.prof):length(par.1.prof)]

  par.1.exponential.CIlow <- approx(prof.lower, prof.par.1.lower, xout = dist.exponential.opt$value + qchisq(0.95, 1)/2)$y
  par.1.exponential.CIupp <- approx(prof.upper, prof.par.1.upper, xout = dist.exponential.opt$value + qchisq(0.95, 1)/2)$y

  par.2.exponential.CIlow <- NA
  par.2.exponential.CIupp <- NA
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
  res <- data.frame(aic.exponential, aicc.exponential, bic.exponential,
                    chi.squared.statistic.exponential, chi.squared.pvalue.exponential, g.max.exponential, KS.exponential,
                    par.1.exponential, par.1.exponential.CIlow, par.1.exponential.CIupp,
                    par.2.exponential, par.2.exponential.CIlow, par.2.exponential.CIupp,
                    mean.exponential, mean.stderr.exponential, stdev.exponential, stdev.stderr.exponential,
                    skewness.exponential, skewness.stderr.exponential, kurtosis.exponential, kurtosis.stderr.exponential)
  exponential.values <- list("opt" = dist.exponential.opt, "res" = res)
}
