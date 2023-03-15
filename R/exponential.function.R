exponential.function <- function (data, chi.res.hist, ks.res.hist, confidence.level) {
  log.dist.exponential <- function (par, r) {
	a <- par[1]
	if(a < 0) return(Inf)

    fexponential <- (1 / (2 * pi * a ^ 2 )) * exp(-r/a) # corrected function, adapted from Nathan 2012
    -sum(log(fexponential))
  }
  dist.exponential <- function (r, a) {
    fexponential <-  2*pi*r*(1 / (2 * pi * a ^ 2 )) * exp(-r/a) # corrected function, adapted from Nathan 2012
  }
  # initial values estimation
  rate <- 1/mean(data)
  # optimization procedure
  dist.opt <- optim (par = rate,
                                 fn = log.dist.exponential,
                                 r = data,
                                 method = "Brent",
                                 lower = 1e-6,
                                 upper = 100000)
  # output values
  # AIC
  aic.exponential <- 2 + 2 * dist.opt$value
  # AICc
  aicc.exponential <- aic.exponential + (2 * length(dist.opt$par)^2 + 2 * length(dist.opt$par))/(length(data) - length(dist.opt$par) - 1 )
  # BIC
  bic.exponential <-  2 * dist.opt$value + length(dist.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.exponential <- dist.exponential(chi.res.hist$mids, dist.opt$par)*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.exponential <- sum((chi.res.hist$counts - chi.expected.values.exponential)^2 / chi.expected.values.exponential)
  chi.squared.pvalue.exponential <- 1-pchisq(chi.squared.statistic.exponential, length(chi.res.hist$counts)-2)
  # Kolmogorov-Smirnov
  ks.expected.values.exponential <- dist.exponential(ks.res.hist$mids, dist.opt$par)*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
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

  CI <- confint.dispfit(dist.opt, log.dist.exponential, data=data, lower=c(1e-6), upper=list(100000), confidence.level=confidence.level)

  # mean dispersal distance
  mean.exponential <- dist.exponential.opt$par*2
  mean.stderr.exponential <- sqrt(diag(solve(numDeriv::hessian(log.dist.exponential, x=dist.opt$par, r=data))))
  # variance
  variance.exponential <- 2*dist.exponential.opt$par^2
  variance.stderr.exponential <- msm::deltamethod(~ x1^2, mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.exponential, x=dist.opt$par, r=data)))
  # standard deviation
  stdev.exponential <- dist.opt$par
  stdev.stderr.exponential <- sqrt(diag(solve(numDeriv::hessian(log.dist.exponential, x=dist.opt$par, r=data))))
  # skewness
  skewness.exponential <- 6 * sqrt(2)
  skewness.stderr.exponential <- NA
  # kurtosis
  kurtosis.exponential <- 30
  kurtosis.stderr.exponential <- NA
  # output
  res <- data.frame(aic.exponential, aicc.exponential, bic.exponential,
                    chi.squared.statistic.exponential, chi.squared.pvalue.exponential, g.max.exponential, KS.exponential,
                    dist.opt$par[1], CI["par1.CIlow"], CI["par1.CIupp"],
                    dist.opt$par[2], CI["par2.CIlow"], CI["par2.CIupp"],
                    mean.exponential, mean.stderr.exponential, stdev.exponential, stdev.stderr.exponential,
                    skewness.exponential, skewness.stderr.exponential, kurtosis.exponential, kurtosis.stderr.exponential)
  exponential.values <- list("opt" = dist.opt, "res" = res)
}
