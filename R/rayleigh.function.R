rayleigh.function <- function (data, chi.res.hist, ks.res.hist, confidence.level) {
  log.dist.rayleigh <- function (par, r) {
	a <- par[1]
	if(a < 0) return(Inf)

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
  dist.opt <- optim (par = sd0, ## initial value
                              fn = log.dist.rayleigh, ## function to minimize
                              r = data, ## dados
                              method = "Brent", ## one-parameter method
                              lower = 1e-6, ## lower-bound for parameter
                              upper = 100000)
  # output values
  # AIC
  aic.rayleigh <- 2 * length(dist.opt$par) + 2 * dist.opt$value
  # AICc
  aicc.rayleigh <- aic.rayleigh + (2 * length(dist.opt$par)^2 + 2 * length(dist.opt$par))/(length(data) - length(dist.opt$par) - 1 )
  # BIC
  bic.rayleigh <-  2 * dist.opt$value + length(dist.opt$par)*log(length(data))

  # Chi-squared - from Press 1992, pp. 621-622
  chi.expected.values.rayleigh <- dist.rayleigh(chi.res.hist$mids,dist.opt$par)*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.rayleigh <- sum((chi.res.hist$counts - chi.expected.values.rayleigh)^2 / chi.expected.values.rayleigh)
  chi.squared.pvalue.rayleigh <- 1 - pchisq(chi.squared.statistic.rayleigh, length(chi.res.hist$counts)-2)
  # Kolmogorov-Smirnov - from Sokal 1995, pp. 223-224
  ks.expected.values.rayleigh <- dist.rayleigh(ks.res.hist$mids,dist.opt$par)*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
  simul.rayleigh <- c()
  for (i in seq_along(ks.res.hist$mids)) {
    simul.rayleigh <- c(simul.rayleigh, rep(ks.res.hist$mids[i], round(ks.expected.values.rayleigh[i], 0)))
  }
  ks.rayleigh <- ks.test(data, simul.rayleigh)
  ks.d.rayleigh <- as.numeric(ks.rayleigh$statistic)
  ks.p.rayleigh <- as.numeric(ks.rayleigh$p.value)

  CI <- confint.dispfit(dist.opt, log.dist.rayleigh, data=data, lower=c(1e-6), upper=list(100000), confidence.level=confidence.level)

  # mean
  mean.rayleigh <- dist.opt$par*sqrt(pi)/2
  mean.stderr.rayleigh <- msm::deltamethod(~ x1 * sqrt(pi) / 2, mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.opt$par, r=data)))
  # variance
  variance.rayleigh <- ((4-pi)*(dist.opt$par^2))/4
  variance.stderr.rayleigh <- msm::deltamethod(~ ((4-pi)*(x1^2))/4, mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.opt$par, r=data)))
  # standard deviation
  stdev.rayleigh <- sqrt(((4-pi)*(dist.opt$par^2))/4)
  stdev.stderr.rayleigh <- msm::deltamethod(~ sqrt(((4-pi)*(x1^2))/4), mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.opt$par, r=data)))
  # skewness
  skewness.rayleigh <- (2*sqrt(pi) * (pi-3)) / ((4 - pi)^(3/2))
  skewness.stderr.rayleigh <- msm::deltamethod(~ (( (3*sqrt(pi)*x1^3) /4)) / (( (4-pi) * (x1^2) )/4)^(3/2), mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.opt$par, r=data)))
  # kurtosis
  kurtosis.rayleigh <- -(6*pi^2-24*pi+16)/(4-pi)^(3/2)
  kurtosis.stderr.rayleigh <- msm::deltamethod(~ (2*x1^4) / (( (4-pi) * (x1^2) )/4)^(2), mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.rayleigh, x=dist.opt$par, r=data)))
  # output
  res <- data.frame(aic.rayleigh, aicc.rayleigh, bic.rayleigh,
                    chi.squared.statistic.rayleigh, chi.squared.pvalue.rayleigh, ks.d.rayleigh, ks.p.rayleigh,
                    dist.opt$par[1], CI["par1.CIlow"], CI["par1.CIupp"],
                    dist.opt$par[2], CI["par2.CIlow"], CI["par2.CIupp"],
                    mean.rayleigh, mean.stderr.rayleigh, stdev.rayleigh, stdev.stderr.rayleigh,
                    skewness.rayleigh, skewness.stderr.rayleigh, kurtosis.rayleigh, kurtosis.stderr.rayleigh)
  rayleigh.values <- list("opt" = dist.opt, "res" = res)
}
