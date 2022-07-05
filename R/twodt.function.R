twodt.function <- function (data, chi.res.hist, ks.res.hist, confidence.level) {
  log.dist.2dt <- function (r, par) {
    a <- par[1] ## scale parameter
    b <- par[2] ## shape parameter
    if(a < 0 || b < 1) return(Inf)
    
    f2dt <- ((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))

    if(any(is.nan(f2dt)) || any(f2dt < 0))
    	return(Inf)
    else
		return(-sum(log(f2dt)))
  }
  dist.2dt <- function (r, a, b) {
    f2dt <- 2*pi*r*((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))
  }
  # initial values estimation
 
    dist.opt <- optim (par = c(1, 1.000001), ## valor inicial para o "a"
                           fn = log.dist.2dt, ## função a minimizar
                           r = data,
                           method = "Nelder-Mead",
                           control = list(maxit = 10000))
  # dist.2dt.opt <- optim (par = c(0.001, 1.001), ## valor inicial para o "a"
  #                        fn = log.dist.2dt, ## função a minimizar
  #                        r = data, ## dados
  #                        control = list(maxit = 10000), ## limite superior do parametro
  #                        hessian = T)
  # output values
  # AIC
  aic.2dt <- 2*length(dist.opt$par) + 2 * dist.opt$value
  # AICc
  aicc.2dt <- aic.2dt + (2 * length(dist.opt$par)^2 + 2 * length(dist.opt$par))/(length(data) - length(dist.opt$par) - 1 )
  # BIC
  bic.2dt <-  2 * dist.opt$value + length(dist.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.2dt <- dist.2dt(chi.res.hist$mids, dist.opt$par[1], dist.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.2dt <- sum((chi.res.hist$counts - chi.expected.values.2dt)^2 / chi.expected.values.2dt)
  chi.squared.pvalue.2dt <- 1-pchisq(chi.squared.statistic.2dt, length(chi.res.hist$counts)-3)
  # Kolmogorov-Smirnov
  ks.expected.values.2dt <- dist.2dt(ks.res.hist$mids, dist.opt$par[1], dist.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
  simul.2dt <- c()
  for (i in seq_along(ks.res.hist$mids)) {
    simul.2dt <- c(simul.2dt, rep(ks.res.hist$mids[i], round(ks.expected.values.2dt[i], 0)))
  }
  ks.2dt <- ks.test(data, simul.2dt)
  g.max.2dt <- as.numeric(ks.2dt$statistic)
  KS.2dt <- as.numeric(ks.2dt$p.value)

  # cumulative.expected.values.2dt <- c(expected.values.2dt[1])
  # for (i in 1+seq_along(expected.values.2dt)) {
  #   cumulative.expected.values.2dt[i] <- cumulative.expected.values.2dt[i-1] + expected.values.2dt[i]
  # }
  # cumulative.expected.values.2dt <- cumulative.expected.values.2dt/sum(expected.values.2dt)
  # cumulative.expected.values.2dt <- cumulative.expected.values.2dt[!is.na(cumulative.expected.values.2dt)]
  # g.max.2dt <- max(abs(cumulative.data - cumulative.expected.values.2dt))
  # if (g.max.2dt < (sqrt(-log(0.01/2)/(2*length(cumulative.data))) * (1/(2*length(cumulative.data))))) {
  #   KS.2dt <- "Accept"
  # } else {KS.2dt <- "Reject"}

  # Confidence intervals
  # limits depend on data and pars, so define functions
  par2.upper.limit <- function(pars, data) {
	-log(.Machine$double.xmin) / log(1 + (max(data)^2)/(pars[1]^2))
  }
  
  CI <- confint.dispfit(dist.opt, log.dist.2dt, data=data, lower=c(1e-6, 1 + 1e-6), upper=list(100000, par2.upper.limit), confidence.level=confidence.level)

  # mean dispersal distance
  if (dist.opt$par[2] >= 3/2) {
    mean.2dt <- dist.opt$par[1] * (sqrt(pi)/2) * (exp(lgamma(dist.opt$par[2]-(3/2))-lgamma(dist.opt$par[2]-1)))
    mean.stderr.2dt <- msm::deltamethod(~ x1 * (sqrt(pi)/2)*(exp(lgamma(x2-(3/2))-lgamma(x2-1))), mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.opt$par, r=data)) )
  } else {
    mean.2dt <- Inf
    mean.stderr.2dt <- Inf
  }

  if (dist.opt$par[2] > 2) {
  # variance
    variance.2dt <- dist.opt$par[1] * sqrt(1/(dist.opt$par[2]-2))
    variance.stderr.2dt <- msm::deltamethod(~ x1 * sqrt(1/(x2-2)), mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.opt$par, r=data)) )
  # standard deviation
    stdev.2dt <- sqrt(dist.opt$par[1] * sqrt(1/(dist.opt$par[2]-2)))
    stdev.stderr.2dt <- msm::deltamethod(~ sqrt(x1 * sqrt(1/(x2-2))), mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.opt$par, r=data)) )
  } else {
    variance.2dt <- Inf
    variance.stderr.2dt <- Inf
    stdev.2dt <- Inf
    stdev.stderr.2dt <- Inf
  }

  # skewness
  if (2*(dist.opt$par[2]-1) > 2.5) {
    skewness.2dt <- (dist.opt$par[1]*(gamma(2.5) * exp(lgamma(dist.opt$par[2]-2.5) - lgamma(dist.opt$par[2]-1)))^(1/3))
    skewness.stderr.2dt <- msm::deltamethod(~ x1 * (gamma(2.5) * exp(lgamma(x2-2.5) - lgamma(x2-1)))^(1/3), mean = dist.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.opt$par, r=data)) )
  } else {
    skewness.2dt <- Inf
    skewness.stderr.2dt <- Inf
  }

  # kurtosis
  if (2*(dist.opt$par[2]-1) > 3) {
    kurtosis.2dt <- 2^(1/4) * dist.opt$par[1] * (1/((dist.opt$par[2] - 3) * (dist.opt$par[2] - 2)))^(1/4)
    kurtosis.stderr.2dt <- msm::deltamethod(~ 2^(1/4) * x1 * (1/((x2 - 3) * (x2 - 2)))^(1/4),
                                            mean = dist.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.opt$par, r=data)) )
  } else {
    kurtosis.2dt <- Inf
    kurtosis.stderr.2dt <- Inf
  }

  # output
  res <- data.frame(aic.2dt, aicc.2dt, bic.2dt,
                    chi.squared.statistic.2dt, chi.squared.pvalue.2dt, g.max.2dt, KS.2dt,
                    dist.opt$par[1], CI["par1.CIlow"], CI["par1.CIupp"],
                    dist.opt$par[2], CI["par2.CIlow"], CI["par2.CIupp"],
                    mean.2dt, mean.stderr.2dt, stdev.2dt, stdev.stderr.2dt,
                    skewness.2dt, skewness.stderr.2dt, kurtosis.2dt, kurtosis.stderr.2dt)
  twodt.values <- list("opt" = dist.opt, "res" = res)
}
