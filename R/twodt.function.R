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
 
    dist.2dt.opt <- optim (par = c(1, 1.000001), ## valor inicial para o "a"
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
  aic.2dt <- 2*length(dist.2dt.opt$par) + 2 * dist.2dt.opt$value
  # AICc
  aicc.2dt <- aic.2dt + (2 * length(dist.2dt.opt$par)^2 + 2 * length(dist.2dt.opt$par))/(length(data) - length(dist.2dt.opt$par) - 1 )
  # BIC
  bic.2dt <-  2 * dist.2dt.opt$value + length(dist.2dt.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.2dt <- dist.2dt(chi.res.hist$mids, dist.2dt.opt$par[1], dist.2dt.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.2dt <- sum((chi.res.hist$counts - chi.expected.values.2dt)^2 / chi.expected.values.2dt)
  chi.squared.pvalue.2dt <- 1-pchisq(chi.squared.statistic.2dt, length(chi.res.hist$counts)-3)
  # Kolmogorov-Smirnov
  ks.expected.values.2dt <- dist.2dt(ks.res.hist$mids, dist.2dt.opt$par[1], dist.2dt.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
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
  # parameter estimate
  par.1.2dt <- dist.2dt.opt$par[1]
  par.2.2dt <- dist.2dt.opt$par[2]
  # parameter estimate standard error
  par.1.se.2dt <- sqrt(diag(solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data))))[1]
  par.2.se.2dt <- sqrt(diag(solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data))))[2]
  # parameter estimate confidence intervals
  log.dist.2dt.ci <- function (r, a, b) {
    # a <- par[1] ## scale parameter
    # b <- par[2] ## shape parameter
    if(a < 0 || b < 1) return(Inf)
    
    f2dt <- ((b-1) / (pi*(a^2))) * ((1 + (r^2)/(a^2))^(-b))

    if(any(is.nan(f2dt)) || any(f2dt < 0))
    	return(Inf)
    else
		return(-sum(log(f2dt)))
  }
  n.se <- 30
  len <- 1000
  par.1.ini <- par.1.2dt - n.se * par.1.se.2dt
  if (par.1.ini <= 0) {
    par.1.ini <- 0.01
  }
  par.1.fin <- par.1.2dt + n.se * par.1.se.2dt
  par.1.est <- seq(par.1.ini, par.1.fin, length.out = len)

  par.1.prof = numeric(len)
  
   
  for (i in 1:len) {
  	# upper limit:
	# (1 + (r^2)/(a^2))^(-b) = M
	# b = -log(M) / log(1 + (r^2)/(a^2))
	  par.1.prof[i] <- optim(log.dist.2dt.ci, par = par.2.2dt, a = par.1.est[i],
                         r = data, lower=1.00001, upper=-log(.Machine$double.xmin) / log(1 + (max(data)^2)/(par.1.est[i]^2)), method = "Brent")$value
  }

  if (length(which(par.1.prof == 0) > 0)) {
    par.1.prof <- par.1.prof[-which(par.1.prof == 0)]
  }

  prof.lower <- par.1.prof[1:which.min(par.1.prof)]
  prof.par.1.lower <- par.1.est[1:which.min(par.1.prof)]
  par.1.2dt.CIlow <- approx(prof.lower, prof.par.1.lower, xout = dist.2dt.opt$value + qchisq(confidence.level, 1)/2)$y

  if(which.min(par.1.prof) == length(par.1.prof)) {
    par.1.2dt.CIupp <- Inf
  } else {
    prof.upper <- par.1.prof[which.min(par.1.prof):length(par.1.prof)]
    prof.par.1.upper <- par.1.est[which.min(par.1.prof):length(par.1.prof)]
    par.1.2dt.CIupp <- approx(prof.upper, prof.par.1.upper, xout = dist.2dt.opt$value + qchisq(confidence.level, 1)/2)$y
  }

  par.2.ini <- par.2.2dt - n.se * par.2.se.2dt
  if (par.2.ini <= 1) {
    par.2.ini <- 1.01
  }
  par.2.fin <- par.2.2dt + n.se * par.2.se.2dt
  par.2.est <- seq(par.2.ini , par.2.fin, length.out = len)

  par.2.prof = numeric(len)
  for (i in 1:len) {
      par.2.prof[i] <- optim(log.dist.2dt.ci, par = par.1.2dt, b = par.2.est[i],
                             r = data, lower=0.00001, upper=10000,
                             method = "Brent")$value
  }

  if (length(which(par.2.prof == 0) > 0)) {
    par.2.prof <- par.2.prof[-which(par.2.prof == 0)]
  }

  prof.lower = par.2.prof[1:which.min(par.2.prof)]
  prof.par.2.lower = par.2.est[1:which.min(par.2.prof)]
  par.2.2dt.CIlow <- approx(prof.lower, prof.par.2.lower, xout = dist.2dt.opt$value + qchisq(confidence.level, 1)/2)$y

  if(which.min(par.2.prof) == length(par.2.prof)) {
    par.2.2dt.CIupp <- Inf
  } else {
    prof.upper <- par.2.prof[which.min(par.2.prof):length(par.2.prof)]
    prof.par.2.upper <- par.2.est[which.min(par.2.prof):length(par.2.prof)]
    par.2.2dt.CIupp <- approx(prof.upper, prof.par.2.upper, xout = dist.2dt.opt$value + qchisq(confidence.level, 1)/2)$y
  }
  # mean dispersal distance
  if (dist.2dt.opt$par[2] >= 3/2) {
    mean.2dt <- dist.2dt.opt$par[1] * (sqrt(pi)/2) * (exp(lgamma(dist.2dt.opt$par[2]-(3/2))-lgamma(dist.2dt.opt$par[2]-1)))
  } else {
    mean.2dt <- "Infinite Value"
  }
  if (dist.2dt.opt$par[2] >= 1.5) {
    mean.stderr.2dt <- msm::deltamethod(~ x1 * (sqrt(pi)/2)*(exp(lgamma(x2-(3/2))-lgamma(x2-1))), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    mean.stderr.2dt <- "Infinite Value"
  }
  # variance
  if (dist.2dt.opt$par[2] > 2) {
    variance.2dt <- dist.2dt.opt$par[1] * sqrt(1/(dist.2dt.opt$par[2]-2))
  } else {
    variance.2dt <- "Infinite Value"
  }
  if (dist.2dt.opt$par[2] > 2) {
    variance.stderr.2dt <- msm::deltamethod(~ x1 * sqrt(1/(x2-2)), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    variance.stderr.2dt <- "Infinite Value"
  }
  # standard deviation
  if (dist.2dt.opt$par[2] > 2) {
    stdev.2dt <- sqrt(dist.2dt.opt$par[1] * sqrt(1/(dist.2dt.opt$par[2]-2)))
  } else {
    stdev.2dt <- "Infinite Value"
  }
  if (dist.2dt.opt$par[2] > 2) {
    stdev.stderr.2dt <- msm::deltamethod(~ sqrt(x1 * sqrt(1/(x2-2))), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    stdev.stderr.2dt <- "Infinite Value"
  }
  # skewness
  if (2*(dist.2dt.opt$par[2]-1) > 2.5) {
    skewness.2dt <- (dist.2dt.opt$par[1]*(gamma(2.5) * exp(lgamma(dist.2dt.opt$par[2]-2.5) - lgamma(dist.2dt.opt$par[2]-1)))^(1/3))
  } else {
    skewness.2dt <- "Infinite Value"
  }
  if (2*(dist.2dt.opt$par[2]-1) > 2.5) {
    skewness.stderr.2dt <- msm::deltamethod(~ x1 * (gamma(2.5) * exp(lgamma(x2-2.5) - lgamma(x2-1)))^(1/3), mean = dist.2dt.opt$par, cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    skewness.stderr.2dt <- "Infinite Value"
  }
  # kurtosis
  if (2*(dist.2dt.opt$par[2]-1) > 3) {
    kurtosis.2dt <- 2^(1/4) * dist.2dt.opt$par[1] * (1/((dist.2dt.opt$par[2] - 3) * (dist.2dt.opt$par[2] - 2)))^(1/4)
  } else {
    kurtosis.2dt <- "Infinite Value"
  }
  if (2*(dist.2dt.opt$par[2]-1) > 3) {
    kurtosis.stderr.2dt <- msm::deltamethod(~ 2^(1/4) * x1 * (1/((x2 - 3) * (x2 - 2)))^(1/4),
                                            mean = dist.2dt.opt$par,
                                            cov = solve(numDeriv::hessian(log.dist.2dt, x=dist.2dt.opt$par, r=data)) )
  } else {
    kurtosis.stderr.2dt <- "Infinite Value"
  }
  # output
  res <- data.frame(aic.2dt, aicc.2dt, bic.2dt,
                    chi.squared.statistic.2dt, chi.squared.pvalue.2dt, g.max.2dt, KS.2dt,
                    par.1.2dt, par.1.2dt.CIlow, par.1.2dt.CIupp,
                    par.2.2dt, par.2.2dt.CIlow, par.2.2dt.CIupp,
                    mean.2dt, mean.stderr.2dt, stdev.2dt, stdev.stderr.2dt,
                    skewness.2dt, skewness.stderr.2dt, kurtosis.2dt, kurtosis.stderr.2dt)
  twodt.values <- list("opt" = dist.2dt.opt, "res" = res)
}
