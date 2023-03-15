generalnormal.function <- function (data, chi.res.hist, ks.res.hist, confidence.level) {
  log.dist.generalnormal <- function (r, par) {
    a <- par[1] ## scale
    b <- par[2] ## shape
    if(a < 0 || b < 0) return(Inf)

    fgeneralnormal <- (b / (2 * pi * (a^2) * gamma(2 / b))) * exp(-(r^b / a^b))

    if(any(is.nan(fgeneralnormal)) || any(fgeneralnormal < 0))
    	return(NaN)
    else
		return(-sum(log(fgeneralnormal)))
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
  dist.opt <- optim (par = c(scale, shape), ## initial values
                                   fn = log.dist.generalnormal, ## minimizing function
                                   r = data,
                                   method = "Nelder-Mead",
                                   control = list(maxit = 10000))

  if (dist.opt$par[1] <= 0 || dist.opt$par[2] <= 0) stop("Unknown error condition, negative parameters")

  #   # output values
  # AIC
  aic.generalnormal <- 2*length(dist.opt$par) + 2 * dist.opt$value ## AIC
  # AICc
  aicc.generalnormal <- aic.generalnormal + (2 * length(dist.opt$par)^2 + 2 * length(dist.opt$par))/(length(data) - length(dist.opt$par) - 1 )
  # BIC
  bic.generalnormal <-  2 * dist.opt$value + length(dist.opt$par)*log(length(data))
  # Chi-squared
  chi.expected.values.generalnormal <- dist.generalnormal(chi.res.hist$mids, dist.opt$par[1], dist.opt$par[2])*length(data)*(chi.res.hist$breaks[2] - chi.res.hist$breaks[1])
  chi.squared.statistic.generalnormal <- sum((chi.res.hist$counts - chi.expected.values.generalnormal)^2 / chi.expected.values.generalnormal)
  chi.squared.pvalue.generalnormal <- 1-pchisq(chi.squared.statistic.generalnormal, length(chi.res.hist$counts)-3)
  # Kolmogorov-Smirnov
  ks.expected.values.generalnormal <- dist.generalnormal(ks.res.hist$mids, dist.opt$par[1], dist.opt$par[2])*length(data)*(ks.res.hist$breaks[2] - ks.res.hist$breaks[1])
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


  # Confidence intervals
  # limits depend on data and pars, so define functions
  par1.upper.limit <- function(pars, data) {
  	sqrt(.Machine$double.xmax / (gamma(2 / pars[2]) * 2 * pi))
  }

  par2.upper.limit <- function(pars, data) {
	# b=par2; a=par.1.est[i]; r=sort(data); (b / (2 * pi * (a^2) * gamma(2 / b))) * exp(-(r^b / a^b))

	# Upper limits for b (minimum of)
	# * exp(-(r^b / a^b)) = m
	#   -(r/a)^b = log(m)
	#   b = log(-log(m))/log(r/a)

	# * r^b = M
	#   b = log(M) / log(r)

	min(
		log(.Machine$double.xmax) / log(max(data)),
		log(-log(.Machine$double.xmin)) / log(max(data) / pars[1])
	)
  }


  CI <- confint.dispfit(dist.opt, log.dist.generalnormal, data=data, lower=c(1e-6, 1e-6), upper=list(par1.upper.limit, par2.upper.limit), confidence.level=confidence.level)

  # mean dispersal distance
  mean.generalnormal <- dist.opt$par[1] * (gamma(3/dist.opt$par[2])/gamma(2/dist.opt$par[2]))
  # mean.stderr.generalnormal <- msm::deltamethod(~ x1 * (gamma(3/x2)/gamma(2/x2)),
  #                                               mean = dist.opt$par,
  #                                               cov = solve(numDeriv::hessian(log.dist.generalnormal,
  #                                                                             x=dist.opt$par,
  #                                                                             r=data)) )
  # variance
  variance.generalnormal <- (dist.opt$par[1]^2 *gamma(4/dist.opt$par[2])) / gamma(2/dist.opt$par[2]) - mean.generalnormal^2
  # variance.stderr.generalnormal <- msm::deltamethod(~ (x1^2 *gamma(4/x2)) / gamma(2/x2) - (x1 * (gamma(3/x2)/gamma(2/x2)))^2,
  #                                                   mean = dist.opt$par,
  #                                                   cov = solve(numDeriv::hessian(log.dist.generalnormal,
  #                                                                                 x=dist.opt$par,
  #                                                                                 r=data)) )
  # standard deviation
  stdev.generalnormal <- sqrt((dist.opt$par[1]^2 *gamma(4/dist.opt$par[2])) / gamma(2/dist.opt$par[2]) - mean.generalnormal^2)
  # stdev.stderr.generalnormal <- msm::deltamethod(~ sqrt((x1^2 *gamma(4/x2)) / gamma(2/x2) - (x1 * (gamma(3/x2)/gamma(2/x2)))^2),
  #                                                mean = dist.opt$par,
  #                                                cov = solve(numDeriv::hessian(log.dist.generalnormal,
  #                                                                              x=dist.opt$par,
  #                                                                              r=data)) )
  # skewness
  skewness.generalnormal <- ((dist.opt$par[1]^3 * gamma(5/dist.opt$par[2])) / gamma(2/dist.opt$par[2])) / (variance.generalnormal^(3/2))
  # skewness.stderr.generalnormal <- msm::deltamethod(~ ((x1^3 * gamma(5/x2)) / gamma(2/x2)) / ((x1^2 *gamma(4/x2)) / gamma(2/x2) - (x1 * (gamma(3/x2)/gamma(2/x2)))^2)^(3/2),
  #                                                   mean = dist.opt$par,
  #                                                   cov = solve(numDeriv::hessian(log.dist.generalnormal,
  #                                                                                 x=dist.opt$par,
  #                                                                                 r=data)) )
  # kurtosis
  kurtosis.generalnormal <- ((dist.opt$par[1]^4 * gamma(6/dist.opt$par[2])) / gamma(2/dist.opt$par[2])) / (variance.generalnormal^2)
  # (gamma(6/dist.opt$par[2]) * gamma(2/dist.opt$par[2])) / (gamma(4/dist.opt$par[2])^2)
  # ((dist.opt$par[1]^4 * gamma(6/dist.opt$par[2])) / gamma(2/dist.opt$par[2])) / ((dist.opt$par[1]^2 *gamma(4/dist.opt$par[2])) / gamma(2/dist.opt$par[2]))^2
  # kurtosis.stderr.generalnormal <- msm::deltamethod(~ ((x1^4 * gamma(6/x2)) / gamma(2/x2)) / ((x1^2 *gamma(4/x2)) / gamma(2/x2) - (x1 * (gamma(3/x2)/gamma(2/x2)))^2)^2,
  #                                                   mean = dist.opt$par,
  #                                                   cov = solve(numDeriv::hessian(log.dist.generalnormal,
  #                                                                                 x=dist.opt$par,
  #                                                                                 r=data)) )
  # output
  res <- data.frame(aic.generalnormal, aicc.generalnormal, bic.generalnormal,
                    chi.squared.statistic.generalnormal, chi.squared.pvalue.generalnormal,g.max.generalnormal, KS.generalnormal,
                    dist.opt$par[1], CI["par1.CIlow"], CI["par1.CIupp"],
                    dist.opt$par[2], CI["par2.CIlow"], CI["par2.CIupp"],
                    mean.generalnormal, mean.stderr.generalnormal, stdev.generalnormal, stdev.stderr.generalnormal,
                    skewness.generalnormal, skewness.stderr.generalnormal, kurtosis.generalnormal, kurtosis.stderr.generalnormal)
  generalnormal.values <- list("opt" = dist.opt, "res" = res)
}
