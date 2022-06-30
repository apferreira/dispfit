confint.dispfit <- function(init.pars, logdistfun, data, lower.limits=c(0, 0), upper.limits=c(100000, 100000), confidence.level) {
	twopars <- length(init.pars$par) > 1
	
	par1 <- init.pars$par[1]
	# parameter estimate standard error
	par.1.se <- sqrt(diag(solve(numDeriv::hessian(logdistfun, x=init.pars$par, r=data))))[1]

	if(twopars) {
		par2 <- init.pars$par[2]
		par.2.se <- sqrt(diag(solve(numDeriv::hessian(logdistfun, x=init.pars$par, r=data))))[2]
	}
	
	# parameter estimate confidence intervals
	log.dist.ci <- function (a, b=NA, r) {
		return(logdistfun(par=c(a, b), r))
	}
	n.se <- 30
	len <- 1000

	# Par 1 CI
	par.1.ini <- par1 - n.se * par.1.se
	if (par.1.ini <= lower.limits[1]) {
		par.1.ini <- lower.limits[1] + 0.01
	}
	par.1.fin <- par1 + n.se * par.1.se
	par.1.est <- seq(par.1.ini, par.1.fin, length.out = len)

	par.1.prof <- numeric(len)
	if(twopars) {
		for (i in 1:len) {
			if(inherits(lower.limits, "list")) {
				if(inherits(lower.limits[[2]], "function"))
					par2.lower <- lower.limits[[2]](c(par.1.est[i], NA), data)
				else
					par2.lower <- lower.limits[[2]] + 0.00001
			} else if(inherits(lower.limits, "numeric")) {
				par2.lower <- lower.limits[2] + 0.00001
			}

			if(inherits(upper.limits, "list")) {
				if(inherits(upper.limits[[2]], "function"))
					par2.upper <- upper.limits[[2]](c(par.1.est[i], NA), data)
				else
					par2.upper <- upper.limits[[2]]
			} else if(inherits(upper.limits, "numeric")) {
				par2.upper <- upper.limits[2]
			}

			if(par2.upper > par2.lower && is.finite(log.dist.ci(par.1.est[i], par2, data))) {
				par.1.prof[i] <- optim(log.dist.ci, par = par2, a = par.1.est[i], r = data, lower=par2.lower, upper=par2.upper, method = "Brent")$value
			}
		}
	} else {
		for (i in 1:len) {
		    par.1.prof[i] = log.dist.ci(a = par.1.est[i], r = data)
	    }
	}

	if (length(which(par.1.prof == 0) > 0)) {
		par.1.prof <- par.1.prof[-which(par.1.prof == 0)]
	}

	prof.lower <- par.1.prof[1:which.min(par.1.prof)]
	prof.par.1.lower <- par.1.est[1:which.min(par.1.prof)]
	par.1.CIlow <- approx(prof.lower, prof.par.1.lower, xout = init.pars$value + qchisq(confidence.level, 1)/2)$y

	if(which.min(par.1.prof) == length(par.1.prof)) {
		par.1.CIupp <- Inf
	} else {
		prof.upper <- par.1.prof[which.min(par.1.prof):length(par.1.prof)]
		prof.par.1.upper <- par.1.est[which.min(par.1.prof):length(par.1.prof)]
		par.1.CIupp <- approx(prof.upper, prof.par.1.upper, xout = init.pars$value + qchisq(confidence.level, 1)/2)$y
	}
	

	if(twopars) {
		# Par 2 CI
		par.2.ini <- par2 - n.se * par.2.se
		if (par.2.ini <= lower.limits[2]) {
			par.2.ini <- lower.limits[2] + 0.01
		}
		par.2.fin <- par2 + n.se * par.2.se
		par.2.est <- seq(par.2.ini , par.2.fin, length.out = len)

		par.2.prof <- numeric(len)
		for (i in 1:len) {
			par1.lower <- lower.limits[1] + 0.00001
			
			if(inherits(lower.limits, "list")) {
				if(inherits(lower.limits[[1]], "function"))
					par1.lower <- lower.limits[[1]](c(NA, par.2.est[i]), data)
				else
					par1.lower <- lower.limits[[1]]
			} else if(inherits(lower.limits, "numeric")) {
				par1.lower <- lower.limits[1]
			}

			if(inherits(upper.limits, "list")) {
				if(inherits(upper.limits[[1]], "function"))
					par1.upper <- upper.limits[[1]](c(NA, par.2.est[i]), data)
				else
					par1.upper <- upper.limits[[1]]
			} else if(inherits(upper.limits, "numeric")) {
				par1.upper <- upper.limits[1]
			}

			if(par1.upper > par1.lower && is.finite(log.dist.ci(par1, par.2.est[i], data))) {
				par.2.prof[i] <- optim(log.dist.ci, par = par1, b = par.2.est[i], r = data, lower=par1.lower, upper=par1.upper, method = "Brent")$value
			}
		}

		if (length(which(par.2.prof == 0) > 0)) {
			par.2.prof <- par.2.prof[-which(par.2.prof == 0)]
		}

		prof.lower <- par.2.prof[1:which.min(par.2.prof)]
		prof.par.2.lower <- par.2.est[1:which.min(par.2.prof)]
		par.2.CIlow <- approx(prof.lower, prof.par.2.lower, xout = init.pars$value + qchisq(confidence.level, 1)/2)$y

		if(which.min(par.2.prof) == length(par.2.prof)) {
			par.2.CIupp <- Inf
		} else {
			prof.upper <- par.2.prof[which.min(par.2.prof):length(par.2.prof)]
			prof.par.2.upper <- par.2.est[which.min(par.2.prof):length(par.2.prof)]
			par.2.CIupp <- approx(prof.upper, prof.par.2.upper, xout = init.pars$value + qchisq(confidence.level, 1)/2)$y
		}
	} else {
	  par.2.CIlow <- NA
	  par.2.CIupp <- NA
	}

	return(c(
		par1.CIlow=par.1.CIlow,
		par1.CIupp=par.1.CIupp,
		par2.CIlow=par.2.CIlow,
		par2.CIupp=par.2.CIupp
	))
}
