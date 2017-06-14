
############################################################################
#         R FUNCTIONS FOR THE INVERSE GAMMA DISTRIBUTION
#
#         (by Jeffrey S. Rosenthal, probability.ca, 2007)
# 
# The base package of the statistical software R does not seem to include
# built-in functions for the inverse gamma distribution.  So, I provide
# them here, using trivial modifications of the corresponding
# built-in functions for the gamma distribution.  (Released under GPL.)
############################################################################

# dinvgamma: density function of the inverse-gamma distribution
# ref wiki alpha = shape, beta = scale 
dinvgamma = function(x, shape = 1, rate = 1, scale = 1/rate, log = FALSE) {
    # return( scale^shape / gamma(shape) * exp(-scale/x) * x^(-shape-1) )
    logval = shape*log(scale) - lgamma(shape) - scale/x - (shape+1)*log(x)
    if (log)
	return(logval)
    else
	return(exp(logval))
}

# pinvgamma: cumulative distribution function of the inverse-gamma distribution. 
# There is error in the version downloaded
pinvgamma = function(q, shape = 1, rate = 1, scale = 1/rate,
			lower.tail = TRUE, log.p = FALSE) {
    return( pgamma(q=1/q,shape=shape,scale=1/scale,lower.tail=!lower.tail, log.p=log.p) )
}

# qinvgamma: quantile function of the inverse-gamma distribution
qinvgamma = function(p, shape = 1, rate = 1, scale = 1/rate,
			lower.tail = TRUE, log.p = FALSE) {
  stopifnot(FALSE)
    return( 1 / qgamma(p=p,shape=shape,scale=1/scale,lower.tail=!lower.tail, log.p=log.p) )
}

# rinvgamma: sample a vector of inverse-gamma random variables
rinvgamma = function(n, shape = 1, rate = 1, scale = 1/rate) {
    return( 1 / rgamma(n,shape=shape,scale=1/scale) )
}

