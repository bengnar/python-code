'''SPARSE STRF CODE'''
import numpy as np
import scipy.stats

def model_pvalue(wts, stim, resp, nboot=1e4, randinds=None):
    """Computes a bootstrap p-value by resampling the [wts] of the model, which
    is [wts] * [stim] ~ [resp].
    """
    origcorr = np.corrcoef(resp, np.dot(stim, wts))[0,1]
    if randinds is None:
        #randinds = np.random.randint(0, len(wts), (len(wts), nboot))
        randinds = make_randinds(len(wts), nboot)
    pwts = wts[randinds]
    pred = np.dot(stim, pwts)
    
    ## Compute correlations using vectorized method and bootstrap p-value
    zpred = (pred-pred.mean(0))/pred.std(0)
    zresp = (resp-resp.mean())/resp.std()
    bootcorrs = np.dot(zpred.T, zresp).ravel()/resp.shape[0]
    #bootcorrs = np.array([np.corrcoef(resp, p.T)[0,1] for p in pred.T])
    bspval = np.mean(bootcorrs>origcorr)
    
    ## Compute parametric p-value based on transformed distribution
    zccs = ztransformccs(bootcorrs)
    zorig = ztransformccs(origcorr)
    ppval = 1-scipy.stats.norm.cdf(zorig, loc=zccs.mean(), scale=zccs.std())
    
    print "Boostrap p-value: %0.3f, parametric p-value: %0.03f"%(bspval, ppval)
    return bspval, ppval

def make_randinds(nwts, nboot, algo="randint"):
    if algo=="randint":
        return np.random.randint(0, nwts, (nwts, nboot))
    
    elif algo=="bytes":
        N = nwts*nboot*2
        return np.mod(np.frombuffer(np.random.bytes(N), dtype=np.uint16), nwts).reshape((nwts, nboot))

def ztransformccs(ccs):
    """Transforms the given correlation coefficients to be vaguely Gaussian.
    """
    return ccs/np.sqrt((1-ccs**2))

def exact_correlation_pvalue(corr, N, alt="greater"):
    """Returns the exact p-value for the correlation between [a] and [b].
    The null hypothesis is that the correlation is zero. The distribution of
    correlation coefficients given that the true correlation is zero and both
    [a] and [b] are gaussian is given at 
    http://en.wikipedia.org/wiki/Pearson_correlation#Exact_distribution_for_Gaussian_data

    Parameters
    ----------
    a : array_like, shape (N,)
    b : array_like, shape (N,)
    alt : string
        The alternative hypothesis, is the correlation 'greater' than zero,
        'less' than zero, or just 'nonzero'.
        
    Returns
    -------
    pval : float
        Probability of sample correlation between [a] and [b] if actual correlation
        is zero.
    """
    #assert a.size == b.size
    #fr = lambda r,rho,n: ((n-2)*scipy.special.gamma(n-1)*(1-rho**2)**((n-1)/2.0)*(1-r**2)**((n-4)/2.0))/(np.sqrt(2*np.pi)*scipy.special.gamma(n-0.5)*(1-r*rho)**(n-1.5))*scipy.special.hyp2f1(0.5, 0.5, (2*n-1)/2.0, (r*rho+1)/2.0)
    f = lambda r,n: (1-r**2)**((n-4.0)/2.0)/scipy.special.beta(0.5, (n-2)/2.0)
    #N = a.size
    pval = scipy.integrate.quad(lambda r: f(r, N), corr, 1)[0]
    if alt=="greater":
        return pval
    elif alt=="less":
        return 1-pval
    elif alt=="nonzero":
        return min(pval, 1-pval)

def correlation_pvalue(a, b, nboot=1e4, confinterval=0.95, method="pearson"):
    """Computes a bootstrap p-value for the correlation between [a] and [b].
    The alternative hypothesis for this test is that the correlation is zero or less.
    This function randomly resamples the timepoints in the [a] and [b] and computes
    the correlation for each sample.

    Parameters
    ----------
    a : array_like, shape (N,)
    b : array_like, shape (N,)
    nboot : int, optional
        Number of bootstrap samples to compute, default 1e4
    conflevel : float, optional
        Confidence interval size, default 0.95
    method : string, optional
        Type of correlation to use, can be "pearson" (default) or "robust"
        
    Returns
    -------
    bspval : float
        The fraction of bootstrap samples with correlation less than zero.
    bsconf : (float, float)
        The [confinterval]-percent confidence interval according to the bootstrap.
    ppval : float
        The probability that the correlation is zero or less according to parametric
        computation using Fisher transform.
    pconf : (float, float)
        The parametric [confinterval]-percent confidence interval according to
        parametric computation using Fisher transform.
    bootcorrs : array_like, shape(nboot,)
        The correlation for each bootstrap sample
    """
    ocorr = np.corrcoef(a, b)[0,1]
    conflims = ((1-confinterval)/2, confinterval/2+0.5)
    confinds = map(int, (conflims[0]*nboot, conflims[1]*nboot))

    N = len(a)
    inds = make_randinds(N, nboot, algo="bytes")
    rsa = a[inds] ## resampled a
    rsb = b[inds] ## resampled b

    if method=="pearson":
        za = (rsa-rsa.mean(0))/rsa.std(0)
        zb = (rsb-rsb.mean(0))/rsb.std(0)
        bootcorrs = np.sum(za*zb, 0)/(N-1) ## The correlation between each pair
    elif method=="robust":
        bootcorrs = np.array([robust_correlation(x,y)[0] for (x,y) in zip(rsa.T, rsb.T)])
    else:
        raise ValueError("Unknown method: %s"%method)

    ## Compute the bootstrap p-value
    #bspval = np.mean(bootcorrs<0) ## Fraction of correlations smaller than zero
    bspval = np.mean(bootcorrs>ocorr)
    bsconf = (np.sort(bootcorrs)[confinds[0]], np.sort(bootcorrs)[confinds[1]])

    ## Compute the parametric bootstrap p-value using Fisher transform
    zccs = np.arctanh(bootcorrs)
    ppval = scipy.stats.norm.cdf(0, loc=zccs.mean(), scale=zccs.std())
    pconf = tuple(map(lambda c: np.tanh(scipy.stats.norm.isf(1-c, loc=zccs.mean(), scale=zccs.std())), conflims))

    ## return things!
    return bspval, bsconf, ppval, pconf, bootcorrs

def rolled_correlation_pvalue(a, b, nboots=1e4):
    N = len(a)
    bscorrs = []
    for bi in range(nboots):
        bscorrs.append(np.corrcoef(np.roll(a, np.random.randint(0,N)), b)[0,1])

    ocorr = np.corrcoef(a, b)[0,1]
    bspval = np.mean(np.array(bscorrs)>ocorr)
    return bspval

def shuffled_correlation_pvalue(a, b, nboots=1e4, blocklen=263):
    N = len(a)
    bscorrs = []
    for bi in range(nboots):
        splits = np.split(a, np.arange(0,N,blocklen))
        np.random.shuffle(splits)
        shuffa = np.hstack(splits)
        bscorrs.append(np.corrcoef(shuffa, b)[0,1])

    ocorr = np.corrcoef(a, b)[0,1]
    bspval = np.mean(np.array(bscorrs)>ocorr)
    print "correlation: %0.05f. p-value: %0.05f" % (ocorr, bspval)
    return bspval


def robust_correlation(a, b, cutoff=2.5):
    """Computes a robust estimate of the correlation between [a] and [b] using the
    least mean squares (LMS) based method defined in: 
    Abdullah, 1990, "On a robust correlation coefficient"

    First, outliers are removed based on the residual of the linear regression of [a]
    on [b] and the [cutoff]. Then the correlation is computed on non-outliers.

    Parameters
    ----------
    a : array_like, shape (N,)
    b : array_like, shape (N,)
    cutoff : float, default=2.5
        The cutoff for outlier detection.
    
    Returns
    -------
    goodcorr : float
        The correlation between non-outliers in [a] and [b]
    """
    assert a.size == b.size
    rho = np.corrcoef(a, b)[0,1]
    zscore = lambda v: (v-v.mean())/v.std()
    res = b - (a*np.linalg.lstsq(np.atleast_2d(zscore(a)).T, b)[0])
    s = 1.4826*(1+5/(a.size-rho))*np.sqrt(np.median(res**2))
    goodpts = np.abs(res/s)<cutoff
    goodcorr = np.corrcoef(a[goodpts], b[goodpts])[0,1]

    return goodcorr, goodpts