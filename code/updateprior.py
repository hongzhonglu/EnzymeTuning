def updateprior(x):
    import math
    from distfit import distfit
    x = x.reshape(-1)
    x = np.log10(np.array(x)/3600)
    dist = distfit(distr = 'norm')
    results = dist.fit_transform(x)
    mu = dist.model['loc']
    sigma = dist.model['scale']
    a = 10**mu*3600
    b = sigma

    return a,b