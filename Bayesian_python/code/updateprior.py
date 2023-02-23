def updateprior(x):
    import math
    from distfit import distfit
    dist = distfit(distr = 'norm')
    dist.fit_transform(x)
    pd = dist(math.log(x/3600,10))
    mu = dist.model['loc']
    sigma = dist.model['scale']
    a = 10**mu*3600
    b = sigma