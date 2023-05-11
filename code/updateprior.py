def updateprior(x):
    from distfit import distfit
    x = np.log10(np.array(x)/3600)
    dist = distfit(distr = 'norm')
    a = []
    b = []
    for i in range(len(x[:,0])):
        results = dist.fit_transform(x[i,:])
        mu = dist.model['loc']
        a_tmp = 10**mu*3600
        sigma = dist.model['scale']
        a.append(a_tmp)
        b.append(sigma)
    a = np.array(a)
    b = np.array(b)
    return a,b