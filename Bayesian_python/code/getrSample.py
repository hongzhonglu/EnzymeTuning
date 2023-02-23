def getrSample(mu, sigma, lb, ub, step, method):
    if lb == ub and lb == 0:
        lb = -2
        ub = 9
    if method == 'normal':
        mutmp = math.log(mu / 3600, 10)
        sigmatmp = sigma
        pdd = np.random.normal(loc=mutmp, scale=sigmatmp)
        r = np.random.normal(loc=1, scale=step)
        if r < lb:
            r = lb
        if r > ub:
            r = ub
        r = 10 ** r * 3600
    elif method == 'uniform':
        mutmp = math.log(mu / 3600, 10)
        sigmatmp = sigma
        pdd = np.random.uniform(low=mutmp - sigmatmp, high=mutmp + sigmatmp)
        t = np.random, uniform(low=-2, high=8)
        r = np.random.normal(loc=1, scale=step)
        r = 10 ** r * 3600
