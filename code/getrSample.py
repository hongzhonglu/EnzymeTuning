def getrSample(mu, sigma, lb, ub, step, method):
    r = []
    for i in range(len(lb)):
        if lb[i] == ub[i] and lb[i] == 0:
            lb[i] = -2
            ub[i] = 9
        if method == 'normal':
            mutmp = np.log10(np.array(mu)/3600)
            sigmatmp = sigma
            pdd = np.random.normal(loc=mutmp, scale=sigmatmp)
            r_tmp = np.random.normal(loc=1, scale=step)
            if r_tmp < lb[i]:
                r_tmp = lb[i]
            if r_tmp > ub[i]:
                r_tmp = ub[i]
            r.append(r_tmp)
        elif method == 'Uniform':
            mutmp = np.log10(np.array(mu)/3600)
            sigmatmp = sigma
            pdd = np.random.uniform(low=mutmp - sigmatmp, high=mutmp + sigmatmp)
            t = np.random.uniform(low=-2, high=8)
            r_tmp = np.random.normal(loc=1, scale=step)
            r.append(r_tmp)
    return r