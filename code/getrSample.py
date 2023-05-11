def getrSample(mu, sigma, lb, ub, step, method = 'Uniform'):
    r = []
    for i in range(len(lb)):
        if lb[i] == ub[i] and lb[i] == 0:
            lb[i] = -2
            ub[i] = 9
        if method == 'normal':
            mutmp = np.log10(np.array(mu)/3600)
            sigmatmp = sigma
            pdd = np.random.normal(loc=mutmp, scale=sigmatmp)
            r_tmp = np.random.choice(pdd, size = (1,int(step)))
            if np.any(r_tmp) < lb[i]:
                r_tmp[r_tmp<lb[i]] = lb[i]
            if np.any(r_tmp) > ub[i]:
                r_tmp[r_tmp>ub[i]] = ub[i]
            r.append(r_tmp)
        elif method == 'Uniform':
            mutmp = np.log10(np.array(mu)/3600)
            sigmatmp = sigma
            pdd = np.random.uniform(low=mutmp - sigmatmp, high=mutmp + sigmatmp)
            t = np.clip(pdd, -2, 8)
            #t = np.random.uniform(pd,low=-2, high=8)
            r_tmp = np.random.choice(t, size = (1,int(step)))
            r.append(r_tmp)
    r = np.array(r)[:,0]
    return r