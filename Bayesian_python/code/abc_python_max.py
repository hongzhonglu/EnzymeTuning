def abc_python_max(model, enzymedata, kcat_random_all, tot_prot_weight, growthdata, max_growth, proc, sample_generation,
                   j, rxn2block):
    nstep = sammple_generation / proc
    rmse_final = np.zeros((1, nstep))
    kcat_sample = kcat_random_all[:, ((j - 1) * nstep + 1):j * nstep]

    # get carbonnum for each exchange rxn to further calculation of error
    if not len(model["excarbon"]):  ##检查model中是否有excarbon, ~is_field
        model = addCarbonNum(model)  ##补写function

    for k in range(nstep):
        print('nstep:' + str(k) + '/' + str(nstep))
        kcat_random = kcat_sample[:, k]
        prot_cost_info['value'] = enzymedata["MW"] / kcat_random
        prot_cost_info['id'] = enzymedata["rxn_list"]

    # first search with substrate constrain
    objective = 'r_2111'
    osenseStr = 'max'
    if len(growthdata) != 0:
        rmse_1 = rmsecal(model, growthdata, true, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block)[0]
        exp_1 = rmsecal(model, growthdata, true, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block)[1]
        simulated_1 = \
        rmsecal(model, growthdata, true, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block)[2]
    else:
        rmse_1 = []
        exp_1 = []
        simulated_1 = []

    # second searcch for maxmial growth rate without constrain
    if len(max_growth) != 0:
        rmse_2 = rmsecal(model, growthdata, true, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block)[0]
        exp_2 = rmsecal(model, growthdata, true, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block)[1]
        simulated_2 = \
        rmsecal(model, growthdata, true, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block)[2]
    else:
        rmse_2 = []
        exp_2 = []
        simulated_2 = []

    exp = np.array([exp_1, exp_2])
    simulated = np.array([simulated_1, simulated_2])
    rmse = np.array([rmse_1, rmse2])
    rmse_final[0, k] = np.nanmean(rmse, 0)

    # only output simulated result for one generation
    if nstep != 1 or sample_generation != 1:
        simulated = []
        exp = []

    return rmse_final, exp, simulated, growthdata, max_growth