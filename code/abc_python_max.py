def abc_python_max(model, enzymedata, kcat_random_all, tot_prot_weight, growthdata, max_growth, proc, sample_generation,
                   j, rxn2block):
    nstep = sample_generation / proc
    nstep = int(nstep)
    rmse_final = np.zeros((1, nstep))
    kcat_sample = kcat_random_all[:, ((j - 1) * nstep ):j * nstep]

    # get carbonnum for each exchange rxn to further calculation of error
    if len(model["excarbon"]) == 0:  ##检查model中是否有excarbon, ~is_field
        model = addCarbonNum(model)  ##补写function

    for k in range(nstep):
        print('nstep:' + str(k + 1) + '/' + str(nstep))
        kcat_random = kcat_sample[:, k]
        prot_cost_info_value = [enzymedata["MW"][i] / kcat_random[i] for i in range(len(kcat_random))]
        prot_cost_info_id = enzymedata["rxn_list"]

        # first search with substrate constrain
        objective = 'r_2111'
        osenseStr = 'max'
        if len(growthdata) != 0:
            result_1 = rmsecal(model, growthdata, constrain= True, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block)
            result_1 = np.array(result_1)
            rmse_1 = result_1[0]
            exp_1 = result_1[1]
            simulated_1 = result_1[2]
        else:
            rmse_1 = []
            exp_1 = []
            simulated_1 = []

        # second searcch for maxmial growth rate without constrain
        if len(max_growth) != 0:
            result_2 = rmsecal(model, max_growth, constrain= False, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block)
            result_2 = np.array(result_2)
            rmse_2 = result_2[0]
            exp_2 = result_2[1]
            simulated_2 = result_2[2]
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

def rmsecal(model, data, constrain, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block):
    data = np.array(data)
    exp_d = []
    sim = []
    sol = []
    rmse_tmp = []
    simulated = np.zeros((len(data[:, 0]), 9))
    rxnNames = model['rxnNames']
    exp = np.array(data[:, 2:11])
    exp = exp * [1, -1, 1, 1, 1, 1, 1, 1, -1]
    exp = exp.astype(np.float)

    for i in range(len(data[0:10, 0])):
        ex_mets = ['growth', data[i, 1][0] + " exchange", 'acetate exchange', 'ethanol exchange',
                   'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
                   'oxygen exchange']
        idx = []
        temp = []
        for k in range(len(rxnNames)):
            if rxnNames[k] in ex_mets:
                temp.append(k)
        idx.append(temp)
        model_tmp = model

        model_tmp = changeMedia(model_tmp, 'D-glucose', data[i, 15])

        id_1634 = np.where(model_tmp['rxns']=='r_1634')
        model_tmp['ub'][id_1634[0][0]] = 0
        model_tmp['lb'][id_1634[0][0]] = 0
        id_1631 = np.where(model_tmp['rxns']=='r_1631')
        model_tmp['ub'][id_1631[0][0]] = 0
        model_tmp['lb'][id_1631[0][0]] = 0

        if data[i, 13]=="anaerobic" or data[i, 13]=="limited":
            model_tmp = anaerobicModel(model_tmp)

        if data[i, 13]=="limited" :
            id_EX_O2 = np.where(model_tmp['rxnNames']=='oxygen exchange')
            model_tmp['lb'][id_1631[0][0]] = -5

        id_1714 = np.where(model_tmp['rxns']=='r_1714')
        model_tmp['lb'][id_1714[0][0]] = 0

        id_substrate = np.where(model_tmp['rxnNames']==(data[i, 1][0] + " exchange"))
        if constrain == False:
            model_tmp['lb'][id_substrate[0][0]] = -1000
        else:
            model_tmp['lb'][id_substrate[0][0]] = exp[:, 1][i]

        io.savemat('modeltemp.mat', {'model': model_tmp})
        # save_matlab_model(model_tmp,"modeltemp.mat")
        modeltemp = load_matlab_model("modeltemp.mat")
        sol_tmp = modeltemp.optimize()
        print("No. " + str(i + 1) + " result is " + str(sol_tmp.objective_value))
        sol = sol_tmp.fluxes  # sol[:,i] = sol_tmp.fluxes

        tmp = np.where(~np.isnan(exp[i]))[0]
        excarbon = model_tmp['excarbon'][:, idx[0]]
        for x in range(len(idx[0])):
            if excarbon[:, x]==0:
                excarbon[:, x] = 1
        exp_tmp = []
        for s in range(len(tmp)):
            exp_tmp.append(exp[0, s] * excarbon[:, s])
        sol_idx = np.transpose(sol[idx[0]])
        simulated_tmp = []
        simulated_tmp.append(np.array(sol_idx)*excarbon[:,tmp])  # normalize the growth rate issue by factor 10

        exp_block = np.zeros((1, 218))  # 218怎么替换
        rxnblockidx_pre = np.setdiff1d(rxn2block, model_tmp['rxns'][idx[0][1]])
        rxnblockidx = []
        for k in range(len(model_tmp['rxns'])):
            if model_tmp['rxns'][k] in rxnblockidx_pre:
                rxnblockidx.append(k)
        simulated_block = []
        soll = np.array(sol)
        for x in range(len(rxnblockidx)):
            simulated_block.append(model_tmp['excarbon'][:, rxnblockidx[x]] * soll[rxnblockidx[x]])
        id_zero = np.where(np.array(simulated_block)!=0)
        exp_block = exp_block[:, id_zero[0]]
        simulated_block = np.array(simulated_block)[id_zero[0]]
        exp_tmp = np.array(exp_tmp)

        if constrain:
            rmse_tmp.append(np.sqrt(mean_squared_error(np.append(exp_tmp,np.transpose(exp_block)),np.append(simulated_tmp,np.transpose(simulated_block)))))
        else:
            if len(exp_tmp) >= 2:
                rmse_tmp.append(np.sqrt(mean_squared_error(exp_tmp[0:2], [simulated_tmp[0][0:2]])))
            else:
                rmse_tmp.append(np.sqrt(mean_squared_error(exp_tmp[0],[simulated_tmp[0][8]])))

        simulated[i, :] = np.transpose(sol[idx, i])

    rmse = sum(rmse_tmp) / len(data[:, 0])

    return rmse, exp, simulated

