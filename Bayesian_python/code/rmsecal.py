def rmsecal(model, data, constrain, objective, osenseStr, prot_cost_info, tot_prot_weight, rxn2block):
    data = np.array(data)
    exp_d = []
    sim = []
    simulated = np.zeros(len(data[:, 0]), 9)

    for i in range(len(data[:, 0])):
        exp = np.array(data[:, 2:11])
        exp = exp * [1, -1, 1, 1, 1, 1, 1, 1, -1]
        ex_mets = ['growth', [data[i, 1], ' exchange'], 'acetate exchange', 'ethanol exchange', 'glycerol exchange',
                   'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange', 'oxygen exchange']
        rxnNames = str(model['rxnNames'])
        idx = rxnNames.find(ex_mets)  ##[Lia,Locb] = ismember(A,B)
        model_temp = model

        model_tmp = changeMedia(model_tmp, 'D-glucose', data[i, 15])
        model_tmp = changeRxnBounds(model_tmp, 'r_1634', 0, 'b')
        model_tmp = changeRxnBounds(model_tmp, 'r_1631', 0, 'b')

        if data[i, 13] == "anaerobic" or data[i, 13] == "limited":
            model_tmp = anaerobicModel(model_tmp)
        if data[i, 13] == "limited":
            if str(model_tmp["rxnNames"]) == 'oxygen exchange':  # if model_tmp["rxnNames"] == 'oxygen exchange'
                model_tmp["lb"] = -5

        for x in range(len(model_tmp["rxns"])):
            if not constrain:
                if model_tmp["rxns"][x] == 'r_1714':
                    model_tmp["lb"][x] = 0
                elif model_tmp["rxns"][x] == model_tmp["rxns"][idx[1]]:
                    model_tmp["lb"][x] = -1000
            else:
                if model_tmp["rxns"][x] == 'r_1714':
                    model_tmp["lb"][x] = 0
                    model_tmp["lb"][idx[1]] = exp[i, 1]

        sol_tmp = solveModel(model_tmp, objective, osenseStr, prot_cost_info, tot_prot_weight, 'ibm_cplex')  # gurobi
        sol[:, i] = sol_tmp["x"]

        if not np.isnan(exp[i]) and not np.isinf(exp[i]):
            tmp = 1
        else:
            tmp = 0

        excarbon = model["excarbon"]["idx"]
        if excarbon == 0:
            excarbon = 1
        exp_tmp = exp[i, tmp] * excarbon[tmp]
        sim_tmp = np.transpose(sol[idx["tmp"], i])
        simulated_tmp = sim_tmp * excarbon[tmp]  # normalize the growth rate issue by factor 10

        setdiff = set(rxn2block).difference(set(model_tmp["rxns"][idx[1]]))
        exp_block = np.zeros((1, len(setdiff)))
        rxnblockidx = model_tmp["rxns"].find(setdiff)
        sim_block = np.transpose(sol[rxnblockidx, i])
        simulated_block = sim_block * model['excarbon'][rxnblockidx]

        sim_index = np.nonzero(simulated_block)
        sim_index = np.mat(sim_index)
        exp_block_new = []
        for i in range(len(exp_block)):
            if i in sim_index:
                exp_block_new.append(exp_block[i])
        exp_block = exp_block_new

        simulated_block_new = []
        for i in range(len(simulated_block)):
            if i in sim_index:
                simulated_block_new.append(simulated_block[i])
        simulated_block = simulated_block_new

        if constrain:
            exp0 = exp_tmp + exp_block
            simulated0 = simulated_tmp + simulated_block
            rmse_tmp[i] = metrics.mean_squared_error(exp0, simulated0) ** 0.5
        else:
            if len(exp_tmp) >= 2:
                rmse_tmp[i] = metrics.mean_squared_error(exp_tmp[0:2], simulated_tmp[0:2]) ** 0.5
            else:
                rmse_tmp[i] = metrics.mean_squared_error(exp_tmp[0], simulated_tmp[0]) ** 0.5

        simulated[i, :] = np.transpose(sol[idx, i])

    rmse = sum(rmse_tmp) / len(data[:, 0])

    return rmse, exp, simulated