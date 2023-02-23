def solveModel(model, objective, osenseStr, prot_cost_info, tot_prot_weight, solver, rxnlst, factor):
    # Determine objective
    for i in range(len(model['c'])):
        model['c'][i] = 1 if model['rxns'][i] == objective
    changerxn = True

    # Add protein cost infomation
    nMets = np.array(model['S']).shape[0]
    nRxns = np.array(model['S']).shape[1]
    cost_list = np.zeros((1, nRxns))
    for i in range(nRxns):
        rxnid = model['rxns'][i]
        if changerxn and len(set(rxnlst) & set(rxnid)) != 0:
            id_tmp = rxnid
            if len(set(prot_cost_info[id]) & set(id_tmp)) != 0:
                cost_prot = np.argwhere(prot_cost_info[id] == id_tmp)
                cost = prot_cost_info['value']['cost_prot'] / factor
            else:
                cost = 0
        else:
            id_tmp = rxnid
            if len(set(prot_cost_info[id]) & set(id_tmp)) != 0:
                cost_prot = np.argwhere(prot_cost_info[id] == id_tmp)  # index怎么修改
                cost = prot_cost_info['value']['cost_prot']
            else:
                cost = 0

        cost_list[0, i] = cost

    # choose solver gurobi/CPLEX


#if solver = 'gurobi':
    #        param.DisplayInterval = 1;
    #        param.FeasibilityTol = 1.0000e-06;
    #        param.OptimalityTol = 1.0000e-06;
    #        param.OutputFlag = 0;
    #        % set constraints
    #        gurobiLP.A = [model.S;cost_list];
    #        gurobiLP.lb = model.lb;
    #        gurobiLP.ub = model.ub;
    #        gurobiLP.rhs = [zeros(nMets,1);tot_prot_weight*1000];
    #        gurobiLP.sense = [repmat('=',nMets,1);'<'];
    #        gurobiLP.modelsense = osenseStr;
    #        gurobiLP.obj = model.c;
    #        % call the solver
    #        resultgurobi = gurobi(gurobiLP,param);
    #        sol = struct();
    #        sol.exitflag = resultgurobi.status;
    #        if strcmp(resultgurobi.status,'OPTIMAL')
    #            sol.x = resultgurobi.x;
    #            sol.obj = sol.x(logical(model.c));
    #            sol.protUsage = cost_list * sol.x / 1000;
    #        else
    #            sol.obj = [];
    #            sol.protUsage = [];
    #            sol.x = zeros(length(model.rxns),1);
    #        end

    return sol