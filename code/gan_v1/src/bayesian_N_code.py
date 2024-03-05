import time
import numpy as np
import os
import cobra
from pathlib import Path
from sklearn.metrics import mean_squared_error
import pandas as pd
from cobra.io import load_matlab_model , save_matlab_model , read_sbml_model , write_sbml_model
import scipy.io as scio
from scipy import sparse
import cobra.util
from cobra import Model , Reaction , Metabolite
from multiprocessing import Pool , cpu_count

# from code0824 import abc_python_max
# from code0824 import parallel_task
# from Bayesian_generation_cobra import anaerobicModel
# from Bayesian_generation_cobra import changeGAM
from src.ec_N_lim_sen import changeMedia
from code_generation_preparation import getrSample
from code_generation_preparation import sumBiomass
from code_generation_preparation import updateprior
def rmsecal(model , model_cobra , data ,constrain ,objective , osenseStr ,prot_cost_info_id , prot_cost_info_value ,prot_cost_info, rxn2block):
    #data = np.array (data)
    rmse_tmp = []
    simulated = np.zeros ((len (data[: , 0]) , 13))
    rxnNames = []
    for rxn in model_cobra.reactions:
        rxnNames.append(rxn.name)
    rxnNames = np.array(rxnNames)
    for i in range (len (data[: , 0])):
        exp = np.array (data)  # u sub ace eth gly pyr ethyl_acetate co2 o2
        exp = exp * [1 , -1 , 1 , 1 , 1 , 1 , 1 , 1 , -1 ,1, -1 , -1 , -1 , -1]
        exp = exp.astype (float)
        ex_mets = ['growth' , "D-glucose exchange" , 'ethanol exchange' ,'acetate exchange' ,
                    'pyruvate exchange' , 'succinate exchange' ,'glycerol exchange' , 'carbon dioxide exchange' ,
                   'oxygen exchange', 'ammonium exchange', 'L-glutamine exchange',
                   'L-phenylalanine exchange', 'L-isoleucine exchange']

        idx = []
        for k in range (len (ex_mets)):
            temp = np.where (rxnNames==ex_mets[k])
            idx.append (temp)
        idx = np.array (idx)
        idx = np.transpose (idx[: , 0])
        # Create a temp model
        start_time1 = time.time ()

        with model_cobra as model_tmp:

            # Suitable for different carbon sources
            #model_tmp = changeMedia (model , model_tmp , data[i , 1][0] , data[i , 15])  # 'D-glucose'
            model_tmp = changeMedia (model , model_tmp , 'D-glucose' , "MIN")
            model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
            model_tmp.reactions.get_by_id ('r_1631').bounds = 0 , 0

            #model_tmp.reactions.get_by_id ('r_1714').lower_bound = 0
            model_tmp.reactions.get_by_id ('r_1761').upper_bound = 1000

            #if constrain==False:
            #    model_tmp.reactions[idx[0][1]].lower_bound = -1000
            #else:
            #    model_tmp.reactions[idx[0][1]].lower_bound = exp[i , 1]
            model_tmp.reactions[idx[0][1]].lower_bound = exp[i , 1] #D-glucose exchange
            model_tmp.reactions.get_by_id("r_1654").bounds = exp[i,-4],0 #NH4 exchange
            #model_tmp.reactions.get_by_id("r_1654").lower_bound = -1000
            model_tmp.reactions.get_by_id("r_1851").bounds = exp[i,-3],0 #L-glutamine exchange
            #model_tmp.reactions.get_by_id("r_1851").lower_bound = -1000
            model_tmp.reactions.get_by_id("r_1903").bounds = exp[i,-2],0 #L-phenylalanine exchange
            #model_tmp.reactions.get_by_id("r_1903").lower_bound = -1000
            model_tmp.reactions.get_by_id("r_1897").bounds = exp[i,-1],0 #L-isoleucine exchange
            #model_tmp.reactions.get_by_id("r_1897").lower_bound = -1000
            #model_tmp.reactions.get_by_id("prot_pool_exchange").upper_bound = data[i,9]*1000 #prot_pool

            tot_prot_weight = data[i,9]
            model_tmp.reactions.get_by_id("SK_cost").upper_bound = tot_prot_weight*1000
            nMets = len (model_tmp.metabolites)
            nRxns = len (model_tmp.reactions)
            cost_list = np.zeros ((nRxns , 1))

            for p , rxnid in enumerate (model['rxns']):
                rxnid = tuple (rxnid)
                cost = prot_cost_info.get (rxnid[0][0] , 0)
                cost_list[p , 0] = cost

            cost_floats = [float (m) for m in cost_list]



            met = model_tmp.metabolites.get_by_id ("cost")
            filtered_reactions = [reaction for reaction in model_tmp.reactions if met in reaction.metabolites]
            for q in range (len (filtered_reactions)-1):
                cost_value = prot_cost_info_value[q]
                coeff = filtered_reactions[q].get_coefficient ('cost')
                filtered_reactions[q].add_metabolites ({"cost": cost_value - coeff})

            # Solve the temp model
            sol_tmp = model_tmp.optimize ()
            sol = sol_tmp.fluxes  # sol[:,i] = sol_tmp.fluxes

        print ("No. " + str (i + 1) + "finish !")

        tmp = np.where (~np.isnan (exp[i]))[0]
        tmp = tmp[tmp != 9]
        tmp = np.where(tmp > 9, tmp - 1, tmp)

        # Normalize the number of carbon
        excarbon = model['excarbon'][: , idx[0]]
        for x in range (len (idx[0 , :])):
            if excarbon[: , x]==0:
                excarbon[: , x] = 1
        exp_tmp = []
        for s in range (len (tmp)):
            #exp_tmp.append (exp[i , tmp[s]] * excarbon[: , tmp[s]])
            exp_tmp.append (exp[i , tmp[s]])
        sol_idx = np.transpose (sol[idx[0]])
        sol_idx = np.array (sol_idx)
        simulated_tmp = []
        for t in range (len (tmp)):
            simulated_tmp.append (
                np.array (sol_idx)[tmp[t]] * excarbon[: , tmp[t]])  # normalize the growth rate issue by factor 10

        # The expected blocked reactions
        exp_block = np.zeros ((1 , 218))  # 218怎么替换
        rxnblockidx_pre = np.setdiff1d (rxn2block , model['rxns'][idx[0][1]])
        rxnblockidx = []
        for k in range (len (model['rxns'])):
            if model['rxns'][k] in rxnblockidx_pre:
                rxnblockidx.append (k)

        simulated_block = []
        for t in range (len (rxnblockidx)):
            simulated_block.append (np.array (sol)[rxnblockidx[t]] * model['excarbon'][: , rxnblockidx[t]])
        id_zero = np.where (np.array (simulated_block)!=0)
        exp_block = exp_block[: , id_zero[0]]
        simulated_block = np.array (simulated_block)[id_zero[0]]
        exp_tmp = np.array (exp_tmp)


        # Calculate the RMSE
        if constrain:
            rmse_tmp.append (np.sqrt (mean_squared_error (np.append (exp_tmp , np.transpose (exp_block)) ,
                                                          np.append (simulated_tmp , np.transpose (simulated_block)))))
        else:
            if len (exp_tmp) >= 2:
                rmse_tmp.append (np.sqrt (mean_squared_error (exp_tmp[0:2] , np.array (simulated_tmp)[0:2])))
            else:
                rmse_tmp.append (np.sqrt (mean_squared_error (exp_tmp[0] , [simulated_tmp[0]])))
        simulated[i , :] = np.transpose (np.array (sol)[idx[0]])
        #print (rmse_tmp)
        #print (simulated[i , :])
    rmse = sum (rmse_tmp) / len (data[: , 0])
    print ("RMSE = " + str (rmse))

    # end_time = time.time ()
    # execution_time = end_time - start_time
    # print ("execution time: " , execution_time , " seconds")

    # return rmse, exp, simulated
    return rmse , exp , simulated
def abc_python_N(model , model_cobra , enzymedata , kcat_random_all , growthdata , proc , sample_generation , j , rxn2block , num):
    nstep = sample_generation / proc
    nstep = int (nstep)
    rmse_final = np.zeros ((1 , nstep))
    if j == 1:
        kcat_sample = kcat_random_all
    else:
        kcat_sample = kcat_random_all[: , ((j - 1) * nstep):(j * nstep)]

    # get carbonnum for each exchange rxn to further calculation of error
    if len (model["excarbon"])==0:  ##检查model中是否有excarbon, ~is_field
        model = addCarbonNum (model)  ##补写function
    for k in range (nstep):
        print ('nstep:' + str (k + 1) + '/' + str (nstep))
        if nstep == 1 :
            kcat_random = kcat_sample
        else:
            kcat_random = kcat_sample[: , k]
        #prot_cost_info_value = [enzymedata["MW"][i] / kcat_random[i] for i in range (len (kcat_random))]
        #prot_cost_info_id = enzymedata["rxn_list"]
        prot_cost_info_value = [enzymedata["MW"][i] / kcat_random[i] for i in range (len (kcat_random))]
        prot_cost_info_value = [item[0] for item in prot_cost_info_value]
        prot_cost_info_id = [str(item[0][0]) for item in enzymedata["rxn_list"]]
        prot_cost_info = dict (zip (prot_cost_info_id , prot_cost_info_value))


        # first search with substrate constrain
        objective = 'r_2111'
        osenseStr = 'max'
        start_time = time.time ()

        if len (growthdata)!=0:
            rmse , exp , simulated = rmsecal (model , model_cobra , growthdata , True , objective , osenseStr ,
                                                    prot_cost_info_id , prot_cost_info_value ,prot_cost_info, rxn2block)
        else:
            rmse = []
            exp = []
            simulated = []
        end_time = time.time ()
        execution_time = end_time - start_time
        print ("execution time: " , execution_time , " seconds")

        rmse_final[0 , k] = np.nanmean (rmse , 0)

        # only output simulated result for one generation
        if nstep!=1 or sample_generation!=1:
            simulated = []
            exp = []
        print ("rmse_final is " , rmse_final)

    return rmse_final , exp , simulated

def parallel_task(i, output, model, model_cobra, enzymedata, kcat_random_all, growthdata, proc, sample_generation, rxn2block):
    print(i+1)
    rmse_final = abc_python_N(model, model_cobra, enzymedata, kcat_random_all, growthdata, int(proc), int(sample_generation), i+1, rxn2block, output)
    print(rmse_final)
    new_tmp_i = np.zeros((1, 7))
    new_tmp_i[0, :] = rmse_final[0][0]
    path = 'result_bayesian/output_rmse'
    os.makedirs(path, exist_ok=True)
    filename = 'output_'+str(output)+'_rmse'+str(i+1)+'.txt'
    full_path = os.path.join(path, filename)
    np.savetxt(full_path, new_tmp_i)
    return new_tmp_i