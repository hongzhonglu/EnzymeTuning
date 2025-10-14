from cobra.io import load_matlab_model, read_sbml_model
from cobra import Reaction, Metabolite
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import scipy.io as scio
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
import time
from src.constrain_ecmodel import *
from src.code_generation_preparation import *

def get_r_max(model,model_dl):
    model.objective = model_dl.objective
    model.objective_direction = 'max'
    return model.optimize()

def rmsecal(model , model_cobra , data , constrain , rxn2block):
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
            model_tmp.reactions.get_by_id("prot_pool_exchange").upper_bound = data[i,9]*1000 #prot_pool

            # Solve the temp model
            sol_tmp = model_tmp.optimize ()
            sol = sol_tmp.fluxes  # sol[:,i] = sol_tmp.fluxes

        #print ("No. " + str (i + 1) + "finish !")

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
        exp_block = np.zeros ((1 , 218))
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

def abc_python_max(model , model_cobra , tmp ,proc , sample_generation ,rxn2block):
    nstep = sample_generation / proc
    nstep = int (nstep)
    rmse_final = np.zeros ((1 , nstep))

    # get carbonnum for each exchange rxn to further calculation of error
    if len (model["excarbon"])==0:  ##检查model中是否有excarbon, ~is_field
        model = addCarbonNum (model)  ##补写function
    for k in range (nstep):
        #print ('nstep:' + str (k + 1) + '/' + str (nstep))

        # first search with substrate constrain
        objective = 'r_2111'
        osenseStr = 'max'
        start_time = time.time ()
        if len (tmp)!=0:
            rmse_1 , exp_1 , simulated_1 = rmsecal (model , model_cobra , tmp , True , rxn2block)
        else:
            rmse_1 = []
            exp_1 = []
            simulated_1 = []
        # print("RMSE_1 finish !")
        # second searcch for maxmial growth rate without constrain
        end_time = time.time ()
        execution_time = end_time - start_time
        #print ("execution time: " , execution_time , " seconds")

        # print("RMSE_2 finish !")

        exp = exp_1
        simulated = simulated_1
        rmse = rmse_1
        rmse_final[0 , k] = np.nanmean (rmse , 0)

        # only output simulated result for one generation
        if nstep!=1 or sample_generation!=1:
            simulated = []
            exp = []
        #print ("rmse_final is " , rmse_final)

    return rmse_final , exp , simulated

def changeMedia(model , model_cobra , c_source , media , anox=False , flux=-1000):
    '''
    Function that modifies the ecModel and makes it suitable for batch growth
    simulations on different carbon sources.

    media:  Media type ('YEP' for complex,
                      'MAA' minimal with Aminoacids,
                      'MIN' for minimal media,
                      'MIN+His' for minimal media with his
                      'MIN+Arg' for minimal media with arg
                      'MIN+Citrate' for minimal media with Citrate)
    anox:   (optional) TRUE if anaerobic conditions are desired, DEFAULT= FALSE
    flux:   (Optional) A cell array with measured uptake fluxes in mmol/gDwh

    '''
    #Give the carbon source (c_source) input variable with the following format:
    #c_source = 'D-glucose exchange'
    c_source = c_source + ' exchange'

    # first block any uptake
    exchangerxn = getExchangeRxns (model , reactionType="both")
    for a in range (len (exchangerxn[0])):
        model_cobra.reactions[exchangerxn[0][a]].lower_bound = 0
    pos = getComponentIndexes (model , c_source)

    # The media will define which rxns to fix:
    if media=='YEP':
        N = 25  # Aminoacids + Nucleotides
    elif media=='MAA':
        N = 21  # Aminoacids
    elif media=='MIN':
        N = 1  # Only the carbon source
    elif media=='MIN+His':
        N = 1  # Only the carbon source
        model_cobra.reactions.get_by_id ('r_1893').lower_bound = -0.08  # Histidine exchange
    elif media=='MIN+Arg':
        N = 1  # Only the carbon source
        model_cobra.reactions.get_by_id ('r_1879').lower_bound = -0.08  # L-arginine exchange
    elif media=='MIN+Citrate':
        N = 1  # Only the carbon source
        model_cobra.reactions.get_by_id ('r_1687').lower_bound = -0.08  # citrate exchange
    # LB parameter (manually optimized for glucose on Min+AA):
    b = -0.08
    # LB parameter (manually optimized for glucose complex media):
    c = -2
    flux = np.array([flux])
    # Define fluxes in case of ec model:
    if N > 1:
        flux = np.append(flux, b * np.ones(N))
        if N > 21:
            flux = np.append(flux, c * np.ones(N - 20))

    # Fix values as LBs:
    i = 0
    for i in range (N):
        model_cobra.reactions[pos[i][0]].lower_bound = flux[i]

    # Allow uptake of essential components
    model_cobra.reactions.get_by_id ('r_1654').lower_bound = -1000  # 'ammonium exchange';
    model_cobra.reactions.get_by_id ('r_2100').lower_bound = -1000  # 'water exchange';
    model_cobra.reactions.get_by_id ('r_1861').lower_bound = -1000  # 'iron(2+) exchange';
    model_cobra.reactions.get_by_id ('r_1992').lower_bound = -1000  # 'oxygen exchange';
    model_cobra.reactions.get_by_id ('r_2005').lower_bound = -1000  # 'phosphate exchange';
    model_cobra.reactions.get_by_id ('r_2060').lower_bound = -1000  # 'sulphate exchange';
    model_cobra.reactions.get_by_id ('r_1832').lower_bound = -1000  # 'H+ exchange';
    model_cobra.reactions.get_by_id ('r_4593').lower_bound = -1000  # 'chloride exchange';
    model_cobra.reactions.get_by_id ('r_4595').lower_bound = -1000  # 'Mn(2+) exchange';
    model_cobra.reactions.get_by_id ('r_4596').lower_bound = -1000  # 'Zn(2+) exchange';
    model_cobra.reactions.get_by_id ('r_4597').lower_bound = -1000  # 'Mg(2+) exchange';
    model_cobra.reactions.get_by_id ('r_2049').lower_bound = -1000  # 'sodium exchange';
    model_cobra.reactions.get_by_id ('r_4594').lower_bound = -1000  # 'Cu(2+) exchange';
    model_cobra.reactions.get_by_id ('r_4600').lower_bound = -1000  # 'Ca(2+) exchange';
    model_cobra.reactions.get_by_id ('r_2020').lower_bound = -1000  # 'potassium exchange';

    # Block some production fluxes
    model_cobra.reactions.get_by_id ('r_1663').upper_bound = 0  # bicarbonate exchange;
    model_cobra.reactions.get_by_id ('r_4062').upper_bound = 0  # lipid backbone exchange;
    model_cobra.reactions.get_by_id ('r_4064').upper_bound = 0  # lipid chain exchange;

    # Allow biomass production
    model_cobra.reactions.get_by_id ('r_2111').upper_bound = 1000  # growth;

    return model_cobra

def getComponentIndexes(model , c_source):
    # c_source = c_source + ' exchange'
    pos = []
    pos.append (np.where (model["rxnNames"]==c_source))
    pos.append (np.where (model["rxnNames"]=='L-alanine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-arginine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-asparagine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-aspartate exchange'))
    pos.append (np.where (model["rxnNames"]=='L-cysteine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-glutamine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-glutamate exchange'))
    pos.append (np.where (model["rxnNames"]=='L-glycine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-histidine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-isoleucine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-leucine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-lysine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-methionine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-phenylalanine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-proline exchange'))
    pos.append (np.where (model["rxnNames"]=='L-serine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-threonine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-tryptophan exchange'))
    pos.append (np.where (model["rxnNames"]=='L-tyrosine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-valine exchange'))
    pos.append (np.where (model["rxnNames"]=="2'-deoxyadenosine exchange"))
    pos.append (np.where (model["rxnNames"]=="2'-deoxyguanosine exchange"))
    pos.append (np.where (model["rxnNames"]=='thymidine exchange'))
    pos.append (np.where (model["rxnNames"]=='deoxycytidine exchange'))
    pos.append (np.where (model["rxnNames"]=='D-glucose exchange'))
    pos = np.array (pos)
    pos = pos[: , 0]
    return pos

def updateEcGEMkcat(ecModel, target_gene, rxnID, kcat_m):
    """
    The function is used to update the kcat of enzyme in specific reaction from ecModel
    Generally, an enzyme and the related reaction determine the corresponding kcat value.

    :param ecGEM: the enzyme constrainted model
    :param target_gene: the gene id of the enzyme
    :param rxnID: the reaction id contains the enzyme
    :param kcat_m: the new value of kcat
    :return:
    ecModel

    """
    #with ecGEM as ecModel:
    coef = 1 / kcat_m
    ss = ecModel.reactions.get_by_id(rxnID).reaction
        # split as coefficient
    ss1 = ss.split(" + ")
    ss2 = []
    for xx in ss1:
        if target_gene + '[' in xx:
            xx1 = xx.split(' ')[1]
            xx2 = str(coef) + ' ' + xx1
            #print('old coefficient', xx)
            #print('old kcat', 1 / float(xx.split(' ')[0]))
            #print('new coefficient', xx2)
            #print('new kcat', kcat_m)
            ss2.append(xx2)
        else:
            ss2.append(xx)
    rxn_update = " + ".join(ss2)
        #print('old rxn:', ss)
        #print('new rxn:', rxn_update)
    ecModel.reactions.get_by_id(rxnID).reaction = rxn_update
    return ecModel



def parallel_task_preprocessing (j , growthdata , eModel , model , kcat_final , proteomicsdata_Clim, proteomicsdata_Nlim , toy_data) :
    if j % 100==0:
        print (f'current set processed: {j}')

    index = 'kcat_value' + str (j)
    toy_data_tmp = pd.DataFrame (columns=['error' , 'corr' , 'pvalue','number'])
    corr_all = []
    ex_mets = ['biomass pseudoreaction' , 'D-glucose exchange' , 'acetate exchange' , 'ethanol exchange' ,
               'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
               'oxygen exchange' , 'prot_pool_exchange']
    # find the related rxnID
    idx = []
    for name0 in ex_mets:
        # print(name0)
        s = getRxnByReactionName (model=eModel , name=name0)
        if len (s) > 1:
            print ("need check")
        elif len (s)==1:
            idx.append (s[0])
    model_tmp = eModel.copy ()
    for i in range (len (growthdata)) :
        dilutionrate = growthdata[i, 0]
        model_tmp.solver = "cplex"
        for k, row in kcat_final.iterrows () :
            gene = row['Gene']
            mean = row[f'vsyn_value{j}']
            kdeg_x = np.maximum(row['kdeg_x'],row['kdeg_y'])
            kdeg_y = np.minimum(row['kdeg_x'],row['kdeg_y'])
            reaction_id = f'draw_prot_{gene}'
            if reaction_id in model_tmp.reactions :
                #model_tmp.reactions.get_by_id (reaction_id).upper_bound = mean/(dilutionrate + kdeg)
                model_tmp.reactions.get_by_id (reaction_id).bounds = mean/(dilutionrate + kdeg_x),mean/(dilutionrate + kdeg_y)
        # for k in range (len (kcat_final)):
        #    target_gene0 = kcat_final['gene'][k]
        #    kcat_m0 = kcat_final[index][k]
        #    rxn0 = kcat_final['rxnID'][k]
        #    model_tmp = updateEcGEMkcat (ecModel=model_tmp ,target_gene=target_gene0 ,rxnID=rxn0 ,kcat_m=kcat_m0)
        model_tmp = changeMedia (model, model_tmp, 'D-glucose', "MIN")
        # model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
        model_tmp.reactions.get_by_id ('r_1631').bounds = 0, 0
        model_tmp.reactions.get_by_id ('r_1761').upper_bound = growthdata[i, 2]  # Ethanol exchange
        model_tmp.reactions.get_by_id ('r_1714').lower_bound = -growthdata[i, 1]  # D-glucose exchange
        model_tmp.reactions.get_by_id ('r_1672').upper_bound = growthdata[i, 3]  # CO2 exchange
        model_tmp.reactions.get_by_id ('r_1992').lower_bound = -growthdata[i, 4]  # oxygen exchange
        model_tmp.reactions.get_by_id ("prot_pool_exchange").upper_bound = growthdata[i, 5]*1000  # protein
        if i <= 6 :
            model_tmp.reactions.get_by_id ("r_1672").lower_bound = -growthdata[i, 6]  # NH4 exchange
        print (model_tmp.optimize ().objective_value)
        model_tmp.reactions.get_by_id ('r_1714').lower_bound = -1000  # glucose uptake
        model_tmp.reactions.get_by_id (idx[0]).lower_bound = dilutionrate
        model_tmp.objective = {model_tmp.reactions.r_1714 : 1}  # minimize the uptake of glucose
        solution2 = model_tmp.optimize ()
        print (solution2.objective_value)
        # then fix glucose uptake and minimize the protein pool
        model_tmp.reactions.get_by_id ('r_1714').lower_bound = solution2.objective_value*1.00001
        print ('Glucose uptake rate: ', solution2.objective_value)
        model_tmp.reactions.get_by_id (idx[5]).lower_bound = -1000  # protein pool

        # minimize the usage of protein pools
        model_tmp.objective = {model_tmp.reactions.prot_pool_exchange : 1}
        solution_f = model_tmp.optimize ()
        flux_max = solution_f.fluxes
        result = pd.DataFrame ({'rxnID' : flux_max.index, 'flux' : flux_max.values})
        result = result[result['rxnID'].str.contains ("draw_prot_")]
        result['geneID'] = result['rxnID'].str.replace ("draw_prot_", "")

        # compare with abundance
        abundance = pd.DataFrame ()
        column_name = f'abundance{i + 1}'
        if i <= 6 :
            column_names = [f'prot.{i + 1}']
            abundance[column_name] = proteomicsdata_Nlim[column_names]
            abundance1 = pd.concat ([abundance, proteomicsdata_Nlim['all_gene']], axis=1)
        else :
            column_names = [f'prot.{i - 6}']
            abundance[column_name] = proteomicsdata_Clim[column_names]
            abundance1 = pd.concat ([abundance, proteomicsdata_Clim['all_gene']], axis=1)
        result['pro_measured'] = singleMapping (abundance1[f"abundance{i + 1}"],
                                                abundance1["all_gene"],
                                                result['geneID'])
        result = result[~result["pro_measured"].isna ()]

        # pearsonr & pvalue
        from scipy.stats import pearsonr
        result1 = result[abs (result['flux']) > 0]
        result1 = result1[abs (result1['pro_measured']) > 0]
        number = len (result1)
        print (len (result1))
        # pearsonr & pvalue
        from scipy.stats import pearsonr, spearmanr
        result1 = result[abs (result['flux']) > 0]
        result1 = result1[abs (result1['pro_measured']) > 0]
        corr, pvalue = spearmanr (np.log10 (result1['pro_measured']), np.log10 (result1['flux']))
        print ("Correlation coefficient:", corr)
        print ("Correlation p_value:", pvalue)
        #corr_all.append (corr)

        # pro_evaluation
        result1['calculated'] = abs (np.log2 ((result1['flux']/result1['pro_measured'])))
        # error = np.sqrt (result1['calculated'].sum () / len (result1['calculated']))
        error = result1['calculated'].sum ()/len (result1['calculated'])

        toy_data_tmp.loc[i] = [error, corr, pvalue,number]
    #print(toy_data_tmp)
    #toy_data.loc[j] = [toy_data_tmp['error'].mean () , toy_data_tmp['corr'].mean () , toy_data_tmp['ss'].mean ()]
    return toy_data_tmp['error'].mean(), toy_data_tmp['corr'].mean(), toy_data_tmp['pvalue'].mean(), toy_data_tmp['number'].mean()

def getRxnByReactionName(model, name):
    '''
    This function is used to extract the rxn id based on rxn name
    It is suitable for ecGEMs as multiple reations could use the same name
    '''
    s = []
    for rxn in model.reactions:
        if name == rxn.name:
            #print(rxn.id)
            s.append(rxn.id)
    return s

def getExchangeRxns(model , reactionType="both"):
    '''
    Retrieves the exchange reactions from a model
    Exchange reactions are defined as reactions which involve only products
    or only reactants. If the unconstrained field is present, then that is
    used instead.

    model: a model structure
    reactionType: retrieve all reactions ('both'), only production ('out'),
    or only consumption ('in') (opt, default 'both')

    exchangeRxns: cell array with the IDs of the exchange reactions
    exchangeRxnsIndexes: vector with the indexes of the exchange reactions
    '''
    hasNoProducts = []
    hasNoReactants = []
    for i in range (len (model['rxns'])):
        isSub = model['S'].getcol (i)
        if reactionType=="both" or reactionType=="out":
            if not (isSub > 0).getnnz ():
                hasNoProducts.append (i)
        if reactionType=="both" or reactionType=="in":
            if not (isSub < 0).getnnz ():
                hasNoReactants.append (i)
    exchangeRxnIndexes = hasNoProducts + hasNoReactants
    exchangeRxns = model["rxns"][exchangeRxnIndexes]
    return exchangeRxnIndexes , exchangeRxns

def singleMapping (description, item1, item2, dataframe=True):
    """get the single description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    # used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    index = [None]*len(item2)
    result = [None]*len(item2)
    tt = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index[i] = item1.index(item2[i])
            result[i] = description[index[i]]
        else:
            index[i] = None
            result[i] = None
    return result

def model_execution (params,kcat_list, model, eModel, growthdata_N, rxn2block) :
    with eModel as ecmodel_tmp:
        for i, kcat_value in enumerate(params):
            rxnID = kcat_list.iloc[i]['rxnID']
            target_gene = kcat_list.iloc[i]['gene']
            ecmodel_tmp = updateEcGEMkcat(ecModel=ecmodel_tmp, target_gene=target_gene, rxnID=rxnID, kcat_m=kcat_value)
        result = abc_python_max(model, ecmodel_tmp, growthdata_N, 1, 1, rxn2block)
    return result[0][0][0], target_gene, rxnID

def parallel_turnover (j , growthdata , eModel , model , kcat_final , proteomicsdata_Clim, proteomicsdata_Nlim , miu) :
    if j % 100==0:
        print (f'current set processed: {j}')

    index = 'kcat_value' + str (j)
    toy_data_tmp = pd.DataFrame (columns=['error' , 'corr' , 'pvalue','number'])
    corr_all = []
    ex_mets = ['biomass pseudoreaction' , 'D-glucose exchange' , 'acetate exchange' , 'ethanol exchange' ,
               'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
               'oxygen exchange' , 'prot_pool_exchange']
    # find the related rxnID
    idx = []
    for name0 in ex_mets:
        # print(name0)
        s = getRxnByReactionName (model=eModel , name=name0)
        if len (s) > 1:
            print ("need check")
        elif len (s)==1:
            idx.append (s[0])
    # model_tmp = eModel.copy ()
    for i in range (len (growthdata)) :
        with eModel as model_tmp:
            try:
                dilutionrate = growthdata[i, 0]
                if miu!=dilutionrate :
                    continue
                model_tmp.solver = "cplex"
                for k, row in kcat_final.iterrows () :
                    gene = row['geneID']
                    mean = row[f'vsyn_value{j}']
                    kdeg_x = np.maximum (row['kdeg_x'], row['kdeg_y'])
                    kdeg_y = np.minimum (row['kdeg_x'], row['kdeg_y'])
                    reaction_id = f'draw_prot_{gene}'
                    if reaction_id in model_tmp.reactions :
                        # model_tmp.reactions.get_by_id (reaction_id).upper_bound = mean/(dilutionrate + kdeg)
                        model_tmp.reactions.get_by_id (reaction_id).bounds = mean/(dilutionrate + kdeg_x), mean/(
                                    dilutionrate + kdeg_y)
                model_tmp = changeMedia (model, model_tmp, 'D-glucose', "MIN")
                # model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
                model_tmp.reactions.get_by_id ('r_1631').bounds = 0, 0
                model_tmp.reactions.get_by_id ('r_1761').upper_bound = growthdata[i, 2]  # Ethanol exchange
                model_tmp.reactions.get_by_id ('r_1714').lower_bound = -growthdata[i, 1]  # D-glucose exchange
                # model_tmp.reactions.get_by_id ('r_1672').upper_bound = growthdata[i, 3]  # CO2 exchange
                # model_tmp.reactions.get_by_id ('r_1992').lower_bound = -growthdata[i, 4]  # oxygen exchange
                model_tmp.reactions.get_by_id ("prot_pool_exchange").upper_bound = growthdata[i, 5]*1000  # protein
                # model_tmp.reactions.get_by_id ("r_1672").lower_bound = -growthdata[i, 6]  # NH4 exchange
                print (model_tmp.optimize ().objective_value)
                model_tmp.reactions.get_by_id ('r_1714').lower_bound = -1000  # glucose uptake
                model_tmp.reactions.get_by_id (idx[0]).lower_bound = dilutionrate
                model_tmp.objective = { model_tmp.reactions.r_1714 : 1 }  # minimize the uptake of glucose
                solution2 = model_tmp.optimize ()
                print (solution2.objective_value)
                # then fix glucose uptake and minimize the protein pool
                model_tmp.reactions.get_by_id ('r_1714').lower_bound = solution2.objective_value*1.00001
                print ('Glucose uptake rate: ', solution2.objective_value)
                model_tmp.reactions.get_by_id (idx[5]).lower_bound = -1000  # protein pool

                # minimize the usage of protein pools
                model_tmp.objective = { model_tmp.reactions.prot_pool_exchange : 1 }
                solution_f = model_tmp.optimize ()
                flux_max = solution_f.fluxes
                result = pd.DataFrame ({ 'rxnID' : flux_max.index, 'flux' : flux_max.values })
                result = result[result['rxnID'].str.contains ("draw_prot_")]
                result['geneID'] = result['rxnID'].str.replace ("draw_prot_", "")

                # compare with abundance
                abundance = pd.DataFrame ()
                column_name = f'abundance{i + 1}'
                column_names = [f'prot.{i + 1}']
                abundance[column_name] = proteomicsdata_Nlim[column_names]
                abundance1 = pd.concat ([abundance, proteomicsdata_Nlim['all_gene']], axis=1)

                result['pro_measured'] = singleMapping (abundance1[f"abundance{i + 1}"],
                                                        abundance1["all_gene"],
                                                        result['geneID'])
                result = result[~result["pro_measured"].isna ()]

                # pearsonr & pvalue
                from scipy.stats import pearsonr
                result1 = result[abs (result['flux']) > 0]
                result1 = result1[abs (result1['pro_measured']) > 0]
                number = len (result1)
                print (len (result1))
                # pearsonr & pvalue
                from scipy.stats import pearsonr, spearmanr
                result1 = result[abs (result['flux']) > 0]
                result1 = result1[abs (result1['pro_measured']) > 0]
                corr, pvalue = spearmanr (np.log10 (result1['pro_measured']), np.log10 (result1['flux']))
                print ("Correlation coefficient:", corr)
                print ("Correlation p_value:", pvalue)
                # corr_all.append (corr)

                # pro_evaluation
                result1['calculated'] = abs (np.log2 ((result1['flux']/result1['pro_measured'])))
                # error = np.sqrt (result1['calculated'].sum () / len (result1['calculated']))
                error = result1['calculated'].sum ()/len (result1['calculated'])

                toy_data_tmp.loc[i] = [error, corr, pvalue, number]
            except Exception as e:
                print (f"Error at index {j}, iteration {i}: {e}")
                continue
    #print(toy_data_tmp)
    #toy_data.loc[j] = [toy_data_tmp['error'].mean () , toy_data_tmp['corr'].mean () , toy_data_tmp['ss'].mean ()]
    return toy_data_tmp['error'].mean(), toy_data_tmp['corr'].mean(), toy_data_tmp['pvalue'].mean(), toy_data_tmp['number'].mean()

def rmsecal_new(model , model_cobra , data , constrain , objective , osenseStr , rxn2block):
    data = np.array (data)
    rmse_tmp = []
    simulated = np.zeros ((len (data[: , 0]) , 9))
    rxnNames = model['rxnNames']
    for i in range (len (data[: , 0])):
        exp = np.array (data[: , 2:11])  # u sub ace eth gly pyr ethyl_acetate co2 o2
        exp = exp * [1 , -1 , 1 , 1 , 1 , 1 , 1 , 1 , -1]
        exp = exp.astype (float)
        ex_mets = ['growth' , data[i , 1][0] + " exchange" , 'acetate exchange' , 'ethanol exchange' ,
                   'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
                   'oxygen exchange']
        idx = []
        temp = []
        for k in range (len (ex_mets)):
            temp = np.where (rxnNames==ex_mets[k])
            idx.append (temp)
        idx = np.array (idx)
        idx = np.transpose (idx[: , 0])
        # Create a temp model
        with model_cobra as model_tmp:
            # Suitable for different carbon sources
           # model_tmp = changeMedia (model , model_tmp , data[i , 1][0] , data[i , 15])  # 'D-glucose'
            model_tmp = changeMedia (model , model_tmp , 'D-glucose' , data[i , 15])
            model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
            model_tmp.reactions.get_by_id ('r_1631').bounds = 0 , 0

            # Suitable for anaerobic
            if data[i , 13]=="anaerobic" or data[i , 13]=="limited":
                model_tmp = anaerobicModel (model , model_tmp)

            if data[i , 13]=="limited":
                model_tmp.reactions.get_by_id('r_1992').lower_bound = -5

            model_tmp.reactions.get_by_id ('r_1714').lower_bound = 0

            id_substrate = np.where (model['rxnNames']==(data[i , 1][0] + " exchange"))
            if constrain==False:
                model_tmp.reactions[idx[0][1]].lower_bound = -1000
            else:
                model_tmp.reactions[idx[0][1]].lower_bound = exp[i , 1]

            # Solve the temp model
            sol_tmp = model_tmp.optimize ()
            print(sol_tmp.objective_value)
            sol = sol_tmp.fluxes  # sol[:,i] = sol_tmp.fluxes

        print ("No. " + str (i + 1) + "finish !")

        tmp = np.where (~np.isnan (exp[i]))[0]

        # Normalize the number of carbon
        excarbon = model['excarbon'][: , idx[0]]
        for x in range (len (idx[0 , :])):
            if excarbon[: , x]==0:
                excarbon[: , x] = 1
        exp_tmp = []
        for s in range (len (tmp)):
            exp_tmp.append (exp[i , tmp[s]] * excarbon[: , tmp[s]])
        sol_idx = np.transpose (sol[idx[0]])
        sol_idx = np.array (sol_idx)
        simulated_tmp = []
        for t in range (len (tmp)):
            simulated_tmp.append (
                np.array (sol_idx)[tmp[t]] * excarbon[: , tmp[t]])  # normalize the growth rate issue by factor 10

        # The expected blocked reactions
        exp_block = np.zeros ((1 , len(rxn2block)))  # 218怎么替换
        rxnblockidx_pre = np.setdiff1d (rxn2block , model['rxns'][idx[0][1]])
        rxnblockidx = []
        for k in range (len (model['rxns'])):
            if model['rxns'][k] in rxnblockidx_pre:
                rxnblockidx.append (k)

        simulated_block = []
        for t in range (len (rxnblockidx)):
            simulated_block.append (np.array (sol)[rxnblockidx[t]] * model['excarbon'][: , rxnblockidx[t]])
        id_zero = np.where (np.array (simulated_block)!=0)
        print (f"id_zero: {id_zero}")

        exp_block = exp_block[: , id_zero[0]]
        print (f"exp_block shape: {exp_block.shape}")
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

def abc_python_max_new(model , model_cobra , growthdata , max_growth , rxn2block ):
    # get carbonnum for each exchange rxn to further calculation of error

    objective = 'r_2111'
    osenseStr = 'max'
    start_time = time.time ()

    if len (growthdata)!=0:
        rmse_1 , exp_1 , simulated_1 = rmsecal_new (model , model_cobra , growthdata , True , objective , osenseStr ,
                                                    rxn2block)
    else:
        rmse_1 = []
        exp_1 = []
        simulated_1 = []
    # second searcch for maxmial growth rate without constrain
    if len (max_growth)!=0:
        rmse_2 , exp_2 , simulated_2 = rmsecal_new (model , model_cobra , max_growth , False , objective , osenseStr ,
                                                    rxn2block)
    else:
        rmse_2 = []
        exp_2 = []
        simulated_2 = []
    end_time = time.time ()
    execution_time = end_time - start_time
    print ("execution time: " , execution_time , " seconds")

    exp = np.array ([exp_1 , exp_2],dtype=object)
    simulated = np.array ([simulated_1 , simulated_2],dtype=object)
    rmse = np.array ([rmse_1 , rmse_2],dtype=object)
    rmse_final = np.nanmean (rmse , 0)

    print ("rmse_final is " , rmse_final)

    return rmse_final , exp , simulated

def with_pro_data(model_tmp, growthdata, model,idx, proteomicsdata_Nlim,proteomicsdata_Clim, toy_data_tmp, i,j):
    try:
        with model_tmp as model_pro_tmp:
            dilutionrate = growthdata[ i, 0 ]
            model_pro_tmp.solver = "cplex"

            model_pro_tmp = changeMedia (model, model_pro_tmp, 'D-glucose', "MIN")
            # model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
            model_pro_tmp.reactions.get_by_id ('r_1631').bounds = 0, 0
            model_pro_tmp.reactions.get_by_id ('r_1761').upper_bound = growthdata[ i, 2 ]  # Ethanol exchange
            model_pro_tmp.reactions.get_by_id ('r_1714').lower_bound = -growthdata[ i, 1 ]  # D-glucose exchange
            model_pro_tmp.reactions.get_by_id ('r_1672').upper_bound = growthdata[ i, 3 ]  # CO2 exchange
            model_pro_tmp.reactions.get_by_id ('r_1992').lower_bound = -growthdata[ i, 4 ]  # oxygen exchange
            model_pro_tmp.reactions.get_by_id ("prot_pool_exchange").upper_bound = growthdata[
                                                                                       i, 5 ] * 1000  # protein
            if i <= 6:
                model_pro_tmp.reactions.get_by_id ("r_1672").lower_bound = -growthdata[ i, 6 ]  # NH4 exchange
            print (model_pro_tmp.optimize ().objective_value)
            model_pro_tmp.reactions.get_by_id ('r_1714').lower_bound = -1000  # glucose uptake
            model_pro_tmp.reactions.get_by_id (idx[ 0 ]).lower_bound = dilutionrate
            model_pro_tmp.objective = { model_pro_tmp.reactions.r_1714: 1 }  # minimize the uptake of glucose
            solution2 = model_pro_tmp.optimize ()
            print (solution2.objective_value)
            # then fix glucose uptake and minimize the protein pool
            model_pro_tmp.reactions.get_by_id ('r_1714').lower_bound = solution2.objective_value * 1.00001
            print ('Glucose uptake rate: ', solution2.objective_value)
            model_pro_tmp.reactions.get_by_id (idx[ 5 ]).lower_bound = -1000  # protein pool

            # minimize the usage of protein pools
            model_pro_tmp.objective = { model_pro_tmp.reactions.prot_pool_exchange: 1 }
            solution_f = model_pro_tmp.optimize ()
            flux_max = solution_f.fluxes
            result = pd.DataFrame ({ 'rxnID': flux_max.index, 'flux': flux_max.values })
            result = result[ result[ 'rxnID' ].str.contains ("draw_prot_") ]
            result[ 'geneID' ] = result[ 'rxnID' ].str.replace ("draw_prot_", "")

            # compare with abundance
            abundance = pd.DataFrame ()
            column_name = f'abundance{i + 1}'
            if i <= 6:
                column_names = [ f'prot.{i + 1}' ]
                abundance[ column_name ] = proteomicsdata_Nlim[ column_names ]
                abundance1 = pd.concat ([ abundance, proteomicsdata_Nlim[ 'all_gene' ] ], axis=1)
            else:
                column_names = [ f'prot.{i - 6}' ]
                abundance[ column_name ] = proteomicsdata_Clim[ column_names ]
                abundance1 = pd.concat ([ abundance, proteomicsdata_Clim[ 'all_gene' ] ], axis=1)
            result[ 'pro_measured' ] = singleMapping (abundance1[ f"abundance{i + 1}" ],
                                                      abundance1[ "all_gene" ],
                                                      result[ 'geneID' ])
            result = result[ ~result[ "pro_measured" ].isna () ]

            # pearsonr & pvalue
            from scipy.stats import pearsonr
            result1 = result[ abs (result[ 'flux' ]) > 0 ]
            result1 = result1[ abs (result1[ 'pro_measured' ]) > 0 ]
            number = len (result1)
            print (len (result1))
            # pearsonr & pvalue
            from scipy.stats import pearsonr, spearmanr
            result1 = result[ abs (result[ 'flux' ]) > 0 ]
            result1 = result1[ abs (result1[ 'pro_measured' ]) > 0 ]
            corr, pvalue = spearmanr (np.log10 (result1[ 'pro_measured' ]), np.log10 (result1[ 'flux' ]))
            print ("Correlation coefficient:", corr)
            print ("Correlation p_value:", pvalue)
            # corr_all.append (corr)

            # pro_evaluation
            result1[ 'calculated' ] = abs (np.log2 ((result1[ 'flux' ] / result1[ 'pro_measured' ])))
            # error = np.sqrt (result1['calculated'].sum () / len (result1['calculated']))
            error = result1[ 'calculated' ].sum () / len (result1[ 'calculated' ])

            flux_ratio = result1[ 'flux' ] / result1[ 'pro_measured' ]
            error_condition = (flux_ratio > 50) | (flux_ratio < 0.02)
            error_num = error_condition.sum()  # 满足条件的数量

            toy_data_tmp.at[ i,'error' ] = error
            toy_data_tmp.at[ i,'corr' ] = corr
            toy_data_tmp.at[ i,'pvalue' ] = pvalue
            toy_data_tmp.at[ i,'number' ] = number
            toy_data_tmp.at[ i,'error_num' ] = error_num
            error_final = toy_data_tmp[ 'error' ].mean()
            corr_final = toy_data_tmp[ 'corr' ].mean()
            pvalue_final = toy_data_tmp[ 'pvalue' ].mean()
            number_final = toy_data_tmp[ 'number' ].mean()
            error_num_final = toy_data_tmp[ 'error_num' ].mean()


    except Exception as e:
        print (f"Error at index {j}, iteration {i}: {e}")
        error_final = np.nan
        corr_final = np.nan
        pvalue_final = np.nan
        number_final = np.nan
        error_num_final = np.nan
    return error_final, corr_final, pvalue_final, number_final,error_num_final

def ecYeastMinimalMedia(model):
    """
    This function is used to define a simple media for ecYeast
    :param model:
    :return: a model with the defined the minimal media
    """
    rxnID = []
    rxnName = []
    for i, x in enumerate(model.reactions):
        rxnID.append(x.id)
        rxnName.append(x.name)

    exchange_rxn =[x for x, y in zip(rxnID, rxnName) if '_rvs' in x and 'exchange' in y]
    # first block any uptake
    for i, x in enumerate(exchange_rxn):
        rxn0 = exchange_rxn[i]
        #print(rxn0)
        model.reactions.get_by_id(rxn0).upper_bound = 0

    #Allow uptake of essential components
    model.reactions.get_by_id("r_1654_REV").upper_bound = 1000 #ammonium exchange (reversible)
    model.reactions.get_by_id("r_1861_REV").upper_bound = 1000 #iron(2+) exchange (reversible)
    model.reactions.get_by_id("r_2100_REV").upper_bound = 1000 #water exchange (reversible)
    model.reactions.get_by_id("r_1992_REV").upper_bound = 1000 #oxygen exchange (reversible)
    model.reactions.get_by_id("r_2005_REV").upper_bound = 1000 #phosphate exchange (reversible)
    model.reactions.get_by_id("r_2060").upper_bound = 1000 #sulphate exchange (reversible)
    model.reactions.get_by_id("r_1832_REV").upper_bound = 1000 #H+ exchange (reversible)
    return model

def ecModelSimulate(model, idx, growthdata,i):
    # model_tmp = model.copy()
    with model as model_tmp:
        # model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
        model_tmp.reactions.get_by_id ('r_1631').bounds = 0 ,0
        model_tmp.reactions.get_by_id ('r_1714').lower_bound = -growthdata[i ,2]  # D-glucose exchange
        model_tmp.reactions.get_by_id ("prot_pool_exchange").upper_bound = growthdata[i ,3] * 1000  # protein
        print ('Growth rate: ' ,model_tmp.optimize ().objective_value)
        model_tmp.reactions.get_by_id ('r_1714').lower_bound = -1000  # glucose uptake
        model_tmp.reactions.get_by_id (idx[0]).lower_bound = growthdata[i ,1]
        model_tmp.objective = { model_tmp.reactions.r_1714: 1 }  # minimize the uptake of glucose
        solution2 = model_tmp.optimize ()
        print (solution2.objective_value)
        # then fix glucose uptake and minimize the protein pool
        model_tmp.reactions.get_by_id ('r_1714').lower_bound = solution2.objective_value * 1.0001
        print ('Glucose uptake rate: ' ,solution2.objective_value)
        model_tmp.reactions.get_by_id (idx[5]).lower_bound = -1000  # protein pool

        # minimize the usage of protein pools
        model_tmp.objective = { model_tmp.reactions.prot_pool_exchange: 1 }
        solution_f = model_tmp.optimize ()
    return solution_f

def updateEcGEMkcat_new(ecModel, target_gene, rxnID, kcat_m):
    coef = 1 / kcat_m
    original_reaction = ecModel.reactions.get_by_id(rxnID).reaction

    if '-->' in original_reaction:
        arrow = '-->'
    elif '<->' in original_reaction:
        arrow = '<->'
    else:
        raise ValueError("No reaction arrow found in the original reaction.")

    substrates, products = original_reaction.split(arrow)
    substrates = substrates.strip()
    products = products.strip()

    substrate_parts = substrates.split(" + ")
    product_parts = products.split(" + ")

    updated_substrates = []
    updated_products = []

    for part in substrate_parts:
        parts = part.split()
        quantity = parts[0]
        species = " ".join(parts[1:])
        if target_gene in species:
            new_part = f"{coef} {species}"
            updated_substrates.append(new_part)
        else:
            updated_substrates.append(part)

    for part in product_parts:
        parts = part.split()
        quantity = parts[0]
        species = " ".join(parts[1:])
        if target_gene in species:
            new_part = f"{coef} {species}"
            updated_products.append(new_part)
        else:
            updated_products.append(part)

    updated_reaction = " + ".join(updated_substrates) + " " + arrow + " " + " + ".join(updated_products)
    ecModel.reactions.get_by_id(rxnID).reaction = updated_reaction

    # print(f"Updated reaction: {updated_reaction}")  # test
    return ecModel


def parallel_task_preprocessing_all (j , growthdata ,growthdata_train, max_growth_train, eModel , model , kcat_final , proteomicsdata_Clim, proteomicsdata_Nlim,  train_type, rxn2block) :
    if j % 100==0:
        print (f'current set processed: {j}')

    index = 'kcat_value' + str (j)
    if 'phenotype' in train_type :
        if "protein" in train_type:
            toy_data_tmp = pd.DataFrame (columns=['error' , 'corr' , 'pvalue','number', 'error_num', 'RMSE'])
        else:
            toy_data_tmp = pd.DataFrame (columns=['RMSE'])
    else:
        toy_data_tmp = pd.DataFrame (columns=['error' , 'corr' , 'pvalue','number','error_num'])
    corr_all = []
    ex_mets = ['biomass pseudoreaction' , 'D-glucose exchange' , 'acetate exchange' , 'ethanol exchange' ,
               'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
               'oxygen exchange' , 'prot_pool_exchange']
    # find the related rxnID
    idx = []
    for name0 in ex_mets:
        # print(name0)
        s = getRxnByReactionName (model=eModel , name=name0)
        if len (s) > 1:
            print ("need check")
        elif len (s)==1:
            idx.append (s[0])
    # model_tmp = eModel.copy ()
    with eModel as model_tmp:
        for k in range (len (kcat_final)):
            target_gene0 = kcat_final[ 'gene' ][ k ]
            kcat_m0 = kcat_final[ index ][ k ]
            rxn0 = kcat_final[ 'rxnID' ][ k ]
            # model_tmp = updateEcGEMkcat_new (ecModel=model_tmp, target_gene=target_gene0, rxnID=rxn0, kcat_m=kcat_m0)
            model_tmp = updateEcGEMkcat (ecModel=model_tmp ,target_gene=target_gene0 ,rxnID=rxn0 ,kcat_m=kcat_m0)
        # model_tmp.objective = {model_tmp.reactions.r_2111 : -1}
        print (model_tmp.optimize ().objective_value)
        error_final = 0
        corr_final = 0
        pvalue_final = 0
        number_final = 0
        RMSE_final = 0
        error_num_final = 0
        # train with proteomics
        # model_pro_tmp = model_tmp.copy()
        if "protein" in train_type:
            print ("yes_protein")
            for i in range (len (growthdata)):
                # result = with_pro_data_iIsor (model_tmp, growthdata, model, idx, proteomicsdata_Clim, toy_data_tmp, i, j)
                result = with_pro_data(model_tmp, growthdata, model,idx, proteomicsdata_Nlim,proteomicsdata_Clim, toy_data_tmp, i,j)
                print (result)
                error_final = result[ 0 ]
                corr_final = result[ 1 ]
                pvalue_final = result[ 2 ]
                number_final = result[ 3 ]
                error_num_final = result[ 4 ]

        # train with growthdata
        # model_phen_tmp = model_tmp.copy()
        if "phenotype" in train_type:
            print ("yes_phen")
            with model_tmp as model_phen_tmp:
                model_phen_tmp.reactions.get_by_id ("prot_pool_exchange").upper_bound = 230  # protein 178.5ilsor
                RMSE = abc_python_max_new (model, model_phen_tmp, growthdata_train, max_growth_train, rxn2block)
                RMSE_final = RMSE[ 0 ]


    return error_final, corr_final, pvalue_final, number_final, error_num_final, RMSE_final

def with_pro_data_other_yeast(model_tmp, model,idx, proteomicsdata_Nlim, toy_data_tmp, i,j,strain,prot):
    try:
        with model_tmp as model_pro_tmp:
        # model_pro_tmp = model_tmp.copy()
            if strain=="kla" :
                dilutionrate = 0.438
            else :
                dilutionrate = 0.1
            model_pro_tmp.solver = "cplex"
            model_pro_tmp = changeMedia( model,model_pro_tmp,'D-glucose',"MIN" )
            # model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
            model_pro_tmp.reactions.get_by_id( 'r_1631' ).bounds = 0,0
            model_pro_tmp.reactions.get_by_id( "prot_pool_exchange" ).upper_bound = prot  # protein
            print( model_pro_tmp.optimize().objective_value )
            model_pro_tmp.reactions.get_by_id( 'r_1714' ).lower_bound = -1000  # glucose uptake
            model_pro_tmp.reactions.get_by_id( idx[ 0 ] ).lower_bound = dilutionrate
            model_pro_tmp.objective = { model_pro_tmp.reactions.r_1714 : 1 }  # minimize the uptake of glucose
            solution2 = model_pro_tmp.optimize()
            print( solution2.objective_value )
            # then fix glucose uptake and minimize the protein pool
            model_pro_tmp.reactions.get_by_id( 'r_1714' ).lower_bound = solution2.objective_value * 1.00001
            print( 'Glucose uptake rate: ',solution2.objective_value )
            model_pro_tmp.reactions.get_by_id( idx[ 9 ] ).upper_bound = 1000  # protein pool

            # minimize the usage of protein pools
            model_pro_tmp.objective = { model_pro_tmp.reactions.prot_pool_exchange : 1 }
            solution_f = model_pro_tmp.optimize()
            flux_max = solution_f.fluxes
            result = pd.DataFrame( { 'rxnID' : flux_max.index,'flux' : flux_max.values } )
            result = result[ result[ 'rxnID' ].str.contains( "draw_prot_" ) ]
            result[ 'geneID' ] = result[ 'rxnID' ].str.replace( "draw_prot_","" )

            # compare with abundance
            abundance = pd.DataFrame()
            column_name = f'abundance1'
            column_names = [ f'{strain}_prot.1' ]
            abundance[ column_name ] = proteomicsdata_Nlim[ column_names ]
            abundance1 = pd.concat( [ abundance,proteomicsdata_Nlim[ 'geneID' ] ],axis=1 )
            result[ 'pro_measured' ] = singleMapping( abundance1[ f"abundance{i + 1}" ],
                                                      abundance1[ "geneID" ],
                                                      result[ 'geneID' ] )
            result = result[ ~result[ "pro_measured" ].isna() ]

            # pearsonr & pvalue
            from scipy.stats import pearsonr
            result1 = result[ abs( result[ 'flux' ] ) > 0 ]
            result1 = result1[ abs( result1[ 'pro_measured' ] ) > 0 ]
            number = len( result1 )
            print( len( result1 ) )
            # pearsonr & pvalue
            from scipy.stats import pearsonr,spearmanr
            result1 = result[ abs( result[ 'flux' ] ) > 0 ]
            result1 = result1[ abs( result1[ 'pro_measured' ] ) > 0 ]
            corr,pvalue = spearmanr( np.log10( result1[ 'pro_measured' ] ),np.log10( result1[ 'flux' ] ) )
            print( "Correlation coefficient:",corr )
            print( "Correlation p_value:",pvalue )
            # corr_all.append (corr)

            # pro_evaluation
            result1[ 'calculated' ] = abs( np.log2( (result1[ 'flux' ] / result1[ 'pro_measured' ]) ) )
            # error = np.sqrt (result1['calculated'].sum () / len (result1['calculated']))
            error = result1[ 'calculated' ].sum() / len( result1[ 'calculated' ] )

            flux_ratio = result1[ 'flux' ] / result1[ 'pro_measured' ]
            error_condition = (flux_ratio > 50) | (flux_ratio < 0.02)
            error_num = error_condition.sum()  # 满足条件的数量

            toy_data_tmp.at[ i,'error' ] = error
            toy_data_tmp.at[ i,'corr' ] = corr
            toy_data_tmp.at[ i,'pvalue' ] = pvalue
            toy_data_tmp.at[ i,'number' ] = number
            toy_data_tmp.at[ i,'error_num' ] = error_num
            error_final = toy_data_tmp[ 'error' ].mean()
            corr_final = toy_data_tmp[ 'corr' ].mean()
            pvalue_final = toy_data_tmp[ 'pvalue' ].mean()
            number_final = toy_data_tmp[ 'number' ].mean()
            error_num_final = toy_data_tmp[ 'error_num' ].mean()

    except Exception as e:
        print (f"Error at index {j}, iteration {i}: {e}")
        error_final = np.nan
        corr_final = np.nan
        pvalue_final = np.nan
        number_final = np.nan
        error_num_final = np.nan
    return error_final, corr_final, pvalue_final, number_final, error_num_final

def parallel_task_other_yeast (j ,growthdata_train, max_growth_train, eModel , model , kcat_final , proteomicsdata_Nlim,  train_type, rxn2block,strain,prot) :
    if j % 100==0:
        print (f'current set processed: {j}')

    index = 'kcat_value' + str (j)
    if 'phenotype' in train_type :
        if "protein" in train_type:
            toy_data_tmp = pd.DataFrame (columns=['error' , 'corr' , 'pvalue','number','error_num','RMSE'])
        else:
            toy_data_tmp = pd.DataFrame (columns=['RMSE'])
    else:
        toy_data_tmp = pd.DataFrame (columns=['error' , 'corr' , 'pvalue','number','error_num'])
    ex_mets = ['biomass pseudoreaction' , 'D-glucose exchange' , 'acetate exchange' , 'ethanol exchange' ,
               'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
               'oxygen exchange' , 'prot_pool_exchange']
    # find the related rxnID
    idx = []
    for name0 in ex_mets:
        # print(name0)
        s = getRxnByReactionName (model=eModel , name=name0)
        if len (s) > 1:
            print ("need check")
        elif len (s)==1:
            idx.append (s[0])
    # model_tmp = eModel.copy ()
    with eModel as model_tmp:
        for k in range( len( kcat_final ) ) :
            target_gene0 = kcat_final[ 'gene' ][ k ]
            kcat_m0 = kcat_final[ index ][ k ]
            rxn0 = kcat_final[ 'rxnID' ][ k ]
            model_tmp = updateEcGEMkcat_new( ecModel=model_tmp ,target_gene=target_gene0 ,rxnID=rxn0 ,kcat_m=kcat_m0 )
            # model_tmp = updateEcGEMkcat (ecModel=model_tmp ,target_gene=target_gene0 ,rxnID=rxn0 ,kcat_m=kcat_m0)
            # model_tmp.objective = {model_tmp.reactions.r_2111 : -1}
        print( model_tmp.optimize().objective_value )
        error_final = 0
        corr_final = 0
        pvalue_final = 0
        number_final = 0
        error_num_final = 0
        RMSE_final = 0
        # train with proteomics
        # model_pro_tmp = model_tmp.copy()
        if "protein" in train_type :
            print( "yes_protein" )
            for i in range( 0 ,1 ) :
                result = with_pro_data_other_yeast( model_tmp ,model ,idx ,proteomicsdata_Nlim ,toy_data_tmp ,i ,j ,strain ,prot )
                # result = with_pro_data(model_tmp, growthdata, model,idx, proteomicsdata_Nlim,proteomicsdata_Clim, toy_data_tmp, i,j)
                print( result )
                error_final = result[ 0 ]
                corr_final = result[ 1 ]
                pvalue_final = result[ 2 ]
                number_final = result[ 3 ]
                error_num_final = result[ 4 ]

        # train with growthdata
        # model_phen_tmp = model_tmp.copy()
        if "phenotype" in train_type :
            print( "yes_phen" )
            with model_tmp as model_phen_tmp :
                model_phen_tmp.reactions.get_by_id( "prot_pool_exchange" ).upper_bound = prot  # protein 178.5ilsor
                RMSE = abc_python_max_new( model ,model_phen_tmp ,growthdata_train ,max_growth_train ,rxn2block )
                RMSE_final = RMSE[ 0 ]


    return error_final, corr_final, pvalue_final, number_final, error_num_final, RMSE_final

def with_pro_data_CN(model_pro_tmp, growthdata, model,idx,proteomicsdata_Clim, toy_data_tmp, k,j):

    dilutionrate = 0.2
    model_pro_tmp.solver = "cplex"

    model_pro_tmp = changeMedia (model, model_pro_tmp, 'D-glucose', "MIN")
            # model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
    model_pro_tmp.reactions.get_by_id ('r_1631').bounds = 0, 0
    model_pro_tmp.reactions.get_by_id ("r_1672").lower_bound = -growthdata[ k, 3 ]  # NH4 exchange
    print (model_pro_tmp.optimize ().objective_value)
    model_pro_tmp.reactions.get_by_id ('r_1714').lower_bound = -1000  # glucose uptake
    model_pro_tmp.reactions.get_by_id (idx[ 0 ]).lower_bound = dilutionrate
    model_pro_tmp.objective = { model_pro_tmp.reactions.r_1714: 1 }  # minimize the uptake of glucose
    solution2 = model_pro_tmp.optimize ()
    print (solution2.objective_value)
            # then fix glucose uptake and minimize the protein pool
    model_pro_tmp.reactions.get_by_id ('r_1714').lower_bound = solution2.objective_value * 1.00001
    print ('Glucose uptake rate: ', solution2.objective_value)
    model_pro_tmp.reactions.get_by_id ('prot_pool_exchange').lower_bound = -1000  # protein pool

            # minimize the usage of protein pools
    model_pro_tmp.objective = { model_pro_tmp.reactions.prot_pool_exchange: 1 }
    solution_f = model_pro_tmp.optimize ()
    flux_max = solution_f.fluxes
    result = pd.DataFrame ({ 'rxnID': flux_max.index, 'flux': flux_max.values })
    result = result[ result[ 'rxnID' ].str.contains ("draw_prot_") ]
    result[ 'geneID' ] = result[ 'rxnID' ].str.replace ("draw_prot_", "")

            # compare with abundance
    abundance = pd.DataFrame ()
    column_name = f'abundance{k}'
    column_names = [ f'C/N_{k+1}' ]
    abundance[ column_name ] = proteomicsdata_Clim[ column_names ]
    abundance1 = pd.concat ([ abundance, proteomicsdata_Clim[ 'gene' ] ], axis=1)
    result[ 'pro_measured' ] = singleMapping (abundance1[ f"abundance{k}" ],
                                                      abundance1[ "gene" ],
                                                      result[ 'geneID' ])
    result = result[ ~result[ "pro_measured" ].isna () ]

            # pearsonr & pvalue
    from scipy.stats import pearsonr
    result1 = result[ abs (result[ 'flux' ]) > 0 ]
    result1 = result1[ abs (result1[ 'pro_measured' ]) > 0 ]
    number = len (result1)
    print (len (result1))
            # pearsonr & pvalue
    from scipy.stats import pearsonr, spearmanr
    result1 = result[ abs (result[ 'flux' ]) > 0 ]
    result1 = result1[ abs (result1[ 'pro_measured' ]) > 0 ]
    corr, pvalue = spearmanr (np.log10 (result1[ 'pro_measured' ]), np.log10 (result1[ 'flux' ]))
    print ("Correlation coefficient:", corr)
    print ("Correlation p_value:", pvalue)
            # corr_all.append (corr)

            # pro_evaluation
    result1[ 'calculated' ] = abs (np.log2 ((result1[ 'flux' ] / result1[ 'pro_measured' ])))
            # error = np.sqrt (result1['calculated'].sum () / len (result1['calculated']))
    error = result1[ 'calculated' ].sum () / len (result1[ 'calculated' ])

    flux_ratio = result1[ 'flux' ] / result1[ 'pro_measured' ]
    error_condition = (flux_ratio > 50) | (flux_ratio < 0.02)
    error_num = error_condition.sum()  # 满足条件的数量

    toy_data_tmp.at[ k,'error' ] = error
    toy_data_tmp.at[ k,'corr' ] = corr
    toy_data_tmp.at[ k,'pvalue' ] = pvalue
    toy_data_tmp.at[ k,'number' ] = number
    toy_data_tmp.at[ k,'error_num' ] = error_num
    error_final = toy_data_tmp[ 'error' ].mean()
    corr_final = toy_data_tmp[ 'corr' ].mean()
    pvalue_final = toy_data_tmp[ 'pvalue' ].mean()
    number_final = toy_data_tmp[ 'number' ].mean()
    error_num_final = toy_data_tmp[ 'error_num' ].mean()
    return error_final, corr_final, pvalue_final, number_final,error_num_final

def parallel_task_preprocessing_CN (j,k , growthdata, eModel , model , kcat_final , proteomicsdata_Clim, rxn2block) :
    if j % 100==0:
        print (f'current set processed: {j}')

    index = 'kcat_value' + str (j)
    toy_data_tmp = pd.DataFrame (columns=['error' , 'corr' , 'pvalue','number','error_num'])
    corr_all = []
    ex_mets = ['biomass pseudoreaction' , 'D-glucose exchange' , 'acetate exchange' , 'ethanol exchange' ,
               'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
               'oxygen exchange' , 'prot_pool_exchange']
    # find the related rxnID
    idx = []
    for name0 in ex_mets:
        # print(name0)
        s = getRxnByReactionName (model=eModel , name=name0)
        if len (s) > 1:
            print ("need check")
        elif len (s)==1:
            idx.append (s[0])
    # model_tmp = eModel.copy ()
    with eModel as model_tmp:
        for m in range (len (kcat_final)):
            target_gene0 = kcat_final[ 'gene' ][ m ]
            kcat_m0 = kcat_final[ index ][ m ]
            rxn0 = kcat_final[ 'rxnID' ][ m ]
            # model_tmp = updateEcGEMkcat_new (ecModel=model_tmp, target_gene=target_gene0, rxnID=rxn0, kcat_m=kcat_m0)
            model_tmp = updateEcGEMkcat (ecModel=model_tmp ,target_gene=target_gene0 ,rxnID=rxn0 ,kcat_m=kcat_m0)
        # model_tmp.objective = {model_tmp.reactions.r_2111 : -1}
        print (model_tmp.optimize ().objective_value)
        error_final = 0
        corr_final = 0
        pvalue_final = 0
        number_final = 0
        error_num_final = 0
        # train with proteomics
        # model_pro_tmp = model_tmp.copy()
        print("start working")
        # result = with_pro_data_iIsor (model_tmp, growthdata, model, idx, proteomicsdata_Clim, toy_data_tmp, i, j)
        result = with_pro_data_CN(model_tmp,growthdata,model,idx,proteomicsdata_Clim,toy_data_tmp,k,j)
        print("finish working")
        print(result)
        error_final = result[ 0 ]
        corr_final = result[ 1 ]
        pvalue_final = result[ 2 ]
        number_final = result[ 3 ]
        error_num_final = result[ 4 ]


    return error_final, corr_final, pvalue_final, number_final, error_num_final

