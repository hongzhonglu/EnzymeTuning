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

def abc_python_max(model , model_cobra , tmp ,proc , sample_generation ,rxn2block):
    nstep = sample_generation / proc
    nstep = int (nstep)
    rmse_final = np.zeros ((1 , nstep))

    # get carbonnum for each exchange rxn to further calculation of error
    if len (model["excarbon"])==0:  ##检查model中是否有excarbon, ~is_field
        model = addCarbonNum (model)  ##补写function
    for k in range (nstep):
        print ('nstep:' + str (k + 1) + '/' + str (nstep))

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
        print ("execution time: " , execution_time , " seconds")

        # print("RMSE_2 finish !")

        exp = exp_1
        simulated = simulated_1
        rmse = rmse_1
        rmse_final[0 , k] = np.nanmean (rmse , 0)

        # only output simulated result for one generation
        if nstep!=1 or sample_generation!=1:
            simulated = []
            exp = []
        print ("rmse_final is " , rmse_final)

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

from multiprocessing import Pool,cpu_count
#def parallel_task (index,row,model,ecModel,tmp_data,rxn2block,kcat_list):
#    print(index)
#    with ecModel as ecmodel_tmp:
#        target_gene0 = row['gene']
#        kcat_m0 = row['kcat']*foldchange
#        rxn0 = row['rxnID']
#        #ecmodel_tmp = updateEcGEMkcat(ecGEM=ecmodel_tmp, target_gene=target_gene0, rxnID=rxn0, kcat_m=kcat_m0)
#        result0 = abc_python_max(model , ecmodel_tmp , tmp_data ,1 , 1 , rxn2block)
#    return result0
def parallel_task(row, foldchange, model, ecModel, tmp_data, rxn2block):
    # 在这个函数内部，我们处理每一行的数据。
    # 'row' 是 kcat_list 的一行数据
    # 'foldchange' 是当前循环的变化系数
    # 'model', 'ecModel', 'tmp_data', 和 'rxn2block' 是其他必要的参数

    # 对每个目标基因、反应和kcat修改进行模拟
    with ecModel as ecmodel_tmp:
        target_gene = row['gene']
        kcat_m = row['kcat'] * foldchange
        rxnID = row['rxnID']
        ecmodel_tmp = updateEcGEMkcat(ecGEM=ecmodel_tmp, target_gene=target_gene, rxnID=rxnID, kcat_m=kcat_m)
        result = abc_python_max(model, ecmodel_tmp, tmp_data, 1, 1, rxn2block)

    # 返回结果和额外的信息
    return result[0][0][0], target_gene, rxnID

def parallel_task_preprocessing (j,growthdata_test,eModel,model,kcat_final,proteo_df,toy_data) :
    if j % 100==0:
        print (f'current set processed: {j}')

    index = 'kcat_value' + str (j)
    toy_data_tmp = pd.DataFrame (columns=['error' , 'corr' , 'ss'])
    for i in range (len (growthdata_test)):
        dilutionrate = growthdata_test[i , 0]
        with eModel as model_tmp:
            model_tmp = changeMedia (model , model_tmp , 'D-glucose' , "MIN")
            # model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
            model_tmp.reactions.get_by_id ('r_1631').bounds = 0 , 0
            model_tmp.reactions.get_by_id ('r_1761').upper_bound = 1000
            model_tmp.reactions.get_by_id ('r_1714').lower_bound = -growthdata_test[i , 1]  # D-glucose exchange
            model_tmp.reactions.get_by_id ("r_1654").bounds = -growthdata_test[i , -4] , 0  # NH4 exchange
            model_tmp.reactions.get_by_id ("prot_pool_exchange").upper_bound = growthdata_test[i , 9] * 1000
            for k in range (len (kcat_final)):
                target_gene0 = kcat_final['gene'][k]
                kcat_m0 = kcat_final[index][k]
                rxn0 = kcat_final['rxnID'][k]
                model_tmp = updateEcGEMkcat (ecModel=model_tmp ,target_gene=target_gene0 ,rxnID=rxn0 ,kcat_m=kcat_m0)
                # print(model_tmp.optimize().objective_value)

            # minimize the usage of protein pools
            model_tmp.objective = {model_tmp.reactions.prot_pool_exchange: 1}
            solution_f = model_tmp.optimize ()
            flux_max = solution_f.fluxes
            result = pd.DataFrame ({'rxnID': flux_max.index , 'flux': flux_max.values})
            result = result[result['rxnID'].str.contains ("prot_")]
            result['geneID'] = result['rxnID'].str.replace ("draw_prot_" , "")

            # compare with abundance
            abundance = pd.DataFrame ()
            column_names = [f'prot.{3 * i + 1}' , f'prot.{3 * i + 2}' , f'prot.{3 * i + 3}']
            column_name = f'abundance{i + 1}'
            abundance[column_name] = proteo_df[column_names].mean (axis=1)
            abundance1 = pd.concat ([abundance , proteo_df['emodel_genes']] , axis=1)
            result['pro_measured'] = singleMapping (abundance1[f"abundance{i + 1}"] ,abundance1["emodel_genes"] ,result['geneID'])
            result = result[~result["pro_measured"].isna ()]

            # pearsonr & pvalue
            from scipy.stats import pearsonr
            result1 = result[result['flux'] > 0]
            result1 = result1[result1['pro_measured'] > 0]
            corr , ss = pearsonr (np.log10 (result1['pro_measured']) , np.log10 (result1['flux']))
            # print("Correlation coefficient:", corr)
            # print("Correlation p_value:", ss)

            # pro_evaluation
            result1['calculated'] = abs (((result1['flux'] - result1['pro_measured']) / result1['pro_measured'])) ** 2
            error = np.sqrt (result1['calculated'].sum () / len (result1['calculated']))

            toy_data_tmp.loc[i] = [error , corr , ss]
    #print(toy_data_tmp)
    #toy_data.loc[j] = [toy_data_tmp['error'].mean () , toy_data_tmp['corr'].mean () , toy_data_tmp['ss'].mean ()]
    return toy_data_tmp['error'].mean(), toy_data_tmp['corr'].mean(), toy_data_tmp['ss'].mean()





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

