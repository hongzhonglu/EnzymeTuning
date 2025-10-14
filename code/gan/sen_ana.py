
import numpy as np
from cobra.io import load_matlab_model, read_sbml_model
from cobra import Reaction, Metabolite
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import scipy.io as scio
# import self function
from src.ec_N_lim_sen import *
from src.mainFunction import *
from src.protein_process import *
from src.model_process import *
from src.constrain_ecmodel import *
from SALib.sample import sobol
from SALib.analyze import sobol as analyze_sobol
from multiprocessing import Pool

def main() :
    #load the dl_model
    enzymedataFile = "data/model/Saccharomyces_cerevisiae_dl.mat"
    z = scio.loadmat(enzymedataFile)
    model = z['model'][0,0]
    enzymedata = z['enzymedata'][0,0]
    rxn2block= z['rxn2block']
    #scio.savemat('data/model/model.mat', {'model': model})
    model_cobra = load_matlab_model('data/model/model.mat')
    model_cobra.solver = "cplex"

    #construnct ecModel
    #kcat_tmp = pd.read_csv("data/kcat_genra4.txt",delimiter=",",header=None)
    kcat_tmp = pd.read_csv("data/kcat_genra100.txt",delimiter=",",header=None)
    kcat_tmp = np.array(kcat_tmp)
    kcat_100 = kcat_tmp[0:-2,:]
    [a,b] = updateprior(kcat_100)
    enzymedata['kcat'] = np.transpose(a)
    enzymedata['kcat_var'] = np.transpose(b)
    kcat_dict = construct_kcat_dict(model,enzymedata)
    eModel = convertToEnzymeModel(model_cobra,kcat_dict)
    MWs = construct_MWs(enzymedata)

    #without proteomic data (initial version)
    with eModel:
        ecModel = constrainPool(eModel, MWs, {}, eModel.enzymes, 230)
       # s2 = get_r_max(ecModel,model_cobra)
       # print("growth rate: ", s2.objective_value)
       # print("Protein_pool: ", ecModel.reactions.get_by_id("prot_pool_exchange").upper_bound)

    #exchange reactions
    ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
                'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
                'oxygen exchange', 'prot_pool_exchange']
    # find the related rxnID
    idx = []
    for name0 in ex_mets:
        #print(name0)
        s = getRxnByReactionName(model=ecModel, name=name0)
        if len(s) > 1:
            print("need check")
        elif len(s) == 1:
            idx.append(s[0])
    genes, rxnIDs = zip(*kcat_dict.keys())
    kcat_values = list(kcat_dict.values())

    # create DataFrame
    kcat_list = pd.DataFrame({
        'gene': genes,
        'rxnID': rxnIDs,
        'kcat': kcat_values
    })
    kcat_list['name'] = kcat_list['gene'] + '+' + kcat_list['rxnID']
    kcat_list = kcat_list.iloc[0:10,:]

    #growth data
    growth_N = pd.read_csv('data/grwothdata_Nlim.csv')
    growthdata_N = growth_N.drop(columns=['qNitrogen (mmol/gDW h)'])
    growthdata_N = np.array(growthdata_N)
    growth_C = pd.read_excel('data/Xia.xlsx',sheet_name='Sheet1')
    growthdata_C = np.array(growth_C)

    #proteomics data (Clim & Nlim)
    eModel_genes = [gene.id for gene in ecModel.genes]
    eModel_genes_series = pd.DataFrame(eModel_genes,columns=['emodel_genes'])
    proteomics_Clim = pd.read_csv('data/xia.csv')
    proteomics_Clim = proteomics_Clim.dropna()
    proteomicsdata_Clim = proteomics_Clim[proteomics_Clim["all_gene"].isin(eModel_genes_series["emodel_genes"])]
    proteomics_Nlim = pd.read_csv('data/Nlim.csv')
    proteomics_Nlim = proteomics_Nlim.dropna()
    proteomicsdata_Nlim = proteomics_Nlim[proteomics_Nlim["all_gene"].isin(eModel_genes_series["emodel_genes"])]

    #define the problem (kcat)
    num_kcats = len(kcat_list)
    print(num_kcats)
    problem = {
        'num_vars' : num_kcats,
        'names' : kcat_list['name'].tolist(),
        'bounds' : [[0.1*kcat,2*kcat] for kcat in kcat_list['kcat']]
    }

    #sampling
    param_values = sobol.sample(problem,2) #num_kcats*2*16 (2^n)

    #run the model
    num_cpus = 10
    pool = Pool(processes=num_cpus)

    results = pool.starmap(model_execution, [(params, kcat_list, model, ecModel, growthdata_N, rxn2block) for params in param_values])
    pool.close()
    pool.join()
    print(results)

    RMSE , genes , rxnIDs = zip (*results)
    RMSE_df = pd.DataFrame ({
        'gene': genes ,
        'rxnID': rxnIDs ,
        'RMSE': RMSE
    })
    #RMSE_list["RMSE{0}".format (foldchange)] = RMSE_df
    RMSE_df.to_csv ("result/rmse2_10.csv" , index=False)
    print(RMSE_df)
    Y = np.array(RMSE)

    Si = analyze_sobol.analyze(problem, Y, print_to_console=True)
    print(Si)

    # save the result
    df_Si = pd.DataFrame({
        'ID': problem['names'],
        'S1': Si['S1'],
        'ST': Si['ST']
    })
    df_Si.to_csv("result/sensitivity_analysis2_100.csv", index=False)
    print(df_Si)

if __name__ == '__main__':
    main()
