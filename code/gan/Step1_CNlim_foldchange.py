from src.ec_N_lim_sen import parallel_task_preprocessing , getRxnByReactionName
from src.code_generation_preparation import getrSample
from src.constrain_ecmodel import updateprior , construct_MWs , construct_kcat_dict , constrainPool , \
    convertToEnzymeModel , get_r_max
import scipy.io as scio
from multiprocessing import Pool , cpu_count
from cobra.io import load_matlab_model , save_matlab_model , read_sbml_model , write_sbml_model
import numpy as np
import pandas as pd
import h5py
import time
import os
import pickle

def main():

    # load the enzyme_id
    #data = pd.read_csv ('data/RMSE2.csv')
    #enzyme_preparation = data[data['ERROR'] < 0.05]
    enzyme_preparation = pd.read_csv ('data/kcat_merge_error_50.csv')
    #enzyme_preparation = pd.read_csv ('data/kcat_merge_sen_50.csv')
    #enzyme_preparation = pd.read_csv ('data/kcat_merge_error_100.csv')
    #enzyme_preparation = pd.read_csv ('data/kcat_merge_sen_100.csv')

    # parameter determination
    number_of_models = 10000  # No of kcat samples
    #parameter_set_dim = 1151  # No of kcat parameters in model
    parameter_set_dim = len(enzyme_preparation)  # No of kcat parameters in model

    # load enzyme information
    enzymedataFile = "data/Saccharomyces_cerevisiae_dl.mat"
    z = scio.loadmat (enzymedataFile)
    model = z['model'][0 , 0]
    enzymedata = z['enzymedata'][0 , 0]
    rxn2block = z['rxn2block']
    kcat_initial = pd.read_csv ("data/kcat_genra100.txt" ,
                                delimiter="," ,
                                header=None)
    kcat_initial = np.array (kcat_initial)
    kcat_100 = kcat_initial[0:-2 , :]
    [a , b] = updateprior (kcat_100)
    enzymedata['kcat'] = np.transpose (a)
    enzymedata['kcat_var'] = np.transpose (b)

    # construct kcat_list
    MWs = construct_MWs (enzymedata)
    kcat_dict = construct_kcat_dict (model , enzymedata)
    kcat_df = pd.DataFrame (list (kcat_dict.items ()) , columns=['gene_rxnID' , 'kcat_value'])
    kcat_df[['gene' , 'rxnID']] = pd.DataFrame (kcat_df['gene_rxnID'].tolist () , index=kcat_df.index)
    merged_df = pd.merge (enzyme_preparation , kcat_df , on=['gene' , 'rxnID'])
    kcat_final = merged_df[['gene' , 'rxnID' , 'kcat_value']]

    # load ecModel
    model_cobra = load_matlab_model ('data/model/model.mat')
    eModel = convertToEnzymeModel (model_cobra , kcat_dict)
    eModel = constrainPool (eModel , MWs , {} , eModel.enzymes , 230)
    s2 = get_r_max (eModel , model_cobra)

    # proteomics data (Clim & Nlim)
    eModel_genes = [gene.id for gene in eModel.genes]
    eModel_genes_series = pd.DataFrame (eModel_genes , columns=['emodel_genes'])
    proteomics_Clim = pd.read_csv ('data/Clim_all.csv')
    proteomicsdata_Clim = proteomics_Clim[proteomics_Clim["all_gene"].isin (eModel_genes_series["emodel_genes"])]
    proteomics_Nlim = pd.read_csv ('data/Nlim_train.csv')
    proteomicsdata_Nlim = proteomics_Nlim[proteomics_Nlim["all_gene"].isin (eModel_genes_series["emodel_genes"])]

    # growth data
    growthdata = pd.read_excel ('data/growthdata.xlsx' , sheet_name="Sheet1")
    growthdata = np.array (growthdata)  # 7N+9C



    # update the kcat
    ex_mets = ['biomass pseudoreaction' , 'D-glucose exchange' , 'acetate exchange' , 'ethanol exchange' ,
               'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
               'oxygen exchange' , 'prot_pool_exchange']
    # find the related rxnID
    idx = []
    for name0 in ex_mets:
        s = getRxnByReactionName (model=eModel , name=name0)
        if len (s) > 1:
            print ("need check")
        elif len (s)==1:
            idx.append (s[0])

    # prepare the whole preprocessing dataframe
    rxn_list = []
    for i in range (len (enzymedata['kcat_var'])):
        rxn_list.append (enzymedata['rxn_list'][i][0][0])
    rxn_df = pd.DataFrame (rxn_list , columns=['rxn_list'])
    kcat_var_df = pd.DataFrame (enzymedata['kcat_var'] , columns=['kcat_var'])
    kcat_range = pd.DataFrame (enzymedata['enzyme_ec_kcat_range'] , columns=['kcat_range_lb' , 'kcat_range_ub'])
    enzymedata_pre = pd.concat ([rxn_df , kcat_var_df , kcat_range] , axis=1)
    kcat_pre = pd.merge (kcat_final , enzymedata_pre , left_on='rxnID' , right_on='rxn_list' , how='left')
    kcat_pre.drop (columns=['rxn_list'] , inplace=True)
    kcat_random_all = getrSample (kcat_pre['kcat_value'] ,
                                  kcat_pre['kcat_var'] ,
                                  kcat_pre['kcat_range_lb'] ,
                                  kcat_pre['kcat_range_ub'] ,
                                  number_of_models ,
                                  method='normal')
    kcat_df = pd.DataFrame (kcat_random_all)

    # kcat_samples
    column_names = [f'kcat_value{i}' for i in range (kcat_df.shape[1])]
    kcat_df.columns = column_names
    kcat_final = pd.concat ([kcat_final , kcat_df] , axis=1)
    kcat_final.drop (['kcat_value'] , axis=1 , inplace=True)

    h5_path = 'result/gan'
    if not os.path.exists (h5_path):
        os.makedirs (h5_path)
    h5_file_path = 'result/gan/kcat_data.h5'
    with h5py.File (h5_file_path , 'w') as h5f:
        for i in range (kcat_random_all.shape[1]):
            h5f.create_dataset (f'kcat{i}' , data=kcat_random_all[: , i])

    # save the parameter names
    parameter = []
    for i in range (len (kcat_final)):
        parameter.append (kcat_final['gene'][i] + '_' + kcat_final['rxnID'][i])
    savepath = 'gan'
    if not os.path.exists (savepath):
        os.makedirs (savepath)
    file_path = os.path.join (savepath , 'parameter.pkl')
    with open (file_path , 'wb') as file:
        pickle.dump (parameter , file)

    # calculate the evalution parameter
    start = time.time ()
    print ('\nSTARTING PREPROCESSING')

    toy_data = pd.DataFrame (columns=['error' , 'corr' , 'ss'])

    num_cpus = 128
    print (num_cpus)
    with Pool (processes=num_cpus) as pool:
        args = [(
            j , growthdata , eModel , model , kcat_final , proteomicsdata_Clim, proteomicsdata_Nlim , toy_data) for j in range (number_of_models)]
        result = pool.starmap (parallel_task_preprocessing , args)
        print(result)
    results = np.vstack(result)
    toy_data = pd.DataFrame(results,columns=['error' , 'corr' , 'ss'])

    toy_data.to_csv ('result/gan/pro_data.csv' , index=False)

    end = time.time ()
    print (f'PROCESSING DONE in {end - start:.05} seconds')

if __name__ == '__main__':
    main()