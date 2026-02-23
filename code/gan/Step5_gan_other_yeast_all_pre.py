from src.ec_N_lim_sen import *
from src.code_generation_preparation import getrSample
from src.constrain_ecmodel import updateprior, construct_MWs, construct_kcat_dict, constrainPool, \
    convertToEnzymeModel, get_r_max
import scipy.io as scio
from multiprocessing import Pool, cpu_count
from cobra.io import load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import numpy as np
import pandas as pd
import h5py
import time
import os
import pickle

def main() :
    strains = {'kma':325,'yli':178.5,'kla':245}
    # parameter determination
    for strain,prot in strains.items() :
        print(strain+" start!")
        number_of_models = 1000  # No of kcat samples
        enzyme_preparation = pd.read_csv (f'data/other_yeast/kcat_all_{strain}.csv')
        rxnlist = enzyme_preparation.drop ([ 'kcat_value1' ] , axis=1 , inplace=False)
        parameter_set_dim = len (enzyme_preparation)  # No of kcat parameters in model

        path = f"data/other_yeast/model_{strain}.xml"
        model_cobra = read_sbml_model (path)
        enzymedataFile = f"data/other_yeast/{strain}_dl.mat"
        z = scio.loadmat (enzymedataFile)
        model = z[ 'model' ][ 0 , 0 ]
        enzymedata = z[ 'enzymedata' ][ 0 , 0 ]
        rxn2block = z[ 'rxn2block' ]
        max_growth_train = z[ 'max_growth' ]
        growthdata_train = z[ 'growthdata' ]
        # strain = z[ 'strain' ]
        kcat_initial = pd.read_csv (f"data/other_yeast/{strain}_kcat_genra.txt" , delimiter="," , header=None)
        kcat_initial = np.array (kcat_initial)
        kcat_100 = kcat_initial[ 0:-2 , : ]
        [ a , b ] = updateprior (kcat_100)
        enzymedata[ 'kcat' ] = np.transpose (a)
        enzymedata[ 'kcat_var' ] = np.transpose (b)

        #replace @ by _
        for gene in model_cobra.genes:
            gene.id = gene.id.replace ('@' , '_')

        for reaction in model_cobra.reactions:
            if '@' in reaction.gene_reaction_rule:
                original_rule = reaction.gene_reaction_rule
                fixed_rule = original_rule.replace ('@' , '_')
                reaction.gene_reaction_rule = fixed_rule

        unique_genes = { gene.id : gene for gene in model_cobra.genes }
        model_cobra.genes = list ( unique_genes.values () )

        for i in range (enzymedata[ 'proteins' ].shape[ 0 ]):
            for j in range (enzymedata[ 'proteins' ][ i ].shape[ 0 ]):
                enzymedata[ 'proteins' ][ i ][ j ][ 0 ] = enzymedata[ 'proteins' ][ i ][ j ][ 0 ].replace ('@' , '_')

        for i in range (enzymedata[ 'enzyme' ].shape[ 0 ]):
            for j in range (enzymedata[ 'enzyme' ][ i ].shape[ 0 ]):
                parts = enzymedata[ 'enzyme' ][ i ][ j ][ 0 ].split (' and ')
                for k in range (len (parts)):
                    parts[ k ] = parts[ k ].replace ('@' , '_')
                enzymedata[ 'enzyme' ][ i ][ j ][ 0 ] = ' and '.join (parts)

        for i in range (enzymedata[ 'subunit' ].shape[ 0 ]):
            for j in range (enzymedata[ 'subunit' ].shape[ 1 ]):
                if enzymedata[ 'subunit' ][ i , j ].size > 0:
                    enzymedata[ 'subunit' ][ i , j ][ 0 ] = enzymedata[ 'subunit' ][ i , j ][ 0 ].replace ('@' , '_')


        # rxnlist = [ i[ 0 ][ 0 ] for i in enzymedata[ 'rxn_list' ] ]
        # savepath = 'data/other_yeast'
        # file_path = os.path.join (savepath , f'rxnlist_{strain}.pkl')
        # with open (file_path , 'wb') as file:
        #     pickle.dump (rxnlist , file)
        kcat_dict = construct_kcat_dict (model , enzymedata)
        kcat_df = pd.DataFrame (list (kcat_dict.items ()) , columns=[ 'gene_rxnID' , 'kcat_value' ])
        kcat_df[ [ 'gene' , 'rxnID' ] ] = pd.DataFrame (kcat_df[ 'gene_rxnID' ].tolist () , index=kcat_df.index)
        MWs = construct_MWs (enzymedata)
        kcat_dict = construct_kcat_dict (model , enzymedata)
        kcat_df = pd.DataFrame (list (kcat_dict.items ()) , columns=[ 'gene_rxnID' , 'kcat_value' ])
        kcat_df[ [ 'gene' , 'rxnID' ] ] = pd.DataFrame (kcat_df[ 'gene_rxnID' ].tolist () , index=kcat_df.index)
        eModel = convertToEnzymeModel (model_cobra , kcat_dict)
        eModel = constrainPool (eModel , MWs , { } , eModel.enzymes , prot)
      #  s2 = get_r_max (eModel , model_cobra)

        merged_df = pd.merge (enzyme_preparation , kcat_df , on=[ 'gene' , 'rxnID' ])
        kcat_final = merged_df[ [ 'gene' , 'rxnID' , 'kcat_value' ] ]
        kcat_all = kcat_df.drop ([ 'gene_rxnID' ] , axis=1 , inplace=False)
        # prepare the whole preprocessing dataframe
        rxn_list = [ ]
        for i in range (len (enzymedata[ 'kcat_var' ])):
            rxn_list.append (enzymedata[ 'rxn_list' ][ i ][ 0 ][ 0 ])
        rxn_df = pd.DataFrame( rxn_list,columns=[ 'rxn_list' ] )
        kcat_var_df = pd.DataFrame (enzymedata[ 'kcat_var' ] , columns=[ 'kcat_var' ])
        kcat_range = pd.DataFrame (enzymedata[ 'enzyme_ec_kcat_range' ] , columns=[ 'kcat_range_lb' , 'kcat_range_ub' ])
        enzymedata_pre = pd.concat ([ rxn_df , kcat_var_df , kcat_range ] , axis=1)
        kcat_pre = pd.merge (kcat_final , enzymedata_pre , left_on='rxnID' , right_on='rxn_list' , how='left')
        kcat_pre.drop (columns=[ 'rxn_list' ] , inplace=True)
        kcat_random_all = getrSample (kcat_pre[ 'kcat_value' ] ,
                                      kcat_pre[ 'kcat_var' ] ,
                                      kcat_pre[ 'kcat_range_lb' ] ,
                                      kcat_pre[ 'kcat_range_ub' ] ,
                                      number_of_models ,
                                      method='normal')
        kcat_df = pd.DataFrame (kcat_random_all)

        # kcat_samples
        column_names = [ f'kcat_value{i}' for i in range (kcat_df.shape[ 1 ]) ]
        kcat_df.columns = column_names
        kcat_final = pd.concat ([ kcat_final , kcat_df ] , axis=1)
        kcat_final.drop ([ 'kcat_value' ] , axis=1 , inplace=True)

        # eModel_genes = [ gene.id for gene in eModel.genes ]
        # eModel_genes_series = pd.read_csv('G:/my_github/EnzymeTuning/data/I.ori/dataset/gene_list.csv')
        # eModel_genes_series = eModel_genes_series.drop_duplicates(subset='Entry')
        proteomics = pd.read_csv (f'data/other_yeast/proteome_{strain}.csv' )
        # growthdata = pd.read_excel ('data/iIsor/growthdata.xlsx' , sheet_name="train")
        # growthdata = np.array (growthdata)

        ex_mets = ['biomass pseudoreaction' , 'D-glucose exchange' , 'acetate exchange' , 'ethanol exchange' ,
                       'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
                       'oxygen exchange' , 'prot_pool_exchange']

        # find the related rxnID
        idx = []
        for name0 in ex_mets:
            s = getRxnByReactionName(model=eModel, name=name0)
            if len(s) > 1:
                print("need check")
            elif len(s) == 1:
                idx.append(s[0])

        savepath = "data/other_yeast"
        os.makedirs (savepath , exist_ok=True)
        h5_file_path = f'data/other_yeast/kcat_data_all_{strain}.h5'
        with h5py.File (h5_file_path , 'w') as h5f:
            for i in range (kcat_random_all.shape[ 1 ]):
                h5f.create_dataset (f'kcat{i}' , data=kcat_random_all[ : , i ])

        # save the parameter names
        parameter = []
        for i in range (len (kcat_final)):
            parameter.append (kcat_final['gene'][i] + '_' + kcat_final['rxnID'][i])
        if not os.path.exists (savepath):
            os.makedirs (savepath)
        file_path = os.path.join (savepath , f'parameter_{strain}.pkl')
        with open (file_path , 'wb') as file:
            pickle.dump (parameter , file)

        # calculate the evalution parameter
        start = time.time ()
        print ('\nSTARTING PREPROCESSING')

        train_type = "protein"
        # train_type = "phenotype"

        # toy_data = pd.DataFrame (columns=['error', 'corr', 'pvalue','number','RMSE'])

        num_cpus = 16
        print (num_cpus)
        eModel.objective = "r_2111"
        eModel.objective_direction = 'max'
        with Pool (processes=num_cpus) as pool :
            args = [(
                j, growthdata_train, max_growth_train, eModel, model, kcat_final, proteomics,
                train_type, rxn2block,strain,prot) for
                j in range (number_of_models)]
            result = pool.starmap (parallel_task_other_yeast, args)
            print (result)
        results = np.vstack (result)
        toy_data = pd.DataFrame (results, columns=['error', 'corr', 'pvalue','number','error_num','RMSE'])

        toy_data.to_csv (f'data/other_yeast/pro_data_{strain}.csv', index=False)

        end = time.time ()
        print (f'PROCESSING DONE in {end - start:.05} seconds')


if __name__=='__main__' :
    main ()
