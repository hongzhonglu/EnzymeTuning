from src.ec_N_lim_sen import *
from src.constrain_ecmodel import *
from src.code_generation_preparation import *
import numpy as np
import pandas as pd
import h5py
import time
import os
import pickle
path = "F:\python\Bayesian_python_ecmodel"
os.chdir(path)

#parameter determination
number_of_models = 10000  # No of kcat samples
parameter_set_dim = 40  # No of kcat parameters in model

#load the enzyme_id
enzyme_preparation = pd.read_csv('kcat_merge_2.csv')

#load enzyme information
enzymedataFile = "data/Saccharomyces_cerevisiae_dl.mat"
z = scio.loadmat(enzymedataFile)
model = z['model'][0,0]
enzymedata = z['enzymedata'][0,0]
rxn2block= z['rxn2block']
kcat_initial= pd.read_csv("data/kcat_genra100.txt",delimiter=",",header=None)
kcat_initial = np.array(kcat_initial)
kcat_100 = kcat_initial[0:-2,:]
[a,b] = updateprior(kcat_100)
enzymedata['kcat'] = np.transpose(a)
enzymedata['kcat_var'] = np.transpose(b)

#construct kcat_list
MWs = construct_MWs(enzymedata)
kcat_dict = construct_kcat_dict(model,enzymedata)
kcat_df = pd.DataFrame(list(kcat_dict.items()), columns=['gene_rxnID', 'kcat_value'])
kcat_df[['gene', 'rxnID']] = pd.DataFrame(kcat_df['gene_rxnID'].tolist(), index=kcat_df.index)
merged_df = pd.merge(enzyme_preparation, kcat_df, on=['gene', 'rxnID'])
kcat_final = merged_df[['gene', 'rxnID', 'kcat_value']]
#kcat_final = data = pd.read_csv('data/merged_data.csv')

#load ecModel
model_cobra= load_matlab_model('data/model/model.mat')
eModel = convertToEnzymeModel(model_cobra,kcat_dict)
eModel = constrainPool(eModel, MWs, {}, eModel.enzymes, 230)
s2 = get_r_max(eModel,model_cobra)

#prepare proteomic data
scer_data = pd.read_excel('data/uniprotkb_Scer.xlsx', sheet_name= "Sheet0")
proteomics_data = pd.read_excel('data/proteomics_Nlim_all.xlsx',sheet_name='Sheet1')
proteo_df = proteomics_data.merge(scer_data[['Entry', 'Gene Names']], left_on='Accession', right_on='Entry', how='left')
prot_cols = [col for col in proteo_df.columns if col.startswith('prot.')]
proteo_df[prot_cols] = proteo_df[prot_cols].apply(lambda x: x * 10**-9)
eModel_genes = [gene.id for gene in eModel.genes]

#growth data
growthdata = pd.read_csv('data/proteomics_test.csv')
growthdata = growthdata.drop(columns=['qNitrogen (mmol/gDW h)'])
growthdata = np.array(growthdata)
growthdata_test = growthdata[0:7] #NH4+

gene_to_entry_mapping = {}
# find the protein id to each gene names (one to several)
for index, row in scer_data.iterrows():
    entry = row['Entry']
    gene_names = row['Gene Names'].split()
    for gene in gene_names: gene_to_entry_mapping[gene] = entry

eModel_genes_series = pd.DataFrame(eModel_genes,columns=['emodel_genes'])
eModel_genes_series['entry_to_gene'] = eModel_genes_series['emodel_genes'].map(gene_to_entry_mapping)
proteo_df = proteo_df.merge(eModel_genes_series[['emodel_genes', 'entry_to_gene']], left_on='Accession', right_on='entry_to_gene', how='inner')
proteo_df.drop(['Accession','Gene','protein length','Entry','Gene Names'],axis = 1,inplace= True)

#update the kcat
ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
               'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
              'oxygen exchange', 'prot_pool_exchange']
# find the related rxnID
idx = []
for name0 in ex_mets:
    s = getRxnByReactionName(model=eModel, name=name0)
    if len(s) > 1:
        print("need check")
    elif len(s) == 1:
        idx.append(s[0])

#prepare the whole preprocessing dataframe
rxn_list = []
for i in range (len (enzymedata['kcat_var'])) :
    rxn_list.append(enzymedata['rxn_list'][i][0][0])
rxn_df = pd.DataFrame(rxn_list,columns=['rxn_list'])
kcat_var_df = pd.DataFrame(enzymedata['kcat_var'],columns=['kcat_var'])
kcat_range = pd.DataFrame(enzymedata['enzyme_ec_kcat_range'],columns=['kcat_range_lb','kcat_range_ub'])
enzymedata_pre = pd.concat([rxn_df, kcat_var_df, kcat_range], axis=1)
kcat_pre = pd.merge(kcat_final,enzymedata_pre,left_on='rxnID',right_on='rxn_list',how = 'left')
kcat_pre.drop(columns=['rxn_list'], inplace=True)
kcat_random_all = getrSample(kcat_pre['kcat_value'], kcat_pre['kcat_var'], kcat_pre['kcat_range_lb'],kcat_pre['kcat_range_ub'], 10000,method = 'normal')
kcat_df = pd.DataFrame(kcat_random_all)

# kcat_samples
column_names = [f'kcat_value{i}' for i in range(kcat_df.shape[1])]
kcat_df.columns = column_names
kcat_final = pd.concat([kcat_final,kcat_df],axis=1)
kcat_final.drop(['kcat_value'],axis=1,inplace=True)

h5_file_path = 'data/gan_input/kcat_data.h5'
with h5py.File(h5_file_path, 'w') as h5f:
    for i in range(kcat_random_all.shape[1]) :
        h5f.create_dataset(f'kcat{i}', data=kcat_random_all[:, i])

#save the parameter names
parameter = []
for i in range(len(kcat_final)) :
    parameter.append(kcat_final['gene'][i]+'_'+kcat_final['rxnID'][i])
savepath = 'data/gan_input'
if not os.path.exists(savepath):
    os.makedirs(savepath)
file_path = os.path.join(savepath, 'parameter.pkl')
with open(file_path, 'wb') as file:
    pickle.dump(parameter, file)

print(kcat_final)
#calculate the evalution parameter
start = time.time()
print('\nSTARTING PREPROCESSING')

toy_data = pd.DataFrame(columns=['error', 'corr', 'ss'])


for j in range(10000) :
    #if j % 100 == 0:
    print(f'current set processed: {j}')
    index = 'kcat_value'+str(j)
    with eModel as model_tmp:
        for k in range (len (kcat_final)):
            target_gene0 = kcat_final['gene'][k]
            kcat_m0 = kcat_final[index][k]
            print(kcat_m0)
            rxn0 = kcat_final['rxnID'][k]
            model_tmp = updateEcGEMkcat (ecGEM=model_tmp , target_gene=target_gene0 , rxnID=rxn0 , kcat_m=kcat_m0)
            # print(model_tmp.optimize().objective_value)


        toy_data_tmp = pd.DataFrame(columns=['error', 'corr', 'ss'])
        for i in range(len(growthdata_test)) :
            dilutionrate = growthdata_test[i,0]

            model_tmp = changeMedia (model , model_tmp , 'D-glucose' , "MIN")
            #model_tmp.reactions.get_by_id ('r_1634').bounds = 0 , 0
            model_tmp.reactions.get_by_id ('r_1631').bounds = 0 , 0
            model_tmp.reactions.get_by_id ('r_1761').upper_bound = 1000
            model_tmp.reactions.get_by_id ('r_1714').lower_bound = -growthdata_test[i , 1] #D-glucose exchange
            model_tmp.reactions.get_by_id("r_1654").bounds = -growthdata_test[i,-4],0 #NH4 exchange
            model_tmp.reactions.get_by_id("prot_pool_exchange").upper_bound = growthdata_test[i,9]*1000 #protein


            # minimize the usage of protein pools
            model_tmp.objective = {model_tmp.reactions.prot_pool_exchange: 1}
            solution_f = model_tmp.optimize()
            flux_max = solution_f.fluxes
            result = pd.DataFrame({'rxnID':flux_max.index, 'flux':flux_max.values})
            result = result[result['rxnID'].str.contains("prot_")]
            result['geneID'] = result['rxnID'].str.replace("draw_prot_", "")

            #compare with abundance
            abundance = pd.DataFrame()
            column_names = [f'prot.{3*i+1}',f'prot.{3*i+2}',f'prot.{3*i+3}']
            column_name = f'abundance{i+1}'
            abundance[column_name] = proteo_df[column_names].mean(axis=1)
            abundance1 = pd.concat([abundance,proteo_df['emodel_genes']],axis=1)
            result['pro_measured'] = singleMapping(abundance1[f"abundance{i+1}"],abundance1["emodel_genes"],result['geneID'])
            result = result[~result["pro_measured"].isna()]

            #pearsonr & pvalue
            from scipy.stats import pearsonr
            result1 = result[result['flux'] > 0]
            result1 = result1[result1['pro_measured'] > 0]
            corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
            print("Correlation coefficient:", corr)
            print("Correlation p_value:", ss)

                #pro_evaluation
            result1['calculated'] = (((result1['flux'] - result1['pro_measured']) / result1['pro_measured']))**2
            error = np.sqrt(result1['calculated'].sum()/len(result1['calculated']))
            error = np.log10(error)
            print("Correlation error:",error)

            toy_data_tmp.loc[i] = [error, corr, ss]
    toy_data.loc[j] = [toy_data_tmp['error'].mean(),toy_data_tmp['corr'].mean(),toy_data_tmp['ss'].mean()]
    print("Final index:",toy_data)


toy_data.to_csv('data/gan_input/pro_data.csv', index=False)

end = time.time()
print(f'PROCESSING DONE in {end - start:.05} seconds')


