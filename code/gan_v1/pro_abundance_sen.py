from cobra.io import load_matlab_model, read_sbml_model
from cobra import Reaction, Metabolite
import sys
import os
path = "F:\python\Bayesian_python_ecmodel"
os.chdir(path)
import matplotlib.pyplot as plt
import seaborn as sns
import cobra
import scipy.io as scio
# import self function
from src.mainFunction import *
from src.protein_process import *
from src.model_process import *
from src.constrain_ecmodel import *

#load the dl_model
enzymedataFile = "data/model/Saccharomyces_cerevisiae_dl.mat"
z = scio.loadmat(enzymedataFile)
model = z['model'][0,0]
enzymedata = z['enzymedata'][0,0]
#scio.savemat('data/model/model.mat', {'model': model})
model_cobra = load_matlab_model('data/model/model.mat')

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
    s2 = get_r_max(ecModel,model_cobra)
    print("growth rate: ", s2.objective_value)
    print("Protein_pool: ", ecModel.reactions.get_by_id("prot_pool_exchange").upper_bound)

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

# 创建DataFrame
kcat_list = pd.DataFrame({
    'gene': genes,
    'rxnID': rxnIDs,
    'kcat': kcat_values
})

dilutionrate = 0.42
correlation_coef = {}
foldchanges = [0.1, 0.9, 1.1, 1.5]
for foldchange in foldchanges:
    Correlation_coefficient = []
    for index, row in kcat_list.iterrows():
        with ecModel as model_tmp:
            #change each kcat in each reaction
            target_gene0 = row['gene']
            kcat_m0 = row['kcat']*foldchange
            rxn0 = row['rxnID']
            model_tmp = updateEcGEMkcat(ecGEM=model_tmp, target_gene=target_gene0, rxnID=rxn0, kcat_m=kcat_m0)

            #phenotype in PNAS2021
            model_tmp.reactions.r_1714.lower_bound = -20
            model_tmp.reactions.r_1634.lower_bound = 0.566
            model_tmp.reactions.r_1761.lower_bound = 27
            model_tmp.reactions.r_1808.lower_bound = 1.469
            print(model_tmp.optimize().objective_value)
            #model_tmp.reactions.get_by_id(idx[1]).lower_bound = -1000  # glucose uptake
            model_tmp.reactions.get_by_id(idx[0]).lower_bound = dilutionrate
            model_tmp.objective = {model_tmp.reactions.r_1714: 1}  # minimize the uptake of glucose
            solution2 = model_tmp.optimize()
            print(solution2.objective_value)
            try:
            # then fix glucose uptake and minimize the protein pool
                model_tmp.reactions.get_by_id(idx[1]).lower_bound = solution2.objective_value * 1.00001

                print('Glucose uptake rate: ',solution2.objective_value)
            except TypeError:
                Correlation_coefficient.append(np.nan)
                continue
            model_tmp.reactions.get_by_id(idx[9]).lower_bound = -1000 # protein pool
            model_tmp.objective = {model_tmp.reactions.prot_pool_exchange: 1}  # minimize the usage of protein pools
            solution_f = model_tmp.optimize()
                #
            flux_max = solution_f.fluxes
            result = pd.DataFrame({'rxnID':flux_max.index, 'flux':flux_max.values})
            result = result[result['rxnID'].str.contains("prot_")]
            result['geneID'] = result['rxnID'].str.replace("draw_prot_", "")
            # input the proteomics under max growth rate
            abundance_ex = pd.read_excel("data/proteomics/data_PNAS_2021.xlsx")
            abundance_ex['g/gDW'] =(abundance_ex['replicate 1 (g gDW-1)']+ abundance_ex['replicate 2 (g gDW-1)']+ abundance_ex['replicate 3 (g gDW-1)'])/3
            abundance_ex=abundance_ex[['Symbol','g/gDW']]
            abundance_ex.columns = ['gene','g/gDW']
            abundance_ex1 = splitAbundance(pro_df=abundance_ex)
                # change the unit from g/gDW as mmol/gDW
                # input the molecular weight
            mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
            mw = mw[["locus","proteins_molecular_weight"]]
            mw.columns = ["gene name", "MW"]
            mw["MW_Kda"] = mw["MW"]/1000
            abundance_ex1["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], abundance_ex1["gene"])
            abundance_ex_check = abundance_ex1[abundance_ex1["MW_Kda"].isna()]
            abundance_ex1=abundance_ex1[~abundance_ex1["MW_Kda"].isna()]
            abundance_ex1["mmol/gDW"] = abundance_ex1["g/gDW"]/abundance_ex1["MW_Kda"]# #mmol/g biomass

            result['pro_measured'] = singleMapping(abundance_ex1["mmol/gDW"],abundance_ex1["gene"],result['geneID'])
            result = result[~result["pro_measured"].isna()]
            result.to_excel("data/data_check.xlsx")

                # change the protein abundance unit from mmol/gDW into protein copy/cell
            coefficient1 = 7.8298e9
            result_unify = result.copy()
            result_unify["pro_measured"] = result['pro_measured']*coefficient1
            result_unify["flux"] = result['flux']*coefficient1

            from scipy.stats import pearsonr
            result1 = result[result['flux'] > 0]
            result1 = result1[result1['pro_measured'] > 0]
            corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
            Correlation_coefficient.append(corr)
            print("Correlation coefficient:", corr)
            print("Correlation p_value:", ss)

    correlation_coef_df = pd.DataFrame({
        'gene': genes,
        'rxnID': rxnIDs,
        'Correlation_coefficient': Correlation_coefficient
    })
    correlation_coef["correlation_coef_{0}".format(foldchange)] = correlation_coef_df


for key, df in correlation_coef.items():
    filename = key + ".csv"
    df.to_csv(filename, index=False)
