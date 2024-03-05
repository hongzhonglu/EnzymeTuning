import numpy as np
import os
import cobra
from pathlib import Path
from sklearn.metrics import mean_squared_error
import pandas as pd
from cobra.io import load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import scipy.io as scio
from scipy import sparse
import cobra.util
from cobra import Model,Reaction,Metabolite

from code0824 import abc_python_max
from code_generation_preparation import getrSample
from code_generation_preparation import sumBiomass
from code_generation_preparation import updateprior

def change_rxn_coeff(rxn,met,new_coeff):
    '''
    # This is based on the rxn.add_metabolites function. If there the metabolite is already in the reaction,
    # new and old coefficients will be added. For example, if the old coeff of metA is 1, use
    # rxn.add_metabolites({metA:2}), After adding, the coeff of metA is 1+2 = 3
    #
    '''

    diff_coeff = new_coeff-rxn.metabolites[met]
    rxn.add_metabolites({met:diff_coeff})

species = "Saccharomyces_cerevisiae"
generation = 150
enzymedataFile = "data/Saccharomyces_cerevisiae_dl.mat"
#enzymedataFile = "data/model.mat"
z = scio.loadmat(enzymedataFile)
enzymedata = z['enzymedata'][0,0]
max_growth = z['max_growth']
growthdata = z['growthdata']
model = z['model'][0,0]
strain = z['strain']
rxn2block = z['rxn2block']

#load the model in cobra format
scio.savemat('data/model.mat', {'model': model})
model_cobra = load_matlab_model('data/model.mat')


import scipy.io as io
F = sumBiomass(model,model_cobra)
#calculate the prot_weight
#F = 0.46  #怎么从sumBiomass中导出结果，matlab版本为0.46
tot_prot_weight = F[0]*0.5
if strain == 'Kluyveromyces_marxianus':
    tot_prot_weight = 0.325
elif strain == 'Kluyveromyces_lactis':
    tot_prot_weight = 0.245

met = cobra.Metabolite ('cost' , name='cost' , compartment='c')
model_cobra.add_metabolites (met)
model_cobra.add_boundary (model_cobra.metabolites.get_by_id ("cost") , type="sink" , ub=tot_prot_weight * 1000 , lb=0)

import time
kcat_random_all = enzymedata["kcat"]
j = 1
proc =1
sample_generation = 1
kcat_random = kcat_random_all
prot_cost_info_value = [enzymedata["MW"][i] / kcat_random[i] for i in range (len (kcat_random))]
prot_cost_info_value = [item[0] for item in prot_cost_info_value]
prot_cost_info_id = [str(item[0][0]) for item in enzymedata["rxn_list"]]
prot_cost_info = dict (zip (prot_cost_info_id , prot_cost_info_value))
nMets = len (model_cobra.metabolites)
nRxns = len (model_cobra.reactions)
cost_list = np.zeros ((nRxns , 1))
for p , rxnid in enumerate (model['rxns']):
    rxnid = tuple (rxnid)
    cost = prot_cost_info.get (rxnid[0][0] , 0)
    cost_list[p , 0] = cost
            # cost_list_new = np.transpose (cost_list)
            # cost_list = [prot_cost_info_value[np.where (prot_cost_info_id==rxnid)[0][0]] if rxnid in prot_cost_info_id else 0 for rxnid in model['rxns']]
            # met = cobra.Metabolite ('cost' , name='cost' , compartment='c')
            # model_tmp.add_metabolites (met)

start_time3 = time.time ()
            #tupll = tuple ([float (m) for m in cost_list])
cost_floats = [float (m) for m in cost_list]
            #reactions_with_cost = [r for r in model_tmp.reactions if 'cost' in r.metabolites]
            #for q , reaction in enumerate (reactions_with_cost):
            #    coeff = reaction.get_coefficient ('cost')
            #    reaction.add_metabolites ({'cost': cost_floats[q] - coeff})
for q , reaction in enumerate (model_cobra.reactions):
    if 'cost' in reaction.metabolites:
        coeff = reaction.get_coefficient ('cost')
    else:
        coeff = 0
    reaction.add_metabolites ({'cost': cost_floats[q] - coeff})
            # model_tmp.add_boundary (model_tmp.metabolites.get_by_id ("cost") , type="sink" , ub=tot_prot_weight * 1000 ,
            #                        lb=0)

end_time3 = time.time ()
execution_time3 = end_time3 - start_time3
print ("execution time1: " , execution_time3 , " seconds")

start_time3 = time.time ()
for q , reaction in enumerate (model_cobra.reactions):
    if 'cost' in reaction.metabolites:
        reaction.add_metabolites ({'cost': cost_floats[q] - reaction.get_coefficient ('cost')})
    else:
        reaction.add_metabolites ({'cost': cost_floats[q]})

end_time3 = time.time ()
execution_time3 = end_time3 - start_time3
print ("execution time2: " , execution_time3 , " seconds")

#model_cobra = load_matlab_model('data/model.mat')
#met = cobra.Metabolite ('cost' , name='cost' , compartment='c')
#model_cobra.add_metabolites (met)
#model_cobra.add_boundary (model_cobra.metabolites.get_by_id ("cost") , type="sink" , ub=tot_prot_weight * 1000 , lb=0)

start_time3 = time.time ()
cost_tmp = model_cobra.metabolites.cost
for q in range(len(model_cobra.reactions)):
    try:
        change_rxn_coeff(model_cobra.reactions[q],cost_tmp,cost_floats[q])
    except KeyError:
        continue
end_time3 = time.time ()
execution_time3 = end_time3 - start_time3
print ("execution time3: " , execution_time3 , " seconds")

met = model_cobra.metabolites.get_by_id("cost")

start_time3 = time.time ()
filtered_reactions = [reaction for reaction in model_cobra.reactions if met in reaction.metabolites]
q = 0
for reaction in filtered_reactions:
    cost_value = prot_cost_info_value[q]
    coeff = reaction.get_coefficient('cost')
    #reaction.add_metabolites({"cost":cost_value-coeff})
    q = q+1
end_time3 = time.time ()
execution_time3 = end_time3 - start_time3
print ("execution time4: " , execution_time3 , " seconds")

start_time3 = time.time ()
filtered_reactions = [reaction for reaction in model_cobra.reactions if met in reaction.metabolites]
q = 0
for q in range(len(filtered_reactions)):
    cost_value = prot_cost_info_value[q]
    coeff = filtered_reactions[q].get_coefficient('cost')
    reaction.add_metabolites({"cost":cost_value-coeff})
    #q = q+1
end_time3 = time.time ()
execution_time3 = end_time3 - start_time3
print ("execution time5: " , execution_time3 , " seconds")

