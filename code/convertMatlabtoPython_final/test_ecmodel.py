import numpy as np
import os
import cobra
import pandas as pd
from cobra.io import load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import scipy.io as scio
from cobra import Model,Reaction,Metabolite
path = "F:\python\Bayesian_python"
os.chdir(path)
from constrain_ecmodel import parse_gr_rule
from constrain_ecmodel import addEnzymesToRxn
from constrain_ecmodel import construct_kcat_dict
from constrain_ecmodel import constrainPool
from constrain_ecmodel import get_r_max

enzymedataFile = "data/Saccharomyces_cerevisiae_dl.mat"
z = scio.loadmat(enzymedataFile)
enzymedata = z['enzymedata'][0,0]
max_growth = z['max_growth']
growthdata = z['growthdata']
model = z['model'][0,0]
strain = z['strain']
rxn2block = z['rxn2block']

model_cobra = load_matlab_model('data/model.mat')
model_dl = cobra.io.read_sbml_model("Saccharomyces_cerevisiae_dl.xml")
kcat_dict = construct_kcat_dict(model,enzymedata)

converted_reaction_list = []
protein_exchange_rxns = {}
number = 0
for rxn in model_cobra.reactions:
    complexes = parse_gr_rule(rxn.gene_reaction_rule)

        # 1. for those reactions without genes
    if len(complexes) <1:
        converted_reaction_list.append(rxn)
    #    #print(rxn,complexes)
        continue

    # 2. for those reactions with genes, but no kcat
    first_gene = [gene.id for gene in rxn.genes][0]
    if kcat_dict.get((first_gene,rxn.id),None) is None:
        converted_reaction_list.append(rxn)
        #print(rxn,complexes)
        continue

    # 3. for those reactions with isoenzymes, add arm reaction
    if len(complexes) >1:
        for i,complx in enumerate(complexes):
            prots = [item.strip() for item in complx.split('and')]
            kcat = kcat_dict[(prots[0],rxn.id)]
            e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn_new, kcat, complx,rxn_index=i+1)

            converted_reaction_list.append(e_rxn)
            for prot_exchange_rxn in prot_exchange_rxns: protein_exchange_rxns[prot_exchange_rxn.id] = prot_exchange_rxn
            print(prot_exchange_rxns)
        continue

    if len(complexes) == 1:
        complx = complexes[0]
        prots = [item.strip() for item in complx.split('and')]
        #kcat = kcat_dict[(prots,rxn.id)]
        e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn, kcat_dict, complx,rxn_index=None)
        converted_reaction_list.append(e_rxn)
        for prot_exchange_rxn in prot_exchange_rxns: protein_exchange_rxns[prot_exchange_rzxn.id] = prot_exchange_rxn

eModel = Model.copy()
eModel.add_reactions(converted_reaction_list)
eModel.add_reactions(protein_exchange_rxns.values())
eModel.enzymes = set([exgrxn.split('_')[1] for exgrxn in protein_exchange_rxns.keys()])
print('Number of enzymes:', len(eModel.enzymes))

MWs = {}
gene = [rxn[0][0] for rxn in enzymedata['proteins']]
mw = [mw[0]/1000 for mw in enzymedata['proteinMW']]
df = pd.DataFrame ({'MW': mw} , index=gene)
for i in range (len (df)):
    MWs[df.index[i]] = df.loc[df.index[i] , "MW"]

with eModel:
    ecModel = constrainPool (eModel , MWs , {} , eModel.enzymes , 230)

    s2 = get_r_max (ecModel , model_dl)
    r2 = s2.objective_value
    print ('Case 2: without omics constraints')
    print ('  Status        :' , s2.status)
    print ('  growth_rate   :' , r2)
    # print('  glucose uptake:',s2.fluxes['Exchange_Glucopyranose'])
    # print('  PHA           :',s2.fluxes['PHA_secretion'])
    print ('  Protein pool  :' , ecModel.reactions.get_by_id ('prot_pool_exchange').upper_bound)

write_sbml_model(ecModel,"test_ecmodel.xml")
