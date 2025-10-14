import scipy.io as scio
from cobra.io import load_matlab_model , save_matlab_model , read_sbml_model , write_sbml_model
import numpy as np
import pandas as pd
import time
import os
import pickle

def construct_MWs(enzymedata):
    MWs = {}
    gene = [rxn[0][0] for rxn in enzymedata['proteins']]
    mw = [mw[0] for mw in enzymedata['proteinMW']]
    df = pd.DataFrame ({'MW': mw} , index=gene)
    for i in range (len (df)):
        MWs[df.index[i]] = df.loc[df.index[i] , "MW"]
    return MWs

def construct_kcat_dict(model,enzymedata):
    '''
    :param model: cobra.Model object (with irreversible reactions)
    :param enzymedata: a mat file contains rxn_list, kcat and kcat_std values
    :return: a dictionary with (enz_id,rxn_id) as keys and kcats as values kcat = kcat/subunitcoef
    E.g.
    kcat_dict = {
                     ('E1','R3'): 10,
                     ('E1','R4'): 100,
                     ('E2','R4'): 0.10,
                     ('E4','R5_REV'): 90,
                  }
    '''
    idx = []
    idx2 = []
    idx3 = []
    kcat_dict = {}
    rxn_list = [rxn[0][0] for rxn in enzymedata['rxn_list']]
    #kcat = [kcat[0] for kcat in enzymedata['kcat']]
    kcat = enzymedata['kcat']
    for i in range(len(enzymedata['rxn_list'])):
        try:
            idx_tmp = np.where(model['rxns'] == enzymedata['rxn_list'][i])
            idx_tmp = idx_tmp[0][0]
        except ValueError:
            idx_tmp = None
        idx.append(idx_tmp)
        subunitlist = enzymedata['subunit'][i,:]
        subunit_num = [bool(x) for x in subunitlist]
        subunitlist = [subunit for subunit, is_non_empty in zip(subunitlist, subunit_num) if is_non_empty]
        subunitlist = np.array(subunitlist)
        subunitlist_tmp = [f"{subunit[0]}" for subunit in subunitlist]
        idx2.append(subunitlist_tmp)
        subunitcoef = enzymedata['subunit_stoichiometry'][i,subunit_num]
        idx3.append(subunitcoef)
    for rxn, subunits, kcat_values,c in zip(rxn_list, idx2, kcat,idx3):
        for subunit, coef_value in zip(subunits, c):
            if subunit:
                kcat_dict[(subunit, rxn)] = kcat_values/coef_value
    return kcat_dict

def updateEcGEMkcat(ecModel, target_gene, rxnID, kcat_m):
    """
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


enzymedataFile = "data/Saccharomyces_cerevisiae_dl.mat" #enzymedata from mat
enzyme_preparation = pd.read_csv ('data/kcat_merge_error_50.csv') #new kcat list
z = scio.loadmat (enzymedataFile)
model = z['model'][0 , 0]
enzymedata = z['enzymedata'][0 , 0]
MWs = construct_MWs (enzymedata)
kcat_dict = construct_kcat_dict (model , enzymedata)
kcat_df = pd.DataFrame (list (kcat_dict.items ()) , columns=['gene_rxnID' , 'kcat_value'])
kcat_df[['gene' , 'rxnID']] = pd.DataFrame (kcat_df['gene_rxnID'].tolist () , index=kcat_df.index)
merged_df = pd.merge (enzyme_preparation , kcat_df , on=['gene' , 'rxnID'])
kcat_final = merged_df[['gene' , 'rxnID' , 'kcat_value']]
eModel = load_matlab_model ('data/model/model.mat') #ecModel

with eModel as model_tmp :
    for k in range (len (kcat_final)) :
        target_gene0 = kcat_final[ 'gene' ][ k ]
        kcat_m0 = kcat_final[ 'kcat_value' ][ k ]
        rxn0 = kcat_final[ 'rxnID' ][ k ]
        model_tmp = updateEcGEMkcat (ecModel=model_tmp,target_gene=target_gene0,rxnID=rxn0,kcat_m=kcat_m0)