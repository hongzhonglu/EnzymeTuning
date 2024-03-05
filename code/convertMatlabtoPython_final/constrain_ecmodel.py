
import numpy as np
import os
import cobra
import pandas as pd
from cobra.io import load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import scipy.io as scio
from cobra import Model,Reaction,Metabolite

def mat_to_xml(species):
    '''
    :param species: the name of the dl_model saved in 'data'
    :return: a dl_model in xml format
    '''
    modelfile = 'data/'+str(species)+'_dl.mat'
    model = cobra.io.load_matlab_model(modelfile)
    cobra.io.write_sbml_model(model,str(species)+'_dl.xml')

def gene_id_to_protein_id(model, data):
    '''
    :param model: a cobra model
    :param data: data = pd.read_excel('data/uniprotkb_Scer.xlsx', sheet_name= "Sheet0") contain all gene and protein id
    :return: a cobra model with protein id
    '''
    # create a dictionary to save Gene Names to Entry (protein id)
    gene_to_entry_mapping = {}

    # find the protein id to each gene names (one to several)
    for index, row in data.iterrows():
        entry = row['Entry']
        gene_names = row['Gene Names'].split()
        for gene in gene_names: gene_to_entry_mapping[gene] = entry

    for gene in model.genes:
        if gene.id in gene_to_entry_mapping:
            gene.id = gene_to_entry_mapping[gene.id]
    return model

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
    kcat = [kcat[0] for kcat in enzymedata['kcat']]
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

def parse_gr_rule(gr):
    '''
    Parse gr rule into a list of components.
    gr: gene reaction rule, defined in cobrapy.

    For example:

    Input         : Output
    A or B        : ["A", "B"]
    (A and B)     : ["A and B"]
    (A and B) or C: ["A and B","C"]

    Usage: complexes = parse_gr_rule(gr)

    Gang Li, last updated 2020-03-04

    '''
    complexes = [item.strip().replace('(','').replace(')','') for item in gr.split('or')]
    if len(complexes) < 2 and len(complexes[0]) < 1: complexes = []

    return complexes

def addEnzymesToRxn(rxn, kcat_dict, protIDs, rxn_index=None):
    '''
    Add each enzyme as one of metabolites in the model, Current version does not support stoichiometric cofficients of subunits
    rxn      : the input Reaction object in cobrapy
    rxn_index: a integer like 1, for isoenzymes
    kcats    : kcat value for the reaction
    protIDs  : a single protein name, like "A", or a complex like "A and B". String
    MWs      : a dictionary with prot_id as key and molecular weight as value

    Usage: e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn, kcat, protIDs,MWs)

    Gang Li, last updated 2020-03-03
    '''

    e_rxn      = rxn.copy()
    if rxn_index is not None:
        e_rxn.id   = e_rxn.id + 'No{0}'.format(rxn_index)
        e_rxn.name = e_rxn.name + ' (No{0})'.format(rxn_index)
    prots = [item.strip() for item in protIDs.split('and')]


    # get compartment
    #comp = None
    #for met in rxn.metabolites:
    #    comp = met.compartment
    #    if rxn.get_coefficient(met)<0: comp = met.compartment

    comp = "c"

    # create Metabolite object for each protein and create exchange reaction
    prot_mets = []
    prot_exchange_rxns = []
    for prot in prots:
        prot_met = Metabolite('prot_{0}[{1}]'.format(prot,comp))
        prot_met.compartment =  comp
        prot_mets.append(prot_met)

        # add excange reaction of protein
        excg_rxn = Reaction('prot_{0}'.format(prot))
        excg_rxn.lower_bound = 0
        excg_rxn.gene_reaction_rule = prot
        excg_rxn.add_metabolites({prot_met:1})
        prot_exchange_rxns.append(excg_rxn)

    # add enzymes into reaction
    # create an empty dict to store the coefficients
    enzyme_coefficients = {}

    for prot in prots:
        kcat_value = kcat_dict.get((prot, rxn.id), None)  # 如果没有找到，则返回 None
        if kcat_value is not None:
            enzyme_coefficients[prot] = -1.0 / kcat_value
    # add enzymes into reaction
    e_rxn.add_metabolites({prot_met: enzyme_coefficients.get(prot_met.id.split('_')[1].split('[')[0], 0) for prot_met in prot_mets})
    e_rxn.gene_reaction_rule = protIDs

    return e_rxn, prot_exchange_rxns

def getArmReaction(rxn):
    '''
    Adapted from addArmReaction.m from geckomat. Add an arm reaction for the selected reaction in the model.

    rxn: the reaction Object in cobrapy

    Original reaction: A + B --> C + D

    Arm reaction    : A + B --> pmet   (no gr rule)
    Change the orginial reaction to:  pmet --> C + D (use old gr rule)

    The arm reaction has a id format of "arm_rxnID" and a name format of "rxnName (arm)"

    The intermediate metabilite has a name format of "pmet_rxnID"

    The arm reaction shares the same lb, ub, gr rules, subsystems with original reaction.

    Compartment: fistly try to use the same compartment as substrates, then products', otherwise None.


    Usage: rxn_new, arm_rxn = addArmReaction(model,rxn_id).

    Gang Li, Last update: 2020-03-03
    '''

    # 1. create intermediate metabilite
    rxnID = rxn.id
    comp = None
    for met in rxn.metabolites:
        comp = met.compartment
        if rxn.get_coefficient(met)<0: comp = met.compartment

    pmet = Metabolite('pmet_{0}'.format(rxnID),compartment=comp)

    # 2. create arm reaction:
    arm_rxn                    = Reaction('arm_{0}'.format(rxnID))
    arm_rxn.name               = rxn.name + ' (arm)'
    arm_rxn.subsystem          = rxn.subsystem
    arm_rxn.lower_bound        = rxn.lower_bound
    arm_rxn.upper_bound        = rxn.upper_bound
    arm_rxn.gene_reaction_rule = ''

    mets = {met:rxn.get_coefficient(met) for met in rxn.metabolites if rxn.get_coefficient(met)<0}
    mets[pmet] = 1

    arm_rxn.add_metabolites(mets)

    # 3. change orignal reaction to pmet --> C + D
    rxn_new = rxn.copy()
    rxn_new.subtract_metabolites({met:rxn_new.get_coefficient(met) for met in rxn_new.metabolites if rxn_new.get_coefficient(met)<0})
    rxn_new.add_metabolites({pmet:-1})

    return rxn_new, arm_rxn

def constrainAbandance(model,measured):
    '''
    model       : eModel from convertToEnzymeModel()
    measured    : a dictionary with measured enzyme abandance, in the unit of mmol/gdw
    g
    # define the upper bound of protein exchange reactions with protein abandance.
    # e.g. the reaction id is in the format of "prot_TD01GL001367_exchange"

    Usage: model = constrainAbandance(model,MWs, non_measured,UB)
    '''

    for prot_id, ab in measured.items():
        rxn_id = 'prot_{0}_exchange'.format(prot_id)
        model.reactions.get_by_id(rxn_id).upper_bound = ab

    return model

def constrainPool(model,MWs, measured, non_measured,UB,copy=True):
    '''

    model       : eModel from convertToEnzymeModel()
    MWs         : a dictionary with molecular weight of enzymes, in the unit of kDa
    non_measured: a list of enzymes without proteomics data
    measured    : a dictionary with measured enzyme abandance, in the unit of mmol/gdw
    UB          : upper bound for the combined pool of those non_measured enzymes
    copy        : if creat a copy of the original model

    Define new rxns: For each enzyme, add a new rxn that draws enzyme from the
    enzyme pool (a new metabolite), and remove previous exchange rxn. The new
    rxns have the following stoichiometry (T is the enzyme pool):
     MW[i]*P[T] -> P[i]

    Usage: model = constrainPool(model,MWs, non_measured,UB)

    Gang Li, last updated 2020-03-04
    '''
    if copy: model = model.copy()
    # create prot_pool metabolite
    prot_pool = Metabolite('prot_pool')
    prot_pool.name = prot_pool.id
    prot_pool.compartment = 'c'

    rxns_to_add  = list()
    rxns_to_drop = list()
    for prot in non_measured:
        #prot_exchange_rxn = model.reactions.get_by_id('prot_{0}_exchange'.format(prot))
        prot_exchange_rxn = model.reactions.get_by_id('prot_{0}'.format(prot))

        draw_rxn = Reaction('draw_prot_{0}'.format(prot))
        draw_rxn.name = draw_rxn.id
        draw_rxn.add_metabolites({prot_pool:-MWs[prot],list(prot_exchange_rxn.metabolites)[0]:1})

        rxns_to_add.append(draw_rxn)
        rxns_to_drop.append(prot_exchange_rxn)

    # add draw reaction into model
    model.add_reactions(rxns_to_add)
    model.remove_reactions(rxns_to_drop)

    # change the upper bound for all reactions as np.inf
    for rxn in model.reactions: rxn.upper_bound = np.inf

    # add prot_pool_exchange rxn
    rxn_prot_pool_exg = Reaction('prot_pool_exchange')
    rxn_prot_pool_exg.name = rxn_prot_pool_exg.id
    rxn_prot_pool_exg.add_metabolites({prot_pool:1})
    rxn_prot_pool_exg.lower_bound = 0
    rxn_prot_pool_exg.upper_bound = UB

    model.add_reactions([rxn_prot_pool_exg])

    # constrain the proteins with measure abandance
    constrainAbandance(model,measured)

    return model

def get_r_max(model,model_dl):
    model.objective = model_dl.objective
    model.objective_direction = 'max'
    return model.optimize()

def convertToEnzymeModel(model,kcats):
    '''
    model .   : irrevModel
    kcats     : a dictionary with kcat values {('protein_id',rxn_id):100,...}

    Usage: eModel = convertToEnzymeModel(model,kcats)

    Gang Li, last updated 2020-03-04
    '''
    converted_reaction_list = []
    protein_exchange_rxns = {}
    for rxn in model.reactions:
        complexes = parse_gr_rule(rxn.gene_reaction_rule)

        # 1. for those reactions without genes
        if len(complexes) <1:
            converted_reaction_list.append(rxn)
            continue

        # 2. for those reactions with genes, but no kcat
        first_gene = [gene.id for gene in rxn.genes][0]
        if kcats.get((first_gene,rxn.id),None) is None:
            converted_reaction_list.append(rxn)
            continue

        # 3. for those reactions with isoenzymes, add arm reaction
        if len(complexes) >1:
            rxn_new, arm_rxn = getArmReaction(rxn)
            converted_reaction_list.append(arm_rxn)

            for i,complx in enumerate(complexes):
                prots = [item.strip() for item in complx.split('and')]
                kcat = kcats[(prots[0],rxn.id)]
                e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn_new, kcat, complx,rxn_index=i+1)

                converted_reaction_list.append(e_rxn)
                for prot_exchange_rxn in prot_exchange_rxns: protein_exchange_rxns[prot_exchange_rxn.id] = prot_exchange_rxn

            continue

        if len(complexes) == 1:
            complx = complexes[0]
            prots = [item.strip() for item in complx.split('and')]
            kcat = kcats[(prots[0],rxn.id)]
            e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn, kcat, complx,rxn_index=1)
            converted_reaction_list.append(e_rxn)
            for prot_exchange_rxn in prot_exchange_rxns: protein_exchange_rxns[prot_exchange_rxn.id] = prot_exchange_rxn

    eModel = Model()
    eModel.add_reactions(converted_reaction_list)
    eModel.add_reactions(protein_exchange_rxns.values())
    eModel.enzymes = set([exgrxn.split('_')[1] for exgrxn in protein_exchange_rxns.keys()])
    print('Number of enzymes:', len(eModel.enzymes))

    return eModel


