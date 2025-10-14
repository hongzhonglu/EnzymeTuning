from cobra.io import read_sbml_model
import pandas as pd
from src.mainFunction import *



def getALLGEMgene():
    """
    The function is used to get the metabolic gene list.
    :return:
    """
    GEM_yeast = read_sbml_model('/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xml')
    gene_yeast = []
    for x in GEM_yeast.genes:
        print(x.id)
        gene_yeast.append(x.id)
    return gene_yeast



def getTransporterCellMembraneGEM():
    """
    Get the transporter proteins of plasma membrane from GEMs
    :return:
    """
    GEM_subsystem = pd.read_csv('data/subsystem/Rxn_unique_subsystem_v2.tsv', sep="\t")
    GEM_subsystem_t = GEM_subsystem[GEM_subsystem['subsystem_unique'].str.contains("Transport")]
    GEM_subsystem_t = GEM_subsystem_t[~GEM_subsystem_t['subsystem_unique'].str.contains("er")]
    GEM_subsystem_t = GEM_subsystem_t[
        GEM_subsystem_t['subsystem_unique'].str.contains("e") | GEM_subsystem_t['subsystem_unique'].str.contains("ce")]
    GEM_subsystem_t = GEM_subsystem_t[~GEM_subsystem_t["GENE ASSOCIATION"].isna()]
    # get the related gene list
    # there 121 genes from GEMs
    rxn_gene = getRXNgeneMapping(rxn0=GEM_subsystem_t['ID'], gpr0=GEM_subsystem_t['GENE ASSOCIATION'])
    m_gene_plasma_membrane = list(set((rxn_gene['gene'])))
    # analyze glucose
    #glucose_transporter = rxn_gene[rxn_gene['rxnID'] == rxnID]
    #m_glucose_transporter = list(set((glucose_transporter['gene'])))
    return m_gene_plasma_membrane



def getProteinForRxnGEM(rxnID):
    """
    Get the proteins based on reaction IDs
    :param rxnID: a list with gene id, like ['r_1166']
    :return:
    """
    GEM_subsystem = pd.read_csv('data/subsystem/Rxn_unique_subsystem_v2.tsv', sep="\t")
    rxn_gene = getRXNgeneMapping(rxn0=GEM_subsystem['ID'], gpr0=GEM_subsystem['GENE ASSOCIATION'])
    # analyze glucose
    df = rxn_gene[rxn_gene['rxnID'].isin(rxnID)]
    pro_list = list(set((df['gene'])))
    return pro_list



def exchange_ecYeast(s1, subystem):
    """
    this function is used to define the exchange reaction
    s1=['a --> b','a <=> c', 'H+ [extracellular] + L-citrulline [extracellular] <=> H+ [cytoplasm] L-citrulline [cytoplasm]', ' a--> ']
    subsystem = ['a','a','b','']

    """
    for i, x in enumerate(s1):
        print(i)
        if ' --> ' in x:
            x0 = x.split(' --> ')
            if len(x0[1]) >=1 and len(x0[0]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        if ' <=> ' in x:
            x0 = x.split(' <=> ')
            if len(x0[1]) >=1 and len(x0[0]) >=1:
                #subystem.append('General')  # exchange
                subystem[i] = subystem[i]
            else:
                subystem[i] ='Exchange reaction' #exchange
                print(subystem[i])
        else:
            subystem[i] = subystem[i]
    return subystem



# function copy from strain_design repo
def ecYeastMinimalMedia(model):
    """
    This function is used to define a simple media for ecYeast
    :param model:
    :return: a model with the defined the minimal media
    """
    rxnID = []
    rxnName = []
    for i, x in enumerate(model.reactions):
        rxnID.append(x.id)
        rxnName.append(x.name)

    exchange_rxn =[x for x, y in zip(rxnID, rxnName) if '_REV' in x and 'exchange' in y]
    # first block any uptake
    for i, x in enumerate(exchange_rxn):
        rxn0 = exchange_rxn[i]
        #print(rxn0)
        model.reactions.get_by_id(rxn0).upper_bound = 0

    #Allow uptake of essential components
    model.reactions.get_by_id("r_1654_REV").upper_bound = 10000 #ammonium exchange (reversible)
    model.reactions.get_by_id("r_1861_REV").upper_bound = 10000 #iron(2+) exchange (reversible)
    model.reactions.get_by_id("r_2100_REV").upper_bound = 10000 #water exchange (reversible)
    model.reactions.get_by_id("r_1992_REV").upper_bound = 10000 #oxygen exchange (reversible)
    model.reactions.get_by_id("r_2005_REV").upper_bound = 10000 #phosphate exchange (reversible)
    model.reactions.get_by_id("r_2060_REV").upper_bound = 10000 #sulphate exchange (reversible)
    model.reactions.get_by_id("r_1832_REV").upper_bound = 10000 #H+ exchange (reversible)
    return model



def chemostatSimulation(model0, D0):
    """
    This funcion is used to simulate the chemostat growth of yeast
    Actually this function is general to solve the ecGEMs
    :param model0: a ecGEMs
    :param D0: a growth rate
    :return: solution of fluxes
    """
    growth = D0
    with model0:
        model0 = ecYeastMinimalMedia(model0)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth, growth)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution2 = model0.optimize()
        GR = solution2.fluxes["r_1714_REV"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution3 = model0.optimize()
    return solution3



def simulationCompare(model_in, growth_in, objective):
    """
    This function is used to compare the result when different objective function is employed.
    :param model_in:
    :param growth_in:
    :param objective: a kind of objective for ecModel_batch, which could minimize the total protein volume
    :return:
    """

    # 1_minimize the glucose uptake and minimize protein pool
    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution00 = model0.optimize()
        GR = solution00.fluxes["r_1714_REV"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
        # then minimization protein usage
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution1 = model0.optimize()

    # 2_minimize the glucose uptake
    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 10)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution2 = model0.optimize()
    # 3_minimize protein pool
    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # only minimization protein usage
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 10)  # open the glucose
        model0.objective = {model0.reactions.prot_pool_exchange: -1}
        solution3 = model0.optimize()
    # 4_minimize protein volume
    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 10)  # open the glucose
        model0.objective = objective
        solution4 = model0.optimize()

    # 5_minimize the glucose uptake and protein volume
    with model_in:
        model0 = ecYeastMinimalMedia(model_in)
        # set growth
        model0.reactions.get_by_id("r_2111").bounds = (growth_in, growth_in)
        # minimization glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose
        model0.objective = {model0.reactions.r_1714_REV: -1}
        solution00 = model0.optimize()
        GR = solution00.fluxes["r_1714_REV"]  # get the glucose uptake rate
        model0.reactions.get_by_id("r_1714_REV").bounds = (GR, GR * 1.001)
        model0.objective = objective
        solution5 = model0.optimize()
    return solution1, solution2, solution3, solution4, solution5




def getRxnByGene(model, gene0):
    """
    :param model: A metabolic model
    :param gene0: A gene
    :return:

    Example:
    getRxnByGene(model, gene0="YGR192C")

    """
    rxn_list2 = []
    for gene in model.genes:
        if gene.id == gene0:
            rxn_list = gene.metabolites
            for x in rxn_list:
                rxn_list2.append(x.id)
    return rxn_list2



"""
def getRxnByReactionName(model, name):
    for rxn in model.reactions:
        if name in rxn.name:
            print(rxn.id)
            return rxn.id
        else:
            return None
"""



def getRxnByReactionName(model, name):
    """
    This function is used to extract the rxn id based on rxn name
    It is suitable for ecGEMs as multiple reations could use the same name

    :param model:
    :param name:
    :return:
    """
    s = []
    for rxn in model.reactions:
        if name == rxn.name:
            #print(rxn.id)
            s.append(rxn.id)
    return s


def DLecModelSimulate(model, dilution_rate):
   """
   This function is used to do simulation with ecModels using kcat value from deep learning.
   :param model: a ecModel
   :param dilution_rate: a dilution rate in range of 0-0.42 /h

   :return: solution_f: fluxes datasets

   """

   dilutionrate = dilution_rate
   ecYeast = model
   if dilutionrate >= 0.4:
       ecYeast.reactions.get_by_id("prot_pool_exchange").bounds = (-167.27 * dilutionrate / 0.4, 0)  # this value is further rescaled by maximal growth rate at 0.42.
   else:
       ecYeast.reactions.get_by_id("prot_pool_exchange").bounds = (-167.27, 0)  # -230/0.55*0.4, this is rescaled by maximal growth rate.

   # refer to bioRxiv
   ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
              'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
              'oxygen exchange', 'prot_pool_exchange']
   # find the related rxnID
   idx = []
   for name0 in ex_mets:
       #print(name0)
       s = getRxnByReactionName(model=ecYeast, name=name0)
       if len(s) > 1:
           print("need check")
       elif len(s) == 1:
           idx.append(s[0])

   model_tmp = ecYeast #.copy() # this is adjusted for the newly added constraints from organelle
   model_tmp.reactions.get_by_id(idx[1]).lower_bound = -1000  # glucose uptake
   model_tmp.reactions.get_by_id(idx[0]).lower_bound = dilutionrate
   model_tmp.objective = {model_tmp.reactions.r_1714: 1}  # minimize the uptake of glucose
   solution2 = model_tmp.optimize()
   # then fix glucose uptake and minimize the protein pool
   model_tmp.reactions.get_by_id(idx[1]).lower_bound = solution2.objective_value * 1.00001
   print('Glucose uptake rate: ',solution2.objective_value)
   model_tmp.reactions.get_by_id(idx[9]).lower_bound = -1000 # protein pool
   model_tmp.objective = {model_tmp.reactions.prot_pool_exchange: 1}  # minimize the usage of protein pools
   solution_f = model_tmp.optimize()
   solution_f.fluxes["prot_pool_exchange"]
   solution_f.fluxes["r_1714"]
   return solution_f


def updateEcGEMkcat(ecGEM, target_gene, rxnID, kcat_m):
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
    ecModel = ecGEM.copy()
    coef = 1 / kcat_m
    ss = ecModel.reactions.get_by_id(rxnID).reaction
    updateEcGEMkcat

