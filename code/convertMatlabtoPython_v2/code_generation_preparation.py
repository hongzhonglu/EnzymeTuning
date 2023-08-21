import numpy as np
import os
import cobra
from pathlib import Path
from sklearn.metrics import mean_squared_error
import pandas as pd
from cobra.io import load_matlab_model , save_matlab_model , read_sbml_model , write_sbml_model
import scipy.io as scio
from scipy import sparse
import cobra.util
from cobra import Model , Reaction , Metabolite
import cProfile
import time

path = "F:\python\Bayesian_python"
os.chdir (path)

def getFraction(model , data , compType , X):
    '''
    X: Biomass fraction without lipids [g/gDW]
    P: Protein fraction [g/gDW]
    C: Carbohydrate fraction [g/gDW]
    R: RNA fraction [g/gDW]
    D: DNA fraction [g/gDW]
    L: Lipid fraction [g/gDW]
    model: dl_model.mat
    data: biomassCompData.mat
    '''
    # define pseudoreaction name
    rxnName = compType + ' pseudoreaction'
    rxnName = rxnName.replace ('P' , 'protein')
    rxnName = rxnName.replace ('C' , 'carbohydrate')
    rxnName = rxnName.replace ('N' , 'biomass')
    rxnName = rxnName.replace ('L' , 'lipid backbone')
    rxnName = rxnName.replace ('R' , 'RNA')
    rxnName = rxnName.replace ('D' , 'DNA')
    rxnName = rxnName.replace ('I' , 'ion')
    rxnName = rxnName.replace ('F' , 'cofactor')

    # add up fraction
    rxnPos = []
    for i in range (len (model['rxnNames'])):
        if rxnName==model['rxnNames'][i]:
            temp = 1
        else:
            temp = 0
        rxnPos.append (temp)
    index = [i for i , e in enumerate (rxnPos) if e!=0]
    if np.nonzero (rxnPos)!=0:
        sub = model['S']
        isSub = sub.getcol (index[0]) < 0  # substrates in pseudo-rxn
        if compType=='L':
            F = - sum (sub[isSub , rxnPoS])  # g/gDW
        else:
            F = 0
            # add up all components:
            for i in range (len (model['mets'])):
                mets = model['mets']
                pos = []
                if mets[i] in data['mets']:
                    temp = np.where (data['mets']==mets[i])
                    pos.append (temp)
                if isSub[i] and len (pos) > 0:
                    if compType=='I' or compType=='F':
                        MW = data['MWs'][pos[0][0]]
                    else:
                        MW = data['MWs'][pos[0][0]] - 18
                    abundance = -sub[i , index] * MW / 1000
                    F = F + abundance
        X = X + F

        print (str (compType) + ' -> ' + str (F) + " g/gDW")
    else:
        print (str (compType) + " do not exist")
        F = 0
        X = X + F
    print ("X ->" + str (X) + " gDW/DW")

    return X , F


def sumBiomass(model , modelcobra):
    '''
    Calculates breakdown of biomass for the yeast model:
    get main fraction (protein) in S.cer 0.46
    '''

    import scipy.io as scio
    biomassCompDataFile = "data/biomassCompData.mat"
    biomassCompData = scio.loadmat (biomassCompDataFile)
    biomassCompData = biomassCompData['data'][0 , 0]

    #Get main fractions:
    P = getFraction(model, biomassCompData, 'P', 0)
    # C = getFraction(model,biomassCompData,'C',X)
    # R = getFraction(model,biomassCompData,'R',X)
    # D = getFraction(model,biomassCompData,'D',X)
    # L = getFraction(model,biomassCompData,'L',X)
    # I = getFraction(model,biomassCompData,'I',X)
    # F = getFraction(model,biomassCompData,'F',X)

    # print("X ->"+str(X)+" gDW/DW")
    sol = modelcobra.optimize ()
    print ("Growth = " + str (sol.objective_value) + " 1/h")

    return P


def getExchangeRxns(model , reactionType="both"):
    '''
    Retrieves the exchange reactions from a model
    Exchange reactions are defined as reactions which involve only products
    or only reactants. If the unconstrained field is present, then that is
    used instead.

    model: a model structure
    reactionType: retrieve all reactions ('both'), only production ('out'),
    or only consumption ('in') (opt, default 'both')

    exchangeRxns: cell array with the IDs of the exchange reactions
    exchangeRxnsIndexes: vector with the indexes of the exchange reactions
    '''
    hasNoProducts = []
    hasNoReactants = []
    for i in range (len (model['rxns'])):
        isSub = model['S'].getcol (i)
        if reactionType=="both" or reactionType=="out":
            if not (isSub > 0).getnnz ():
                hasNoProducts.append (i)
        if reactionType=="both" or reactionType=="in":
            if not (isSub < 0).getnnz ():
                hasNoReactants.append (i)
    exchangeRxnIndexes = hasNoProducts + hasNoReactants
    exchangeRxns = model["rxns"][exchangeRxnIndexes]
    return exchangeRxnIndexes , exchangeRxns


def changeMedia(model , model_cobra , c_source , media , anox=False , flux=-1000):
    '''
    Function that modifies the ecModel and makes it suitable for batch growth
    simulations on different carbon sources.

    media:  Media type ('YEP' for complex,
                      'MAA' minimal with Aminoacids,
                      'MIN' for minimal media,
                      'MIN+His' for minimal media with his
                      'MIN+Arg' for minimal media with arg
                      'MIN+Citrate' for minimal media with Citrate)
    anox:   (optional) TRUE if anaerobic conditions are desired, DEFAULT= FALSE
    flux:   (Optional) A cell array with measured uptake fluxes in mmol/gDwh

    '''
    #Give the carbon source (c_source) input variable with the following format:
    #c_source = 'D-glucose exchange'
    c_source = c_source + ' exchange'

    # first block any uptake
    exchangerxn = getExchangeRxns (model , reactionType="both")
    for a in range (len (exchangerxn[0])):
        model_cobra.reactions[exchangerxn[0][a]].lower_bound = 0
    pos = getComponentIndexes (model , c_source)

    # The media will define which rxns to fix:
    if media=='YEP':
        N = 25  # Aminoacids + Nucleotides
    elif media=='MAA':
        N = 21  # Aminoacids
    elif media=='MIN':
        N = 1  # Only the carbon source
    elif media=='MIN+His':
        N = 1  # Only the carbon source
        model_cobra.reactions.get_by_id ('r_1893').lower_bound = -0.08  # Histidine exchange
    elif media=='MIN+Arg':
        N = 1  # Only the carbon source
        model_cobra.reactions.get_by_id ('r_1879').lower_bound = -0.08  # L-arginine exchange
    elif media=='MIN+Citrate':
        N = 1  # Only the carbon source
        model_cobra.reactions.get_by_id ('r_1687').lower_bound = -0.08  # citrate exchange
    # LB parameter (manually optimized for glucose on Min+AA):
    b = -0.08
    # LB parameter (manually optimized for glucose complex media):
    c = -2
    flux = np.array (flux)
    # Define fluxes in case of ec model:
    if N > 1:
        flux = b * np.ones ((1 , N))
        if N > 21:
            flux[21:25] = c

    # Fix values as LBs:
    i = 0
    for i in range (N):
        model_cobra.reactions[pos[i][0]].lower_bound = -1000  # flux[i]

    # Allow uptake of essential components
    model_cobra.reactions.get_by_id ('r_1654').lower_bound = -1000  # 'ammonium exchange';
    model_cobra.reactions.get_by_id ('r_2100').lower_bound = -1000  # 'water exchange';
    model_cobra.reactions.get_by_id ('r_1861').lower_bound = -1000  # 'iron(2+) exchange';
    model_cobra.reactions.get_by_id ('r_1992').lower_bound = -1000  # 'oxygen exchange';
    model_cobra.reactions.get_by_id ('r_2005').lower_bound = -1000  # 'phosphate exchange';
    model_cobra.reactions.get_by_id ('r_2060').lower_bound = -1000  # 'sulphate exchange';
    model_cobra.reactions.get_by_id ('r_1832').lower_bound = -1000  # 'H+ exchange';
    model_cobra.reactions.get_by_id ('r_4593').lower_bound = -1000  # 'chloride exchange';
    model_cobra.reactions.get_by_id ('r_4595').lower_bound = -1000  # 'Mn(2+) exchange';
    model_cobra.reactions.get_by_id ('r_4596').lower_bound = -1000  # 'Zn(2+) exchange';
    model_cobra.reactions.get_by_id ('r_4597').lower_bound = -1000  # 'Mg(2+) exchange';
    model_cobra.reactions.get_by_id ('r_2049').lower_bound = -1000  # 'sodium exchange';
    model_cobra.reactions.get_by_id ('r_4594').lower_bound = -1000  # 'Cu(2+) exchange';
    model_cobra.reactions.get_by_id ('r_4600').lower_bound = -1000  # 'Ca(2+) exchange';
    model_cobra.reactions.get_by_id ('r_2020').lower_bound = -1000  # 'potassium exchange';

    # Block some production fluxes
    model_cobra.reactions.get_by_id ('r_1663').upper_bound = 0  # bicarbonate exchange;
    model_cobra.reactions.get_by_id ('r_4062').upper_bound = 0  # lipid backbone exchange;
    model_cobra.reactions.get_by_id ('r_4064').upper_bound = 0  # lipid chain exchange;

    # Allow biomass production
    model_cobra.reactions.get_by_id ('r_2111').upper_bound = 1000  # growth;

    return model_cobra

def getComponentIndexes(model , c_source):
    # c_source = c_source + ' exchange'
    pos = []
    pos.append (np.where (model["rxnNames"]==c_source))
    pos.append (np.where (model["rxnNames"]=='L-alanine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-arginine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-asparagine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-aspartate exchange'))
    pos.append (np.where (model["rxnNames"]=='L-cysteine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-glutamine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-glutamate exchange'))
    pos.append (np.where (model["rxnNames"]=='L-glycine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-histidine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-isoleucine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-leucine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-lysine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-methionine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-phenylalanine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-proline exchange'))
    pos.append (np.where (model["rxnNames"]=='L-serine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-threonine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-tryptophan exchange'))
    pos.append (np.where (model["rxnNames"]=='L-tyrosine exchange'))
    pos.append (np.where (model["rxnNames"]=='L-valine exchange'))
    pos.append (np.where (model["rxnNames"]=="2'-deoxyadenosine exchange"))
    pos.append (np.where (model["rxnNames"]=="2'-deoxyguanosine exchange"))
    pos.append (np.where (model["rxnNames"]=='thymidine exchange'))
    pos.append (np.where (model["rxnNames"]=='deoxycytidine exchange'))
    pos.append (np.where (model["rxnNames"]=='D-glucose exchange'))
    pos = np.array (pos)
    pos = pos[: , 0]
    return pos


def anaerobicModel(model , model_cobra):
    #Convert model to anaerobic
    # 1th change: Refit GAM and NGAM to exp. data, change biomass composition
    GAM = 58.1988  # Data from Nissen et al. 1997
    strain = str (model['id']).replace (' specific model genereted from panYeast' , '')
    if strain=='Candida_glabrata' or strain=='Candida_parapsilosis':
        GAM = 30

    #P = 0.461 #Data from Nissen et al. 1997
    NGAM = 0  # Refit done in Jouthen et al. 2012
    model_cobra = changeGAM (model , model_cobra , GAM , NGAM)

    # 2nd change: Removes the requirement of heme a in the biomass equation
    mets = ['s_3714[c]' , 's_1198[c]' , 's_1203[c]' , 's_1207[c]' , 's_1212[c]' , 's_0529[c]']
    for i in range (len (mets)):
        model_cobra.reactions.get_by_id ('r_4598').add_metabolites (
            {mets[i]: 0 - model_cobra.reactions.get_by_id ('r_4598').get_coefficient (mets[i])})

    # 3st change: Changes media to anaerobic
    # no O2 uptake and allows sterol and fatty acid exchanges
    model_cobra.reactions.get_by_id ('r_1992').lower_bound = 0      #O2
    model_cobra.reactions.get_by_id ('r_1757').lower_bound = -1000  #ergosterol
    model_cobra.reactions.get_by_id ('r_1994').lower_bound = -1000  #lanosterol
    model_cobra.reactions.get_by_id ('r_1915').lower_bound = -1000  #palmitoleate
    model_cobra.reactions.get_by_id ('r_2106').lower_bound = -1000  #zymosterol
    model_cobra.reactions.get_by_id ('r_2134').lower_bound = -1000  #14-demethyllanosterol
    model_cobra.reactions.get_by_id ('r_2137').lower_bound = -1000  #ergosta-5,7,22,24(28)-tetraen-3beta-ol
    model_cobra.reactions.get_by_id ('r_2189').lower_bound = -1000  #oleate

    # extra media set up for ura1 original but also the anaerobic growth media
    if len (model["id"]):
        strain = str (model['id']).replace (' specific model genereted from panYeast' , '')
    species_onlyura9 = 'Alloascoidea_hylecoeti;Ambrosiozyma_kashinagacola;Ambrosiozyma_monospora;Arxula_adeninivorans;Ascoidea_asiatica;Ascoidea_rubescens;Ashbya_aceri;Aspergillus_nidulans;Babjeviella_inositovora;Brettanomyces_anomalus;Candida_albicans;Candida_apicola;Candida_arabinofermentans;Candida_auris;Candida_boidinii_JCM9604;Candida_carpophila;Candida_dubliniensis;Candida_glabrata;Candida_homilentoma;Candida_infanticola;Candida_intermedia;Candida_orthopsilosis;Candida_parapsilosis;Candida_sorboxylosa;Candida_succiphila;Candida_tanzawaensis;Candida_tenuis;Candida_tropicalis;Candida_versatilis;Clavispora_lusitaniae;Cyberlindnera_fabianii_JCM3601;Cyberlindnera_jadinii;Debaryomyces_hansenii;Dekkera_bruxellensis;Eremothecium_coryli;Eremothecium_cymbalariae;Eremothecium_gossypii;Eremothecium_sinecaudum;Geotrichum_candidum;Hanseniaspora_uvarum;Hanseniaspora_valbyensis;Hanseniaspora_vinae;Hyphopichia_burtonii;Komagataella_pastoris;Kuraishia_capsulata;Lipomyces_starkeyi;Lodderomyces_elongisporus;Metschnikowia_aberdeeniae;Metschnikowia_arizonensis;Metschnikowia_bicuspidata;Metschnikowia_borealis;Metschnikowia_bowlesiae;Metschnikowia_cerradonensis;Metschnikowia_continentalis;Metschnikowia_dekortum;Metschnikowia_drakensbergensis;Metschnikowia_hamakuensis;Metschnikowia_hawaiiensis;Metschnikowia_hibisci;Metschnikowia_ipomoeae;Metschnikowia_kamakouana;Metschnikowia_kipukae;Metschnikowia_lockheadii;Metschnikowia_matae;Metschnikowia_matae_maris;Metschnikowia_mauinuiana;Metschnikowia_proteae;Metschnikowia_santaceciliae;Metschnikowia_shivogae;Metschnikowia_similis;Meyerozyma_guilliermondii;Millerozyma_acaciae;Nadsonia_fulvescens_var_elongata;Nakaseomyces_bracarensis;Nakaseomyces_castellii;Nakaseomyces_delphensis;Nakaseomyces_nivariensis;Nakazawaea_peltata;Ogataea_methanolica;Ogataea_parapolymorpha;Ogataea_polymorpha;Pachysolen_tannophilus;Pichia_membranifaciens;Priceomyces_haplophilus;Saccharomycopsis_malanga;Saprochaete_clavata;Scheffersomyces_lignosus;Scheffersomyces_stipitis;Schizosaccharomyces_pombe;Spathaspora_arborariae;Spathaspora_girioi;Spathaspora_gorwiae;Spathaspora_hagerdaliae;Spathaspora_passalidarum;Sporopachydermia_quercuum;Starmerella_bombicola_JCM9596;Sugiyamaella_lignohabitans;Tortispora_caseinolytica;Vanderwaltozyma_polyspora;Wickerhamia_fluorescens;Wickerhamiella_domercqiae;Wickerhamomyces_anomalus;Wickerhamomyces_ciferrii;Yarrowia_deformans;Yarrowia_keelungensis;Yarrowia_lipolytica;yHMPu5000026124_Ogataea_henricii;yHMPu5000026137_Ambrosiozyma_ambrosiae;yHMPu5000026142_Citeromyces_matritensis;yHMPu5000026145_Ambrosiozyma_vanderkliftii;yHMPu5000026197_Brettanomyces_custersianus;yHMPu5000026274_Komagataella_populi;yHMPu5000034594_Starmera_quercuum;yHMPu5000034597_Candida_stellimalicola;yHMPu5000034604_Sporopachydermia_lactativora;yHMPu5000034605_Spencermartinsiella_europaea;yHMPu5000034606_Priceomyces_medius;yHMPu5000034607_Saccharomycopsis_capsularis;yHMPu5000034610_Saturnispora_hagleri;yHMPu5000034611_Saturnispora_mendoncae;yHMPu5000034612_Saturnispora_saitoi;yHMPu5000034613_Saturnispora_serradocipensis;yHMPu5000034614_Saturnispora_silvae;yHMPu5000034615_Saturnispora_zaruensis;yHMPu5000034622_Pichia_occidentalis;yHMPu5000034623_Pichia_norvegensis;yHMPu5000034624_Pichia_nakasei;yHMPu5000034625_Pichia_kudriavzevii;yHMPu5000034627_Pichia_heedii;yHMPu5000034629_Pichia_exigua;yHMPu5000034631_Martiniozyma_abiesophila;yHMPu5000034632_Candida_athensensis;yHMPu5000034635_Nadsonia_fulvescens;yHMPu5000034636_Ogataea_nitratoaversa;yHMPu5000034637_Ogataea_populiabae;yHMPu5000034643_Candida_schatavii;yHMPu5000034646_Wickerhamiella_cacticola;yHMPu5000034648_Candida_restingae;yHMPu5000034654_Aciculoconidium_aculeatum;yHMPu5000034655_Botryozyma_nematodophila;yHMPu5000034660_Diddensiella_caesifluorescens;yHMPu5000034661_Dipodascus_albidus;yHMPu5000034665_Kodamaea_laetipori;yHMPu5000034667_Blastobotrys_serpentis;yHMPu5000034669_Blastobotrys_raffinofermentans;yHMPu5000034670_Blastobotrys_proliferans;yHMPu5000034671_Blastobotrys_peoriensis;yHMPu5000034673_Blastobotrys_nivea;yHMPu5000034674_Blastobotrys_muscicola;yHMPu5000034675_Blastobotrys_mokoenaii;yHMPu5000034681_Blastobotrys_americana;yHMPu5000034742_Lipomyces_suomiensis;yHMPu5000034748_Lipomyces_oligophaga;yHMPu5000034749_Lipomyces_mesembrius;yHMPu5000034754_Lipomyces_arxii;yHMPu5000034760_Lipomyces_kononenkoae;yHMPu5000034761_Lipomyces_lipofer;yHMPu5000034883_Peterozyma_xylosa;yHMPu5000034884_Peterozyma_toletana;yHMPu5000034885_Ogataea_zsoltii;yHMPu5000034886_Ogataea_trehalophila;yHMPu5000034887_Ogataea_trehaloabstinens;yHMPu5000034890_Ogataea_ramenticola;yHMPu5000034891_Ogataea_pini;yHMPu5000034892_Ogataea_pilisensis;yHMPu5000034893_Ogataea_philodendra;yHMPu5000034897_Ogataea_glucozyma;yHMPu5000034899_Ogataea_kodamae;yHMPu5000034901_Ogataea_methylivora;yHMPu5000034902_Ogataea_minuta;yHMPu5000034903_Ogataea_naganishii;yHMPu5000034904_Ogataea_nonfermentans;yHMPu5000034918_Nakazawaea_holstii;yHMPu5000034933_Kuraishia_molischiana;yHMPu5000034939_Komagataella_pseudopastoris;yHMPu5000034946_Ambrosiozyma_oregonensis;yHMPu5000034947_Ambrosiozyma_philentoma;yHMPu5000034950_Citeromyces_hawaiiensis;yHMPu5000034952_Citeromyces_siamensis;yHMPu5000034957_Hanseniaspora_osmophila;yHMPu5000034963_Hanseniaspora_clermontiae;yHMPu5000034967_Candida_freyschussii;yHMPu5000034973_Danielozyma_ontarioensis;yHMPu5000034974_Deakozyma_indianensis;yHMPu5000034978_Cyberlindnera_mrakii;yHMPu5000034979_Cyberlindnera_misumaiensis;yHMPu5000034986_Candida_oregonensis;yHMPu5000034988_Candida_fructus;yHMPu5000034990_Candida_corydali;yHMPu5000034998_Cephaloascus_albidus;yHMPu5000034999_Cephaloascus_fragrans;yHMPu5000035011_Candida_pyralidae;yHMPu5000035018_Candida_canberraensis;yHMPu5000035022_Candida_emberorum;yHMPu5000035031_Candida_kruisii;yHMPu5000035032_Candida_gatunensis;yHMPu5000035033_Candida_cretensis;yHMPu5000035037_Candida_montana;yHMPu5000035040_Ambrosiozyma_maleeae;yHMPu5000035041_Ambrosiozyma_pseudovanderkliftii;yHMPu5000035044_Barnettozyma_californica;yHMPu5000035045_Barnettozyma_hawaiiensis;yHMPu5000035046_Barnettozyma_populi;yHMPu5000035047_Barnettozyma_pratensis;yHMPu5000035048_Barnettozyma_salicaria;yHMPu5000035242_Zygoascus_ofunaensis;yHMPu5000035243_Zygoascus_meyerae;yHMPu5000035244_Candida_incommunis;yHMPu5000035252_Yamadazyma_nakazawae;yHMPu5000035261_Candida_ponderosae;yHMPu5000035268_Wickerhamomyces_hampshirensis;yHMPu5000035271_Wickerhamomyces_bovis;yHMPu5000035274_Wickerhamomyces_alni;yHMPu5000035279_Tortispora_starmeri;yHMPu5000035282_Trigonopsis_vinaria;yHMPu5000035286_Candida_azyma;yHMPu5000035296_Priceomyces_carsonii;yHMPu5000035297_Priceomyces_castillae;yHMPu5000035301_Pichia_terricola;yHMPu5000035302_Candida_fragi;yHMPu5000035318_Hyphopichia_heimii;yHMPu5000035325_Cyberlindnera_petersonii;yHMPu5000035335_Candida_blattae;yHMPu5000035629_Yueomyces_sinensis;yHMPu5000035633_Candida_hispaniensis;yHMPu5000035639_Wickerhamomyces_canadensis;yHMPu5000035640_Yamadazyma_philogaea;yHMPu5000035641_Yamadazyma_scolyti;yHMPu5000035643_Yarrowia_bubula;yHMPu5000035645_Yarrowia_divulgata;yHMPu5000035650_Trigonopsis_variabilis;yHMPu5000035654_Tortispora_ganteri;yHMPu5000035658_Starmera_amethionina;yHMPu5000035659_Saturnispora_dispora;yHMPu5000035662_Meyerozyma_caribbica;yHMPu5000035665_Middelhovenomyces_tepae;yHMPu5000035667_Kurtzmaniella_cleridarum;yHMPu5000035670_Phaffomyces_opuntiae;yHMPu5000035671_Phaffomyces_antillensis;yHMPu5000035672_Phaffomyces_thermotolerans;yHMPu5000035673_Candida_orba;yHMPu5000035674_Kregervanrija_delftensis;yHMPu5000035675_Kregervanrija_fluxuum;yHMPu5000035677_Kodamaea_ohmeri;yHMPu5000035679_Candida_rhagii;yHMPu5000035681_Candida_gotoi;yHMPu5000035684_Kloeckera_hatyaiensis;yHMPu5000035686_Cyberlindnera_saturnus;yHMPu5000035687_Cyberlindnera_suaveolens;yHMPu5000035688_Cyberlindnera_xylosilytica;yHMPu5000035689_Candida_mycetangii;yHMPu5000035690_Candida_vartiovaarae;yHMPu5000035691_Candida_salmanticensis;yHMPu5000035695_Hanseniaspora_pseudoguilliermondii;yHMPu5000035696_Hanseniaspora_singularis;yHMPu5000035699_Cyberlindnera_maclurae;yHMPu5000035703_Cyberlindnera_americana;yHMPu5000035707_Candida_heveicola;yHMPu5000041678_Debaryomyces_prosopidis;yHMPu5000041693_Debaryomyces_nepalensis;yHMPu5000041713_Debaryomyces_maramus;yHMPu5000041743_Candida_hawaiiana;yHMPu5000041818_Magnusiomyces_tetrasperma;yHMPu5000041822_Dipodascus_geniculatus;yHMPu5000041824_Debaryomyces_subglobosus;yHMPu5000041829_Debaryomyces_fabryi;yHMPu5000041833_Candida_tammaniensis;yHMPu5000041840_Candida_wancherniae;yHMPu5000041855_Candida_ascalaphidarum;yHMPu5000041862_Candida_golubevii;yHMPu5000041863_Candida_gorgasii'
    species_onlyura9 = species_onlyura9.split (';')
    anaerobic = ['Sugiyamaella_lignohabitans' , 'Dekkera_bruxellensis' , 'yHMPu5000034625_Pichia_kudriavzevii' ,
                 'yHMPu5000026142_Citeromyces_matritensis' , 'Candida_albicans' , 'Candida_parapsilosis' ,
                 'Candida_tropicalis' ,
                 'Clavispora_lusitaniae' , 'Spathaspora_passalidarum' , 'Wickerhamia_fluorescens' ,
                 'Wickerhamomyces_anomalus' , 'yHMPu5000035686_Cyberlindnera_saturnus' , 'Hanseniaspora_uvarum' ,
                 'Hanseniaspora_valbyensis' ,
                 'Hanseniaspora_vinae' , 'yHMPu5000034957_Hanseniaspora_osmophila' , 'Ashbya_aceri' ,
                 'Candida_glabrata' ,
                 'Eremothecium_coryli' , 'Kluyveromyces_lactis' , 'Kluyveromyces_marxianus' , 'Lachancea_fermentati' ,
                 'Lachancea_kluyveri' , 'Lachancea_thermotolerans' , 'Lachancea_waltii' , 'Nakaseomyces_bacillisporus' ,
                 'Nakaseomyces_castellii' , 'Nakaseomyces_delphensis' , 'Naumovozyma_castellii' ,
                 'Naumovozyma_dairenensis' ,
                 'Saccharomyces_cerevisiae' , 'Saccharomyces_eubayanus' , 'Saccharomyces_paradoxus' ,
                 'Saccharomyces_uvarum' , 'Tetrapisispora_blattae' , 'Tetrapisispora_phaffii' ,
                 'Torulaspora_delbrueckii' ,
                 'Vanderwaltozyma_polyspora' ,
                 'Zygosaccharomyces_bailii' , 'yHAB154_Kazachstania_transvaalensis' ,
                 'yHMPu5000034881_Torulaspora_pretoriensis' , 'yHMPu5000034876_Tetrapisispora_iriomotensis' ,
                 'yHMPu5000034862_Zygotorulaspora_florentina' ,
                 'yHMPu5000026152_Torulaspora_franciscae' , 'Schizosaccharomyces_pombe']
    species_intersect = [val for val in anaerobic if val in species_onlyura9]
    set_strain_species = [val for val in species_intersect if val in strain]

    if len (set_strain_species):
        model_cobra.reactions.get_by_id ('r_2090').lower_bound = -1000 #uracil uptake

        # 4rd change: Blocked pathways for proper glycerol production
        # Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
        id_0713 = []
        for i in range (len (model['rxns'])):
            if str (model['rxns'][i][0][0]).startswith ('r_0713'):
                id_0713.append (i)
        model_cobra.reactions[id_0713].lower_bound = 0  # Mithocondria
        model_cobra.reactions.get_by_id ('r_0714').lower_bound = 0  # Cytoplasm
        model_cobra.reactions.get_by_id ('r_0713_rvs').lower_bound = 0  # Mithocondria
        model_cobra.reactions.get_by_id ('r_0714_rvs').lower_bound = 0  # Cytoplasm

        # Block glycerol dehydroginase (only acts in microaerobic conditions)
        model_cobra.reactions.get_by_id ('r_0487').upper_bound = 0
        model_cobra.reactions.get_by_id ('r_0487_rvs').upper_bound = 0

        # Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
        model_cobra.reactions.get_by_id ('r_0472').upper_bound = 0
        model_cobra.reactions.get_by_id ('r_0472_fwd').upper_bound = 0

    # 4th change: Blocked pathways for proper glycerol production
    # Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
    idx = []
    for t in range (len (model['rxns'])):
        if str (model['rxns'][t][0][0]).find ('r_0713') >= 0 and str (model['rxns'][t][0][0]).find ('rvs') >= 0:
            idx.append (t)
    for abc in range (len (idx)):
        model_cobra.reactions[idx[abc]].bounds = 0 , 0
    id_0713 = np.where (model['rxns']=='r_0713')
    if len (id_0713[0]) > 0:
        model_cobra.reactions.get_by_id (
            'r_0713').lower_bound = 0  # Mithocondria % in case this one does not have any grRule

    idx2 = []
    for t in range (len (model['rxns'])):
        if str (model['rxns'][t]).find ('r_0714') >= 0 and str (model['rxns'][t]).find ('rvs') >= 0:
            idx2.append (t)
    for cc in range (len (idx2)):
        model_cobra.reactions[idx2[cc]].bounds = 0 , 0
    id_0714 = np.where (model['rxns']=='r_0714')
    if len (id_0714[0]) > 0:
        model_cobra.reactions.get_by_id ('r_0714').lower_bound = 0  # Cytoplasm

    # Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
    idx3 = []
    for s in range (len (model['rxns'])):
        if str (model['rxns'][s][0][0]).startswith ('r_0472_'):
            idx3.append (s)
    for dd in range (len (idx3)):
        model_cobra.reactions[idx3[dd]].bounds = 0 , 0
    model_cobra.reactions.get_by_id ('r_0472').upper_bound = 0
    return model_cobra


def changeGAM(model , model_cobra , GAM , NGAM):
    bioPos = np.where (model['rxnNames']=='biomass pseudoreaction')
    bioPos = bioPos[0][0]

    isGAM_id = ['ATP [cytoplasm]' , 'ADP [cytoplasm]' , 'H2O [cytoplasm]' , 'H+ [cytoplasm]' , 'phosphate [cytoplasm]']
    for i in range (len (isGAM_id)):
        temp = np.where (model['metNames']==isGAM_id[i])
        temp = temp[0][0]
        if model_cobra.reactions[bioPos].get_coefficient (str (model['mets'][temp][0][0]))!=0:
            model_cobra.reactions[bioPos].add_metabolites ({
                                                               str (model['mets'][temp][0][0]): np.sign (
                                                                   model_cobra.reactions[bioPos].get_coefficient (
                                                                       str (model['mets'][temp][0][0]))) * GAM -
                                                                                                model_cobra.reactions[
                                                                                                    bioPos].get_coefficient (
                                                                                                    str (model['mets'][
                                                                                                             temp][0][
                                                                                                             0]))})

    id_non_growth = np.where (model['rxnNames']=='non-growth associated maintenance reaction')
    model_cobra.reactions[id_non_growth[0][0]].bounds = NGAM , NGAM

    return model_cobra


def getrSample(mu, sigma, lb, ub, step, method = 'Uniform'):
    '''
    mu: enzymedata['kcat']
    sigma: enzymedata['kcat_var']
    lb: the lb of kcat
    ub: the ub of kcat
    step: the length of the column of kcat_random_all
    method: way to get kcat
    '''
    from scipy.stats import truncnorm,uniform,norm
    r = np.zeros((len(lb),int(step)))
    for i in range(len(lb)):
        if lb[i] == ub[i] and lb[i] == 0:
            lb[i] = -2
            ub[i] = 8
        if method == 'normal':
            mutmp = np.log10(np.array(mu[i])/3600)
            sigmatmp = sigma[i]
            samples = norm.rvs(loc = mutmp, scale = sigmatmp,size = int(step))
            samples[samples<lb[i]] = lb[i]
            samples[samples>ub[i]] = ub[i]
            samples = 10**samples*3600
            r[i,:] = samples

        elif method == 'Uniform':
            mutmp = np.log10(np.array(mu[i])/3600)
            sigmatmp = sigma[i]
            samples = uniform.rvs(loc = mutmp, scale = sigmatmp,size = int(step))
            samples[samples<lb[i]] = lb[i]
            samples[samples>ub[i]] = ub[i]
            samples = 10**samples*3600
            r[i,:] = samples

    return r

def updateprior(x):
    '''

    x: kcat_100
    :return: a:enzymedata['kcat'](mean), b:enzymedata['kcat_var'](std)
    '''
    from scipy.stats import norm
    x = np.log10(np.array(x)/3600)
    a = []
    b = []
    for i in range (len (x[: , 0])):
        loc , scale = norm.fit (x[i , :])
        a_tmp = 10 ** loc * 3600
        a.append (a_tmp)
        b.append (scale)
    a = np.array(a)
    b = np.array(b)
    return a,b
