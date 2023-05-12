def getFraction(model, data, compType, X):
    import numpy as np
    from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, \
        write_sbml_model
    # define pseudoreaction name
    rxnName = compType+' pseudoreaction'
    rxnName = rxnName.replace('P', 'protein')
    rxnName = rxnName.replace('C', 'carbohydrate')
    rxnName = rxnName.replace('N', 'biomass')
    rxnName = rxnName.replace('L', 'lipid backbone')
    rxnName = rxnName.replace('R', 'RNA')
    rxnName = rxnName.replace('D', 'DNA')
    rxnName = rxnName.replace('I', 'ion')
    rxnName = rxnName.replace('F', 'cofactor')

    # add up fraction
    rxnPos = []
    for i in range(len(model['rxnNames'])):
        if rxnName == model['rxnNames'][i]:
            temp = 1
        else:
            temp = 0
        rxnPos.append(temp)
    index = [i for i, e in enumerate(rxnPos) if e!=0]
    if np.nonzero(rxnPos) != 0:
        sub = model['S']
        isSub = sub.getcol(index[0])<0 #substrates in pseudo-rxn
        if compType == 'L':
            F = - sum(sub[isSub, rxnPoS])  #g/gDW   算出来是空？应为4.2336
        else:
            F = 0
            # add up all components:
            for i in range(len(model['mets'])):
                mets = model['mets']
                pos = []
                if mets[i] in data['mets']:
                    temp = np.where(data['mets']==mets[i])
                    pos.append(temp)
                if isSub[i] and len(pos) > 0:
                    if compType == 'I' or compType == 'F':
                        MW = data['MWs'][pos[0][0]]
                    else:
                        MW = data['MWs'][pos[0][0]] - 18
                    abundance = -sub[i,index]*MW/1000
                    F = F + abundance
        X = X + F

        print(str(compType) + ' -> ' + str(F) + " g/gDW")
    else:
        print(str(compType) + " do not exist")
        F = 0
        X = X + F
    print("X ->"+str(X)+" gDW/DW")

    return X, F

def sumBiomass(model):
    import scipy.io as scio
    from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, \
        write_sbml_model
    import numpy as np
    biomassCompDataFile = "data/biomassCompData.mat"
    biomassCompData = scio.loadmat(biomassCompDataFile)
    biomassCompData = biomassCompData['data'][0, 0]

    P = getFraction(model,biomassCompData,'P',0)
    #C = getFraction(model,biomassCompData,'C',X)
    #R = getFraction(model,biomassCompData,'R',X)
    #D = getFraction(model,biomassCompData,'D',X)
    #L = getFraction(model,biomassCompData,'L',X)
    #I = getFraction(model,biomassCompData,'I',X)
    #F = getFraction(model,biomassCompData,'F',X)

    #print("X ->"+str(X)+" gDW/DW")
    scio.savemat('model__.mat', {'model__': model})
    model__ = load_matlab_model("model__.mat")
    sol = model__.optimize()
    print("Growth = " + str(sol.objective_value) + " 1/h")

    return P

def changeMedia(model, c_source, media, anox = False, flux = -1000):
    import numpy as np
    c_source = c_source + ' exchange'
    # first block any uptake

    # model = setParam(model,'eq',model.rxns(rxnidx),0);
    #exchange = getExchangeRxns(model)
    #idx = model['rxnNames'].find('EX_protein_pool')
    #exchange = set(exchange).difference(set(idx))
    #model['lb'][exchange] = 0
    exchangerxn = np.loadtxt("data/exchangerxns.csv",dtype=str)
    exchangerxn[0] = 'r_1542'
    xx = np.where(model['rxns'] == exchangerxn)
    model['lb'][xx[0]] = 0
    pos = getComponentIndexes(model,c_source)

    # The media will define which rxns to fix:
    if media == 'YEP':
        N = 25  # Aminoacids + Nucleotides
    elif media == 'MAA':
        N = 21  # Aminoacids
    elif media == 'MIN':
        N = 1  # Only the carbon source
    elif media == 'MIN+His':
        N = 1  # Only the carbon source
        id_1893 = np.where(model['rxns'] == 'r_1893')
        model['lb'][id_1893[0][0]] = -0.08 # Histidine exchange
    elif media == 'MIN+Arg':
        N = 1  # Only the carbon source
        id_1879 = np.where(model['rxns'] == 'r_1879')
        model['lb'][id_1879[0][0]] = -0.08 # L-arginine exchange
    elif media == 'MIN+Citrate':
        N = 1  # Only the carbon source
        id_1687 = np.where(model['rxns'] == 'r_1687')
        model['lb'][id_1687[0][0]] = -0.08 # citrate exchange
    # LB parameter (manually optimized for glucose on Min+AA):
    b = -0.08
    # LB parameter (manually optimized for glucose complex media):
    c = -2
    flux = np.array(flux)
    # Define fluxes in case of ec model:
    if N > 1:
        flux = b * np.ones((1, N))
        if N > 21:
            flux[21:25] = c

    #flux = -1000 #flux[0] = -1000

    # Fix values as LBs:
    for i in range(N):
        model['lb'][pos[i]] = -1000 #flux[i]

    # Allow uptake of essential components
    id_1654 = np.where(model['rxns'] == 'r_1654')
    model['lb'][id_1654[0][0]] = -1000 # 'ammonium exchange';
    id_2100 = np.where(model['rxns'] == 'r_2100')
    model['lb'][id_2100[0][0]] = -1000 # 'water exchange';
    id_1861 = np.where(model['rxns'] == 'r_1861')
    model['lb'][id_1861[0][0]] = -1000 # 'iron(2+) exchange';
    id_1992 = np.where(model['rxns'] == 'r_1992')
    model['lb'][id_1992[0][0]] = -1000 # 'oxygen exchange';
    id_2005 = np.where(model['rxns'] == 'r_2005')
    model['lb'][id_2005[0][0]] = -1000 # 'phosphate exchange';
    id_2060 = np.where(model['rxns'] == 'r_2060')
    model['lb'][id_2060[0][0]] = -1000 # 'sulphate exchange';
    id_1832 = np.where(model['rxns'] == 'r_1832')
    model['lb'][id_1832[0][0]] = -1000 # 'H+ exchange';
    id_4593 = np.where(model['rxns'] == 'r_4593')
    model['lb'][id_4593[0][0]] = -1000 # 'chloride exchange';
    id_4595 = np.where(model['rxns'] == 'r_4595')
    model['lb'][id_4595[0][0]] = -1000 # 'Mn(2+) exchange';
    id_4596 = np.where(model['rxns'] == 'r_4596')
    model['lb'][id_4596[0][0]] = -1000 # 'Zn(2+) exchange';
    id_4597 = np.where(model['rxns'] == 'r_4597')
    model['lb'][id_4597[0][0]] = -1000 # 'Mg(2+) exchange';
    id_2049 = np.where(model['rxns'] == 'r_2049')
    model['lb'][id_2049[0][0]] = -1000 # 'sodium exchange';
    id_4594 = np.where(model['rxns'] == 'r_4594')
    model['lb'][id_4594[0][0]] = -1000 # 'Cu(2+) exchange';
    id_4600 = np.where(model['rxns'] == 'r_4600')
    model['lb'][id_4600[0][0]] = -1000 # 'Ca(2+) exchange';
    id_2020 = np.where(model['rxns'] == 'r_2020')
    model['lb'][id_2020[0][0]] = -1000 # 'potassium exchange';

    # Block some production fluxes
    id_1663 = np.where(model['rxns'] == 'r_1663')
    model['ub'][id_1663[0][0]] = 0 # bicarbonate exchange;
    id_4062 = np.where(model['rxns'] == 'r_4062')
    model['ub'][id_4062[0][0]] = 0 # lipid backbone exchange;
    id_4064 = np.where(model['rxns'] == 'r_4064')
    model['ub'][id_4064[0][0]] = 0 # lipid chain exchange;


    # Allow biomass production
    id_2111 = np.where(model['rxns'] == 'r_2111')
    model['ub'][id_2111[0][0]] = 1000 # growth;


    if anox == 'anaerobic':
        1
        model = anerobicModel(model)

    return model

def getComponentIndexes(model,c_source):
    import numpy as np
    c_source = c_source + ' exchange'
    pos = []
    pos.append(np.where(model["rxnNames"] == c_source))
    pos.append(np.where(model["rxnNames"] == 'L-alanine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-arginine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-asparagine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-aspartate exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-cysteine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-glutamine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-glutamate exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-glycine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-histidine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-isoleucine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-leucine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-lysine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-methionine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-phenylalanine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-proline exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-serine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-threonine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-tryptophan exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-tyrosine exchange'))
    pos.append(np.where(model["rxnNames"] == 'L-valine exchange'))
    pos.append(np.where(model["rxnNames"] == '2''-deoxyadenosine exchange'))
    pos.append(np.where(model["rxnNames"] == '2''-deoxyguanosine exchange'))
    pos.append(np.where(model["rxnNames"] == 'thymidine exchange'))
    pos.append(np.where(model["rxnNames"] == 'deoxycytidine exchange'))
    pos.append(np.where(model["rxnNames"] == 'D-glucose exchange'))
    pos = np.array(pos)
    pos = pos[:,0]
    return pos

def anaerobicModel(model):
    import numpy as np
    # 1th change: Refit GAM and NGAM to exp. data, change biomass composition
    GAM = 58.1988
    strain = str(model['id']).replace(' specific model genereted from panYeast', '')
    if strain == 'Candida_glabrata' or strain == 'Candida_parapsilosis':
        GAM = 30

    NGAM = 0
    model = changeGAM(model, GAM, NGAM)
    # 2nd change: Removes the requirement of heme a in the biomass equation
    mets = ['s_3714[c]', 's_1198[c]', 's_1203[c]', 's_1207[c]', 's_1212[c]', 's_0529[c]']
    met_index = []
    for i in range(len(model['mets'])):
        if model['mets'][i] in mets:
            met_index.append(i)
    rxn_index = np.where(model['rxns'] == 'r_4598')
    rxn_index = rxn_index[0][0]
    model['S'][met_index, rxn_index] = 0
    # 3st change: Changes media to anaerobic
    id_1992 = np.where(model['rxns'] == 'r_1992')
    model["lb"][id_1992[0][0]] = 0
    id_1757 = np.where(model['rxns'] == 'r_1757')
    model["lb"][id_1757[0][0]] = -1000
    id_1994 = np.where(model['rxns'] == 'r_1994')
    model["lb"][id_1994[0][0]] = -1000
    id_1915 = np.where(model['rxns'] == 'r_1915')
    model["lb"][id_1915[0][0]] = -1000
    id_2106 = np.where(model['rxns'] == 'r_2106')
    model["lb"][id_2106[0][0]] = -1000
    id_2134 = np.where(model['rxns'] == 'r_2134')
    model["lb"][id_2134[0][0]] = -1000
    id_2137 = np.where(model['rxns'] == 'r_2137')
    model["lb"][id_2137[0][0]] = -1000
    id_2189 = np.where(model['rxns'] == 'r_2189')
    model["lb"][id_2189[0][0]] = -1000

    # extra media set up for ura1 original but also the anaerobic growth media
    if len(model["id"]):
        strain = str(model['id']).replace(' specific model genereted from panYeast', '')
    species_onlyura9 = 'Alloascoidea_hylecoeti;Ambrosiozyma_kashinagacola;Ambrosiozyma_monospora;Arxula_adeninivorans;Ascoidea_asiatica;Ascoidea_rubescens;Ashbya_aceri;Aspergillus_nidulans;Babjeviella_inositovora;Brettanomyces_anomalus;Candida_albicans;Candida_apicola;Candida_arabinofermentans;Candida_auris;Candida_boidinii_JCM9604;Candida_carpophila;Candida_dubliniensis;Candida_glabrata;Candida_homilentoma;Candida_infanticola;Candida_intermedia;Candida_orthopsilosis;Candida_parapsilosis;Candida_sorboxylosa;Candida_succiphila;Candida_tanzawaensis;Candida_tenuis;Candida_tropicalis;Candida_versatilis;Clavispora_lusitaniae;Cyberlindnera_fabianii_JCM3601;Cyberlindnera_jadinii;Debaryomyces_hansenii;Dekkera_bruxellensis;Eremothecium_coryli;Eremothecium_cymbalariae;Eremothecium_gossypii;Eremothecium_sinecaudum;Geotrichum_candidum;Hanseniaspora_uvarum;Hanseniaspora_valbyensis;Hanseniaspora_vinae;Hyphopichia_burtonii;Komagataella_pastoris;Kuraishia_capsulata;Lipomyces_starkeyi;Lodderomyces_elongisporus;Metschnikowia_aberdeeniae;Metschnikowia_arizonensis;Metschnikowia_bicuspidata;Metschnikowia_borealis;Metschnikowia_bowlesiae;Metschnikowia_cerradonensis;Metschnikowia_continentalis;Metschnikowia_dekortum;Metschnikowia_drakensbergensis;Metschnikowia_hamakuensis;Metschnikowia_hawaiiensis;Metschnikowia_hibisci;Metschnikowia_ipomoeae;Metschnikowia_kamakouana;Metschnikowia_kipukae;Metschnikowia_lockheadii;Metschnikowia_matae;Metschnikowia_matae_maris;Metschnikowia_mauinuiana;Metschnikowia_proteae;Metschnikowia_santaceciliae;Metschnikowia_shivogae;Metschnikowia_similis;Meyerozyma_guilliermondii;Millerozyma_acaciae;Nadsonia_fulvescens_var_elongata;Nakaseomyces_bracarensis;Nakaseomyces_castellii;Nakaseomyces_delphensis;Nakaseomyces_nivariensis;Nakazawaea_peltata;Ogataea_methanolica;Ogataea_parapolymorpha;Ogataea_polymorpha;Pachysolen_tannophilus;Pichia_membranifaciens;Priceomyces_haplophilus;Saccharomycopsis_malanga;Saprochaete_clavata;Scheffersomyces_lignosus;Scheffersomyces_stipitis;Schizosaccharomyces_pombe;Spathaspora_arborariae;Spathaspora_girioi;Spathaspora_gorwiae;Spathaspora_hagerdaliae;Spathaspora_passalidarum;Sporopachydermia_quercuum;Starmerella_bombicola_JCM9596;Sugiyamaella_lignohabitans;Tortispora_caseinolytica;Vanderwaltozyma_polyspora;Wickerhamia_fluorescens;Wickerhamiella_domercqiae;Wickerhamomyces_anomalus;Wickerhamomyces_ciferrii;Yarrowia_deformans;Yarrowia_keelungensis;Yarrowia_lipolytica;yHMPu5000026124_Ogataea_henricii;yHMPu5000026137_Ambrosiozyma_ambrosiae;yHMPu5000026142_Citeromyces_matritensis;yHMPu5000026145_Ambrosiozyma_vanderkliftii;yHMPu5000026197_Brettanomyces_custersianus;yHMPu5000026274_Komagataella_populi;yHMPu5000034594_Starmera_quercuum;yHMPu5000034597_Candida_stellimalicola;yHMPu5000034604_Sporopachydermia_lactativora;yHMPu5000034605_Spencermartinsiella_europaea;yHMPu5000034606_Priceomyces_medius;yHMPu5000034607_Saccharomycopsis_capsularis;yHMPu5000034610_Saturnispora_hagleri;yHMPu5000034611_Saturnispora_mendoncae;yHMPu5000034612_Saturnispora_saitoi;yHMPu5000034613_Saturnispora_serradocipensis;yHMPu5000034614_Saturnispora_silvae;yHMPu5000034615_Saturnispora_zaruensis;yHMPu5000034622_Pichia_occidentalis;yHMPu5000034623_Pichia_norvegensis;yHMPu5000034624_Pichia_nakasei;yHMPu5000034625_Pichia_kudriavzevii;yHMPu5000034627_Pichia_heedii;yHMPu5000034629_Pichia_exigua;yHMPu5000034631_Martiniozyma_abiesophila;yHMPu5000034632_Candida_athensensis;yHMPu5000034635_Nadsonia_fulvescens;yHMPu5000034636_Ogataea_nitratoaversa;yHMPu5000034637_Ogataea_populiabae;yHMPu5000034643_Candida_schatavii;yHMPu5000034646_Wickerhamiella_cacticola;yHMPu5000034648_Candida_restingae;yHMPu5000034654_Aciculoconidium_aculeatum;yHMPu5000034655_Botryozyma_nematodophila;yHMPu5000034660_Diddensiella_caesifluorescens;yHMPu5000034661_Dipodascus_albidus;yHMPu5000034665_Kodamaea_laetipori;yHMPu5000034667_Blastobotrys_serpentis;yHMPu5000034669_Blastobotrys_raffinofermentans;yHMPu5000034670_Blastobotrys_proliferans;yHMPu5000034671_Blastobotrys_peoriensis;yHMPu5000034673_Blastobotrys_nivea;yHMPu5000034674_Blastobotrys_muscicola;yHMPu5000034675_Blastobotrys_mokoenaii;yHMPu5000034681_Blastobotrys_americana;yHMPu5000034742_Lipomyces_suomiensis;yHMPu5000034748_Lipomyces_oligophaga;yHMPu5000034749_Lipomyces_mesembrius;yHMPu5000034754_Lipomyces_arxii;yHMPu5000034760_Lipomyces_kononenkoae;yHMPu5000034761_Lipomyces_lipofer;yHMPu5000034883_Peterozyma_xylosa;yHMPu5000034884_Peterozyma_toletana;yHMPu5000034885_Ogataea_zsoltii;yHMPu5000034886_Ogataea_trehalophila;yHMPu5000034887_Ogataea_trehaloabstinens;yHMPu5000034890_Ogataea_ramenticola;yHMPu5000034891_Ogataea_pini;yHMPu5000034892_Ogataea_pilisensis;yHMPu5000034893_Ogataea_philodendra;yHMPu5000034897_Ogataea_glucozyma;yHMPu5000034899_Ogataea_kodamae;yHMPu5000034901_Ogataea_methylivora;yHMPu5000034902_Ogataea_minuta;yHMPu5000034903_Ogataea_naganishii;yHMPu5000034904_Ogataea_nonfermentans;yHMPu5000034918_Nakazawaea_holstii;yHMPu5000034933_Kuraishia_molischiana;yHMPu5000034939_Komagataella_pseudopastoris;yHMPu5000034946_Ambrosiozyma_oregonensis;yHMPu5000034947_Ambrosiozyma_philentoma;yHMPu5000034950_Citeromyces_hawaiiensis;yHMPu5000034952_Citeromyces_siamensis;yHMPu5000034957_Hanseniaspora_osmophila;yHMPu5000034963_Hanseniaspora_clermontiae;yHMPu5000034967_Candida_freyschussii;yHMPu5000034973_Danielozyma_ontarioensis;yHMPu5000034974_Deakozyma_indianensis;yHMPu5000034978_Cyberlindnera_mrakii;yHMPu5000034979_Cyberlindnera_misumaiensis;yHMPu5000034986_Candida_oregonensis;yHMPu5000034988_Candida_fructus;yHMPu5000034990_Candida_corydali;yHMPu5000034998_Cephaloascus_albidus;yHMPu5000034999_Cephaloascus_fragrans;yHMPu5000035011_Candida_pyralidae;yHMPu5000035018_Candida_canberraensis;yHMPu5000035022_Candida_emberorum;yHMPu5000035031_Candida_kruisii;yHMPu5000035032_Candida_gatunensis;yHMPu5000035033_Candida_cretensis;yHMPu5000035037_Candida_montana;yHMPu5000035040_Ambrosiozyma_maleeae;yHMPu5000035041_Ambrosiozyma_pseudovanderkliftii;yHMPu5000035044_Barnettozyma_californica;yHMPu5000035045_Barnettozyma_hawaiiensis;yHMPu5000035046_Barnettozyma_populi;yHMPu5000035047_Barnettozyma_pratensis;yHMPu5000035048_Barnettozyma_salicaria;yHMPu5000035242_Zygoascus_ofunaensis;yHMPu5000035243_Zygoascus_meyerae;yHMPu5000035244_Candida_incommunis;yHMPu5000035252_Yamadazyma_nakazawae;yHMPu5000035261_Candida_ponderosae;yHMPu5000035268_Wickerhamomyces_hampshirensis;yHMPu5000035271_Wickerhamomyces_bovis;yHMPu5000035274_Wickerhamomyces_alni;yHMPu5000035279_Tortispora_starmeri;yHMPu5000035282_Trigonopsis_vinaria;yHMPu5000035286_Candida_azyma;yHMPu5000035296_Priceomyces_carsonii;yHMPu5000035297_Priceomyces_castillae;yHMPu5000035301_Pichia_terricola;yHMPu5000035302_Candida_fragi;yHMPu5000035318_Hyphopichia_heimii;yHMPu5000035325_Cyberlindnera_petersonii;yHMPu5000035335_Candida_blattae;yHMPu5000035629_Yueomyces_sinensis;yHMPu5000035633_Candida_hispaniensis;yHMPu5000035639_Wickerhamomyces_canadensis;yHMPu5000035640_Yamadazyma_philogaea;yHMPu5000035641_Yamadazyma_scolyti;yHMPu5000035643_Yarrowia_bubula;yHMPu5000035645_Yarrowia_divulgata;yHMPu5000035650_Trigonopsis_variabilis;yHMPu5000035654_Tortispora_ganteri;yHMPu5000035658_Starmera_amethionina;yHMPu5000035659_Saturnispora_dispora;yHMPu5000035662_Meyerozyma_caribbica;yHMPu5000035665_Middelhovenomyces_tepae;yHMPu5000035667_Kurtzmaniella_cleridarum;yHMPu5000035670_Phaffomyces_opuntiae;yHMPu5000035671_Phaffomyces_antillensis;yHMPu5000035672_Phaffomyces_thermotolerans;yHMPu5000035673_Candida_orba;yHMPu5000035674_Kregervanrija_delftensis;yHMPu5000035675_Kregervanrija_fluxuum;yHMPu5000035677_Kodamaea_ohmeri;yHMPu5000035679_Candida_rhagii;yHMPu5000035681_Candida_gotoi;yHMPu5000035684_Kloeckera_hatyaiensis;yHMPu5000035686_Cyberlindnera_saturnus;yHMPu5000035687_Cyberlindnera_suaveolens;yHMPu5000035688_Cyberlindnera_xylosilytica;yHMPu5000035689_Candida_mycetangii;yHMPu5000035690_Candida_vartiovaarae;yHMPu5000035691_Candida_salmanticensis;yHMPu5000035695_Hanseniaspora_pseudoguilliermondii;yHMPu5000035696_Hanseniaspora_singularis;yHMPu5000035699_Cyberlindnera_maclurae;yHMPu5000035703_Cyberlindnera_americana;yHMPu5000035707_Candida_heveicola;yHMPu5000041678_Debaryomyces_prosopidis;yHMPu5000041693_Debaryomyces_nepalensis;yHMPu5000041713_Debaryomyces_maramus;yHMPu5000041743_Candida_hawaiiana;yHMPu5000041818_Magnusiomyces_tetrasperma;yHMPu5000041822_Dipodascus_geniculatus;yHMPu5000041824_Debaryomyces_subglobosus;yHMPu5000041829_Debaryomyces_fabryi;yHMPu5000041833_Candida_tammaniensis;yHMPu5000041840_Candida_wancherniae;yHMPu5000041855_Candida_ascalaphidarum;yHMPu5000041862_Candida_golubevii;yHMPu5000041863_Candida_gorgasii'
    species_onlyura9 = species_onlyura9.split(';')
    anaerobic = ['Sugiyamaella_lignohabitans', 'Dekkera_bruxellensis', 'yHMPu5000034625_Pichia_kudriavzevii',
                 'yHMPu5000026142_Citeromyces_matritensis', 'Candida_albicans', 'Candida_parapsilosis',
                 'Candida_tropicalis',
                 'Clavispora_lusitaniae', 'Spathaspora_passalidarum', 'Wickerhamia_fluorescens',
                 'Wickerhamomyces_anomalus', 'yHMPu5000035686_Cyberlindnera_saturnus', 'Hanseniaspora_uvarum',
                 'Hanseniaspora_valbyensis',
                 'Hanseniaspora_vinae', 'yHMPu5000034957_Hanseniaspora_osmophila', 'Ashbya_aceri', 'Candida_glabrata',
                 'Eremothecium_coryli', 'Kluyveromyces_lactis', 'Kluyveromyces_marxianus', 'Lachancea_fermentati',
                 'Lachancea_kluyveri', 'Lachancea_thermotolerans', 'Lachancea_waltii', 'Nakaseomyces_bacillisporus',
                 'Nakaseomyces_castellii', 'Nakaseomyces_delphensis', 'Naumovozyma_castellii',
                 'Naumovozyma_dairenensis',
                 'Saccharomyces_cerevisiae', 'Saccharomyces_eubayanus', 'Saccharomyces_paradoxus',
                 'Saccharomyces_uvarum', 'Tetrapisispora_blattae', 'Tetrapisispora_phaffii', 'Torulaspora_delbrueckii',
                 'Vanderwaltozyma_polyspora',
                 'Zygosaccharomyces_bailii', 'yHAB154_Kazachstania_transvaalensis',
                 'yHMPu5000034881_Torulaspora_pretoriensis', 'yHMPu5000034876_Tetrapisispora_iriomotensis',
                 'yHMPu5000034862_Zygotorulaspora_florentina',
                 'yHMPu5000026152_Torulaspora_franciscae', 'Schizosaccharomyces_pombe']
    species_intersect = [val for val in anaerobic if val in species_onlyura9]
    set_strain_species = [val for val in species_intersect if val in strain]

    if len(set_strain_species):
        id_2090 = np.where(model['rxns'] == 'r_2090')
        model["lb"][id_2090[0][0]] = -1000

    # 4rd change: Blocked pathways for proper glycerol production
        # Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
        id_0713 = []
        for i in range(len(model['rxns'])):
            if str(model['rxns'][i][0][0]).startswith('r_0713'):
                id_0713.append(i)
        model['lb'][id_0713] = 0 # Mithocondria
        id_0714 = np.where(model['rxns'] == 'r_0714')
        model["lb"][id_0714[0][0]] = 0 # Cytoplasm
        id_0713_rvs = np.where(model['rxns'] == 'r_0713_rvs')
        model["lb"][id_0713_rvs[0][0]] = 0 # Mithocondria
        id_0714_rvs = np.where(model['rxns'] == 'r_0714_rvs')
        model["lb"][id_0714_rvs[0][0]] = 0 # Cytoplasm

        # Block glycerol dehydroginase (only acts in microaerobic conditions)
        id_0487 = np.where(model['rxns'] == 'r_0487')
        model["ub"][id_0487[0][0]] = 0
        id_0487_rvs = np.where(model['rxns'] == 'r_0487_rvs')
        model["ub"][id_0487_rvs[0][0]] = 0

        # Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
        id_0472 = np.where(model['rxns'] == 'r_0472')
        model["ub"][id_0472[0][0]] = 0
        id_0472_fwd = np.where(model['rxns'] == 'r_0472_fwd')
        model["ub"][id_0472_fwd[0][0]] = 0

    # 4th change: Blocked pathways for proper glycerol production
    # Block oxaloacetate-malate shuttle (not present in anaerobic conditions)
    idx = []
    for t in range(len(model['rxns'])):
        if str(model['rxns'][t][0][0]).find('r_0713') >= 0 and str(model['rxns'][t][0][0]).find('rvs') >= 0:
            idx.append(t)
    model['ub'][idx] = 0
    model['lb'][idx] = 0
    id_00713 = np.where(model['rxns'] == 'r_0713')
    if len(id_00713[0]) != 0:
        model['lb'][id_00713[0][0]] = 0 #Mithocondria % in case this one does not have any grRule

    idx2 = []
    for t in range(len(model['rxns'])):
        if str(model['rxns'][t]).find('r_0714') >= 0 and str(model['rxns'][t]).find('rvs') >= 0:
            idx2.append(t)
    model['ub'][idx2] = 0
    model['lb'][idx2] = 0
    id_00714 = np.where(model['rxns'] == 'r_0714')
    if len(id_00714[0]) != 0:
        model['lb'][id_00714[0][0]] = 0 #Cytoplasm

    # %Block 2-oxoglutarate + L-glutamine -> 2 L-glutamate (alternative pathway)
    idx3 = []
    for s in range(len(model['rxns'])):
        if str(model['rxns'][s][0][0]).startswith('r_0472_'):
            idx3.append(s)
    model['ub'][idx3] = 0
    model['lb'][idx3] = 0
    id_00472 = np.where(model['rxns'] == 'r_0472')
    model['ub'][id_00472[0][0]] = 0
    return model

def changeGAM(model, GAM, NGAM):
    import numpy as np
    bioPos = np.where(model['rxnNames'] == 'biomass pseudoreaction')
    bioPos = bioPos[0][0]

    for i in range(len(model['metNames'])):
        S_ix = model['S'][i, bioPos]
        isGAM_id = ['ATP [cytoplasm]', 'ADP [cytoplasm]', 'H2O [cytoplasm]', 'H+ [cytoplasm]', 'phosphate [cytoplasm]']
        if str(model['metNames'][i]) in isGAM_id:
            isGAM = True
        else:
            isGAM = False

        if S_ix!=0 and isGAM:
            if S_ix > 0:
                model['S'][i, bioPos] = GAM
            else:
                model['S'][i, bioPos] = -GAM

    id_non_growth = np.where(model['rxnNames'] == 'non-growth associated maintenance reaction')
    model["lb"][id_non_growth[0][0]] = NGAM
    model["ub"][id_non_growth[0][0]] = NGAM

    return model

def rmsecal(model, data, constrain, objective, osenseStr,prot_cost_info_id,prot_cost_info_value, tot_prot_weight, rxn2block):
    import numpy as np
    import scipy.io as scio
    import cobra
    from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, \
        write_sbml_model
    from sklearn.metrics import mean_squared_error
    data = np.array(data)
    rmse_tmp = []
    simulated = []
    rxnNames = model['rxnNames']
    exp = np.array(data[:, 2:11])
    exp = exp * [1, -1, 1, 1, 1, 1, 1, 1, -1]
    exp = exp.astype(np.float)

    for i in range(len(data[:, 0])):
        ex_mets = ['growth', data[i, 1][0] + " exchange", 'acetate exchange', 'ethanol exchange',
                   'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
                   'oxygen exchange']
        idx = []
        temp = []
        for k in range(len(rxnNames)):
            if rxnNames[k] in ex_mets:
                temp.append(k)
        idx.append(temp)
        model_tmp = model

        model_tmp = changeMedia(model_tmp, data[i, 1][0], data[i, 15]) #'D-glucose'

        id_1634 = np.where(model_tmp['rxns']=='r_1634')
        model_tmp['ub'][id_1634[0][0]] = 0
        model_tmp['lb'][id_1634[0][0]] = 0
        id_1631 = np.where(model_tmp['rxns']=='r_1631')
        model_tmp['ub'][id_1631[0][0]] = 0
        model_tmp['lb'][id_1631[0][0]] = 0

        if data[i, 13]=="anaerobic" or data[i, 13]=="limited":
            model_tmp = anaerobicModel(model_tmp)

        if data[i, 13]=="limited" :
            id_EX_O2 = np.where(model_tmp['rxnNames']=='oxygen exchange')
            model_tmp['lb'][id_EX_O2[0][0]] = -5

        id_1714 = np.where(model_tmp['rxns']=='r_1714')
        model_tmp['lb'][id_1714[0][0]] = 0

        id_substrate = np.where(model_tmp['rxnNames']==(data[i, 1][0] + " exchange"))
        if constrain == False:
            model_tmp['lb'][id_substrate[0][0]] = -1000
        else:
            model_tmp['lb'][id_substrate[0][0]] = exp[:, 1][i]
        print("No. " + str(i + 1) + "finish !")
        scio.savemat('modeltemp.mat', {'model': model_tmp})
        # save_matlab_model(model_tmp,"modeltemp.mat")
        modeltemp = load_matlab_model("modeltemp.mat")
        # Add protein cost information
        id_obj = np.where(model_tmp['rxns'] == objective)
        model_tmp['c'][id_obj[0][0]] = 1
        nMets = model_tmp['S'].shape[0]
        nRxns = model_tmp['S'].shape[1]
        cost_list = np.zeros((1,nRxns))
        for p in range(nRxns):
            rxnid = model_tmp['rxns'][p]
            id_tmp = rxnid
            id_temp = np.where(prot_cost_info_id == id_tmp)
            if len(np.array(id_temp)[0]) > 0:
                cost = prot_cost_info_value[id_temp[0][0]]
            else:
                cost = 0
            cost_list[0,p] = cost
        cost_list_new = np.transpose(cost_list)
        met = cobra.Metabolite('cost',name = 'cost',compartment= 'c')
        modeltemp.add_metabolites(met)
        tupll = tuple([float(m) for m in cost_list_new])
        for q,reaction in enumerate(modeltemp.reactions):
            reaction.add_metabolites({'cost': tupll[q]})
        modeltemp.add_boundary(modeltemp.metabolites.get_by_id("cost"),type="sink",ub = tot_prot_weight*1000)


        sol_tmp = modeltemp.optimize()
        print("No. " + str(i + 1) + " result is " + str(sol_tmp.objective_value))
        sol = sol_tmp.fluxes  # sol[:,i] = sol_tmp.fluxes

        tmp = np.where(~np.isnan(exp[i]))[0]
        excarbon = model_tmp['excarbon'][:, idx[0]]
        for x in range(len(idx[0])):
            if excarbon[:, x]==0:
                excarbon[:, x] = 1
        exp_tmp = []
        for s in range(len(tmp)):
            exp_tmp.append(exp[i, tmp[s]] * excarbon[:, tmp[s]])
        print("exp_tmp is " + str(exp_tmp))
        sol_idx = np.transpose(sol[idx[0]])
        sol_idx = np.array(sol_idx)
        simulated_tmp = []
        for t in range(len(tmp)):
            simulated_tmp.append(np.array(sol_idx)[tmp[t]]*excarbon[:,tmp[t]]) # normalize the growth rate issue by factor 10

        exp_block = np.zeros((1, 218))  # 218怎么替换
        rxnblockidx_pre = np.setdiff1d(rxn2block, model_tmp['rxns'][idx[0][1]])
        rxnblockidx = []
        for k in range(len(model_tmp['rxns'])):
            if model_tmp['rxns'][k] in rxnblockidx_pre:
                rxnblockidx.append(k)
        simulated_block = []
        soll = np.array(sol)
        for x in range(len(rxnblockidx)):
            simulated_block.append(model_tmp['excarbon'][:, rxnblockidx[x]] * soll[rxnblockidx[x]])
        id_zero = np.where(np.array(simulated_block)!=0)
        exp_block = exp_block[:, id_zero[0]]
        simulated_block = np.array(simulated_block)[id_zero[0]]
        exp_tmp = np.array(exp_tmp)

        if constrain:
            rmse_tmp.append(np.sqrt(mean_squared_error(np.append(exp_tmp,np.transpose(exp_block)),np.append(simulated_tmp,np.transpose(simulated_block)))))
        else:
            if len(exp_tmp) >= 2:
                rmse_tmp.append(np.sqrt(mean_squared_error(exp_tmp[0:2], np.array(simulated_tmp)[0:2])))
            else:
                rmse_tmp.append(np.sqrt(mean_squared_error(exp_tmp[0],[simulated_tmp[0]])))
        print(rmse_tmp)
        simulated.append(np.transpose(soll[idx[0]]))

    rmse = sum(rmse_tmp) / len(data[:, 0])
    print("RMSE = " + str(rmse))

    #return rmse, exp, simulated
    return rmse,exp,simulated

def abc_python_max(model, enzymedata, kcat_random_all, tot_prot_weight, growthdata, max_growth, proc, sample_generation,
                   j, rxn2block):
    import numpy as np
    nstep = sample_generation / proc
    nstep = int(nstep)
    rmse_final = np.zeros((1, nstep))
    kcat_sample = kcat_random_all[:, ((j-1)* nstep ):(j * nstep)]


    # get carbonnum for each exchange rxn to further calculation of error
    if len(model["excarbon"]) == 0:  ##检查model中是否有excarbon, ~is_field
        model = addCarbonNum(model)  ##补写function
    for k in range(nstep):
        print('nstep:' + str(k + 1) + '/' + str(nstep))
        kcat_random = kcat_sample[:,k]
        prot_cost_info_value = [enzymedata["MW"][i] / kcat_random[i] for i in range(len(kcat_random))]
        prot_cost_info_id = enzymedata["rxn_list"]

        # first search with substrate constrain
        objective = 'r_2111'
        osenseStr = 'max'
        if len(growthdata) != 0:
            rmse_1,exp_1,simulated_1 = rmsecal (model, growthdata, True, objective, osenseStr,prot_cost_info_id,prot_cost_info_value, tot_prot_weight, rxn2block)
            #result_1 = np.array(result_1)
            print(rmse_1)
            #rmse_1 = result_1[0]
           # exp_1 = result_1[1]
            #simulated_1 = result_1[2]
        else:
            rmse_1 = []
            exp_1 = []
            simulated_1 = []
        print("RMSE_1 finish !")
        # second searcch for maxmial growth rate without constrain
        if len(max_growth) != 0:
            rmse_2,exp_2,simulated_2 = rmsecal (model, max_growth, False, objective, osenseStr,prot_cost_info_id,prot_cost_info_value, tot_prot_weight, rxn2block)
            #result_2 = np.array(result_2)
            #rmse_2 = result_2[0]
            #exp_2 = result_2[1]
            #simulated_2 = result_2[2]
            print(rmse_2)
        else:
            rmse_2 = []
            exp_2 = []
            simulated_2 = []

        print("RMSE_2 finish !")

        exp = np.array([exp_1, exp_2])
        simulated = np.array([simulated_1, simulated_2])
        rmse = np.array([rmse_1, rmse_2])
        rmse_final = np.nanmean(rmse, 0)

        # only output simulated result for one generation
        if nstep != 1 or sample_generation != 1:
            simulated = []
            exp = []

    #return rmse_final, exp, simulated, growthdata, max_growth,result_1,result_2
    return rmse_final,exp,simulated

def getrSample(mu, sigma, lb, ub, step, method = 'Uniform'):
    import numpy as np
    r = []
    for i in range(len(lb)):
        if lb[i] == ub[i] and lb[i] == 0:
            lb[i] = -2
            ub[i] = 9
        if method == 'normal':
            mutmp = np.log10(np.array(mu)/3600)
            sigmatmp = sigma
            pdd = np.random.normal(loc=mutmp, scale=sigmatmp)
            r_tmp = np.random.choice(pdd, size = (1,int(step)))
            if np.any(r_tmp) < lb[i]:
                r_tmp[r_tmp<lb[i]] = lb[i]
            if np.any(r_tmp) > ub[i]:
                r_tmp[r_tmp>ub[i]] = ub[i]
            r.append(r_tmp)
        elif method == 'Uniform':
            mutmp = np.log10(np.array(mu)/3600)
            sigmatmp = sigma
            pdd = np.random.uniform(low=mutmp - sigmatmp, high=mutmp + sigmatmp)
            t = np.clip(pdd, -2, 8)
            #t = np.random.uniform(pd,low=-2, high=8)
            r_tmp = np.random.choice(t, size = (1,int(step)))
            r.append(r_tmp)
    r = np.array(r)[:,0]
    return r

def updateprior(x):
    from distfit import distfit
    import numpy as np
    x = np.log10(np.array(x)/3600)
    dist = distfit(distr = 'norm')
    a = []
    b = []
    for i in range(len(x[:,0])):
        results = dist.fit_transform(x[i,:])
        mu = dist.model['loc']
        a_tmp = 10**mu*3600
        sigma = dist.model['scale']
        a.append(a_tmp)
        b.append(sigma)
    a = np.array(a)
    b = np.array(b)
    return a,b
