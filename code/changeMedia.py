def changeMedia(model, c_source, media, anox = False, flux = -1000):
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