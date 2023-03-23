def changeMedia(model, c_source, media, anox = False, flux = -1000):
    c_source = c_source + ' exchange'
    # first block any uptake

    # model = setParam(model,'eq',model.rxns(rxnidx),0);
    #exchange = getExchangeRxns(model)
    #idx = model['rxnNames'].find('EX_protein_pool')
    #exchange = set(exchange).difference(set(idx))
    #model['lb'][exchange] = 0
    pos = []
    c_source = 'D-glucose exchange'
    for i in range(len(model['rxnNames'])):
        if model['rxnNames'][i] == c_source:
            pos.append(i)

    # The media will define which rxns to fix:
    if media == 'YEP':
        N = 25  # Aminoacids + Nucleotides
    elif media == 'MAA':
        N = 21  # Aminoacids
    elif media == 'MIN':
        N = 1  # Only the carbon source
    elif media == 'MIN+His':
        N = 1  # Only the carbon source
        model.reactions.get_by_id('r_1893').lower_bound = -0.08  # Histidine exchange
    elif media == 'MIN+Arg':
        N = 1  # Only the carbon source
        model.reactions.get_by_id('r_1879').lower_bound = -0.08  # L-arginine exchange
    elif media == 'MIN+Citrate':
        N = 1  # Only the carbon source
        model.reactions.get_by_id('r_1687').lower_bound = -0.08  # citrate exchange
    # LB parameter (manually optimized for glucose on Min+AA):
    b = -0.08
    # LB parameter (manually optimized for glucose complex media):
    c = -2
    flux = np.array(flux)
    # Define fluxes in case of ec model:
    if N > 1:
        flux = b * np.ones((1, N))
        if N > 21:
            flux[21:25] = c;
    flux = -1000 #flux[0] = -1000

    # Fix values as LBs:
    for i in range(N):
        model['lb'][pos[i]] = flux #flux[i]

    # Allow uptake of essential components
    for i in range(len(model['rxns'])):
        if model['rxns'][i] == 'r_1654':
            model['lb'][i] = -1000 # 'ammonium exchange';
        if model['rxns'][i] == 'r_2100':
            model['lb'][i] = -1000 # 'water exchange' ;
        if model['rxns'][i] == 'r_1861':
            model['lb'][i] = -1000 # 'iron(2+) exchange';
        if model['rxns'][i] == 'r_1992':
            model['lb'][i] = -1000 # 'oxygen exchange';
        if model['rxns'][i] == 'r_2005':
            model['lb'][i] = -1000 # 'phosphate exchange';
        if model['rxns'][i] == 'r_2060':
            model['lb'][i] = -1000 # 'sulphate exchange';
        if model['rxns'][i] == 'r_1832':
            model['lb'][i] = -1000 # 'H+ exchange' ;
        if model['rxns'][i] == 'r_4593':
            model['lb'][i] = -1000 # 'chloride exchange' ;
        if model['rxns'][i] == 'r_4595':
            model['lb'][i] = -1000 # Mn(2+) exchange
        if model['rxns'][i] == 'r_4596':
            model['lb'][i] = -1000 # Zn(2+ exchange
        if model['rxns'][i] == 'r_4597':
            model['lb'][i] = -1000 # Mg(2+) exchange
        if model['rxns'][i] == 'r_2049':
            model['lb'][i] = -1000 # sodium exchange
        if model['rxns'][i] == 'r_4594':
            model['lb'][i] = -1000 # Cu(2+) exchange
        if model['rxns'][i] == 'r_4600':
            model['lb'][i] = -1000 # Ca(2+) exchange
        if model['rxns'][i] == 'r_2020':
            model['lb'][i] = -1000 # potassium exchange

    # Block some production fluxes
        if model['rxns'][i] == 'r_1663':
            model['ub'][i] = 0 # bicarbonate exchange
        if model['rxns'][i] == 'r_4062':
            model['ub'][i] = 0 # lipid backbone exchange
        if model['rxns'][i] == 'r_4064':
            model['ub'][i] = 0 # lipid chain exchange

    # Allow biomass production
        if model['rxns'][i] == 'r_2111':
            model['ub'][i] = 1000 # growth

    if anox == 'anaerobic':
        1
        model = anerobicModel(model)

    return model
