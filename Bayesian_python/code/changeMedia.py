def changeMedia(model, c_source, media, anox, flux):
    c_source = c_source + ' exchange'
    # first block any uptake

    # model = setParam(model,'eq',model.rxns(rxnidx),0);
    exchange = getExchangeRxns(model)
    idx = model['rxnNames'].find('EX_protein_pool')
    exchange = set(exchange).difference(set(idx))
    model['lb'][exchange] = 0
    pos = getComponentIndexes(model, c_source)

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

    # Define fluxes in case of ec model:
    if N > 1:
        flux = b * np.ones((1, N))
        if N > 21:
            flux[21:25] = c;
    flux[0] = -1000

    # Fix values as LBs:
    for i in range(N):
        model['lb'][pos[i]] = flux[i]

    # Allow uptake of essential components
    model.reactions.get_by_id('r_1654').lower_bound = -1000  # 'ammonium exchange';
    model.reactions.get_by_id('r_2100').lower_bound = -1000  # 'water exchange' ;
    model.reactions.get_by_id('r_1861').lower_bound = -1000  # 'iron(2+) exchange';
    model.reactions.get_by_id('r_1992').lower_bound = -1000  # 'oxygen exchange';
    model.reactions.get_by_id('r_2005').lower_bound = -1000  # 'phosphate exchange';
    model.reactions.get_by_id('r_2060').lower_bound = -1000  # 'sulphate exchange';
    model.reactions.get_by_id('r_1832').lower_bound = -1000  # 'H+ exchange' ;
    model.reactions.get_by_id('r_4593').lower_bound = -1000  # 'chloride exchange' ;
    model.reactions.get_by_id('r_4595').lower_bound = -1000  # Mn(2+) exchange
    model.reactions.get_by_id('r_4596').lower_bound = -1000  # Zn(2+ exchange
    model.reactions.get_by_id('r_4597').lower_bound = -1000  # Mg(2+) exchange
    model.reactions.get_by_id('r_2049').lower_bound = -1000  # sodium exchange
    model.reactions.get_by_id('r_4594').lower_bound = -1000  # Cu(2+) exchange
    model.reactions.get_by_id('r_4600').lower_bound = -1000  # Ca(2+) exchange
    model.reactions.get_by_id('r_2020').lower_bound = -1000  # potassium exchange
    # Block some production fluxes
    model.reactions.get_by_id('r_1663').upper_bound = 0  # bicarbonate exchange
    model.reactions.get_by_id('r_4062').upper_bound = 0  # lipid backbone exchange
    model.reactions.get_by_id('r_4064').upper_bound = 0  # lipid chain exchange
    # Allow biomass production
    model.reactions.get_by_id('r_2111').upper_bound = 1000  # growth

    if anox == 'anaerobic':
        1
        model = anerobicModel(model)