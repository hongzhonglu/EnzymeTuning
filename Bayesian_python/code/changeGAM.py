def changeGAM(model, GAM, NGAM):
    bioPos = str(model['rxnNames']).find("biomass pseudoreaction")
    isGAM = ['ATP [cytoplasm]', 'ADP [cytoplasm]', 'H2O [cytoplasm]', 'H+ [cytoplasm]', 'phosphate [cytoplasm]']
    for i in range(len(model['mets'])):
        S_ix = model['S'][i, bioPos]
        if S_ix != 0 and isGAM in model['metNames[i]']:
            if S_ix > 0:
                s_ix = 1
            else:
                s_ix = -1
            model['S'][i, bioPos] = GAM * s_ix

    pos = str(model['rxnNames']).find("non-growth associated maintenance reaction")
    model['lb']['pos'] = NGAM
    model['Ub']['pos'] = NGAM