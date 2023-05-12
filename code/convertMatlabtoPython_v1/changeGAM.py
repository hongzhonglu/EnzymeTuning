def changeGAM(model, GAM, NGAM):
    bioPos = []
    for i in range(len(rxnNames)):
        if rxnNames[i]=='biomass pseudoreaction':
            bioPos.append(i)

    for i in range(len(model['mets'])):
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

    for i in range(len(rxnNames)):
        if rxnNames[i]=='non-growth associated maintenance reaction':
            model_tmp["lb"][i] = NGAM
            model_tmp["ub"][i] = NGAM
