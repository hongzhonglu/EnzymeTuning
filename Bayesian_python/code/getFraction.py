def getFraction(model, data, compType, X):
    # define pseudoreaction name
    rxnName = compType['pseudoreaction']
    rxnName = rxnName.replace('P', 'protein')
    rxnName = rxnName.replace('C', 'carbohydrate')
    rxnName = rxnName.replace('N', 'biomass')
    rxnName = rxnName.replace('L', 'lipid backbone')
    rxnName = rxnName.replace('R', 'RNA')
    rxnName = rxnName.replace('D', 'DNA')
    rxnName = rxnName.replace('I', 'ion')
    rxnName = rxnName.replace('F', 'cofactor')

    # add up fraction
    rxnPos = model.rxnNames == rxnName
    if rxnPos.all() != 0:
        sub = model['S']
        isSub = sub[:, rxnPos](sub[:, rxnPos] < 0)
        if compType == 'L':
            F = - sum(sub[isSub, rxnPoS])
        else:
            F = 0
            # add up all components:
            for i in model['mets']:
                mets = model['mets']
                if data[mets] == mets[i]:
                    pos = 1
                if isSub[i] and sum(pos) == 1:
                    if compType == 'I' or compType == 'F':
                        MW = data['MWs'][pos]
                    else:
                        MW = data['MWs'][pos] - 18
                    abundance = -sub[i, rxnPos] * MW / 1000
                    F = F + abundance
        X = X + F

        print(str(compType) + ' -> ' + str(F) + " g/gDW")
    else:
        print(str(compType) + " do not exist")
        F = 0
        X = X + F
    return X, F