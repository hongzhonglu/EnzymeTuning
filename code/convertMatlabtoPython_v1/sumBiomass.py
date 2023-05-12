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

