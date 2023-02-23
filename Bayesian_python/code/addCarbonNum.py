def addCarbonNum(model):
    model["excarbon"] = np.zeros((1,len(model['rxns'])))
    #EXrxn = getExchangeRxns(model)[0]
    #EXrxnIdx = getExchangeRxns(model)[1]
    CarbonNum = getcarbonnum(model,EXrxn)[0]
    EXfors = getcarbonnum(model,EXrxn)[1]
    model['excarbon']['EXrxnIdx'] = CarbonNum
    for i in range(len(model['rxnNames'])):
        if model['rxnNames'] == 'growth':
            model['excarbon'] = 41