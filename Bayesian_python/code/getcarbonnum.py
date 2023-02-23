def getcarbonnum(model, exrxn):
    # exrxn must be an exchange rxn and only contain one mets in the equation
    idx = model['rxns'].find(exrxn)  ##ismember
    EXmets = model['S'][:, idx]

    # find all exchange mets
    EXmetsIdx = np.zeros((len(exrxn), 1))
    for k in range(len(EXmets[0, :])):
        EXmetIdx[k] = np.nonzero(EXmets[:, k])

    EXfors = model['metFormulas'][EXmetsIdx]
    # Ematrix = getElementalComposition(EXfors,{'C'})[0]  #cobratoolbox
    # elements = getElementalComposition(EXfors,{'C'})[1]
    Ematrix = Ematrix[:, 0]
    changedmets = exrxn[np.where(np.isnan(Ematrix))]
    Ematrix[np.where(np.isnan(Ematrix))] = 1
    CarbonNum = np.transpose(Ematrix)
