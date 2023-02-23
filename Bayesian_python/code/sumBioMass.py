##sumbiomass
#   model    (struct) Metabolic model in COBRA format
#
#   X         (float) Total biomass fraction [gDW/gDW]
#   P         (float) Protein fraction [g/gDW]
#   C         (float) Carbohydrate fraction [g/gDW]
#   R         (float) RNA fraction [g/gDW]
#   D         (float) DNA fraction [g/gDW]
#   L         (float) Lipid fraction [g/gDW]
#   F         (float) cofactor [g/gDW]
#   I         (float) ion [g/gDW]
def sumBioMass(model):
    biomassdatafile = "F:/2022/SJTU/DLKcat/github/DLKcat-master/BayesianApproach/Code/common/biomassCompData.mat"
    data = scio.loadmat(biomassdatafile)

    [P, X] = getFraction(model, data, 'P', 0)
    [C, X] = getFraction(model, data, 'C', X)
    [R, X] = getFraction(model, data, 'R', X)
    [D, X] = getFraction(model, data, 'D', X)
    [L, X] = getFraction(model, data, 'L', X)
    [I, X] = getFraction(model, data, 'I', X)
    [F, X] = getFraction(model, data, 'F', X)

    print("X ->" + str(X) + 'gDW/DW')
    biomass_file = np.array([[P, X], [C, X], [R, X], [D, X], [L, X], [I, X], [F, X]])
    # simulate growth
    sol = model.optimize()
    print('Growth = ' + str(sol) + ' 1/h')
    return biomass_file