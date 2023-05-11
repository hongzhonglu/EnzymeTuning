import numpy as np
import os
import cobra
from pathlib import Path
from sklearn.metrics import mean_squared_error
import pandas as pd
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import scipy.io as scio
import importlib
from scipy import sparse
import sys
sys.path.append(r'F:\python\Bayesian_python\code')
from bayesian_generation import abc_python_max
#import anaerobicModel
#import changeGAM
#import changeMedia
from bayesian_generation import getrSample
from bayesian_generation import sumBiomass
from bayesian_generation import updateprior

species = "Saccharomyces_cerevisiae"
generation = 150
enzymedataFile = "data/Saccharomyces_cerevisiae_dl.mat"
#enzymedataFile = "data/model.mat"
z = scio.loadmat(enzymedataFile)
enzymedata = z['enzymedata'][0,0]
max_growth = z['max_growth']
growthdata = z['growthdata']
model = z['model'][0,0]
strain = z['strain']
rxn2block = z['rxn2block']

#split the training data
growthdata_train = growthdata[1:17:2,:]
growthdata_test = growthdata[0:17:2,:]
max_growth_train = max_growth[1:8:2,:]
max_growth_test = max_growth[0:8:2,:]

import scipy.io as io
F = sumBiomass(model)
#calculate the prot_weight
#F = 0.46  #怎么从sumBiomass中导出结果，matlab版本为0.46
tot_prot_weight = F[0]*0.5
if strain == 'Kluyveromyces_marxianus':
    tot_prot_weight = 0.325
elif strain == 'Kluyveromyces_lactis':
    tot_prot_weight = 0.245

proc = 18
numPerGeneration = 126 #126/18 = 7
rejectnum = 0.2
generation = 100

if len(max_growth) == 0 and len(growthdata) == 0:
    max_growth = [strain,'D-glucose',0.2,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,'aerobic','Batch','MIN']
    simulated_3 = abc_python_max(model,enzymedata,enzymedata['kcat'],tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block)
    simulated_3_tmp = np.array(np.array(simulated_3)[2][0])[0,1]
    if simulated_3_tmp >0.2:
        max_growth = [strain,'D-glucose',simulated_3_tmp,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,'aerobic','Batch','MIN']

D = abc_python_max(model,enzymedata,enzymedata["kcat"],tot_prot_weight,growthdata,max_growth,1,1,1,rxn2block)
D_100 = D[0]
theta_100 = []
kcat_100 = []
sampledgeneration = 1

#j = 1
#while j <= 1:
while D[0] > rejectnum and D_100 > 0.5:
    nfound = len(os.listdir(r'F:\2022\SJTU\DLKcat\Saccharomyces_cerevisiae2'))
    if nfound > 0:
        tmp = pd.read_csv('F:/2022/SJTU/DLKcat/Saccharomyces_cerevisiae2/'+'kcat_genra' + str(nfound) + '.txt', delimiter=',',header= None)
        tmp = np.array(tmp)
        theta_100 = tmp[-1, :]
        kcat_100 = tmp[0:-3, :]
        tot_prot_weight = tmp[-2, 0]
        sampledgeneration = nfound + 1
        # recalculate the sigma and mu
        [a, b] = updateprior(kcat_100)
        enzymedata['kcat'] = np.transpose(a)
        enzymedata['kcat_var'] = np.transpose(b)

    if sampledgeneration <= generation:
        print("No. " + str(sampledgeneration) + ' generation')

        # generate a
        old = theta_100
        kcat_old_100 = kcat_100

        # repeat a generation
        if sampledgeneration == 1:
            sample_generation = 144
        else:
            sample_generation = numPerGeneration

        # generate one generation sample of kcats
        kcat_random_all = getrSample(enzymedata['kcat'], enzymedata['kcat_var'], enzymedata['enzyme_ec_kcat_range'][:, 0],
                                     enzymedata['enzyme_ec_kcat_range'][:, 1], np.tile(sample_generation, (1, 1)),
                                     method = 'Uniform')
        print("kcat random finish !")
        #print(kcat_random_all)
        # start sampling
        new_tmp = []
        for i in range(proc):
            print(i+1)
            rmse_final = abc_python_max(model, enzymedata, kcat_random_all, tot_prot_weight, growthdata, max_growth,proc, sample_generation, i+1, rxn2block)
            new_tmp.append(rmse_final[0])

        print("RMSE calculation finish !")

        # find the best 100 samples
        new = new_tmp
        theta = new + old
        kcat = kcat_random_all + kcat_old_100

        # initialize an empty set to store the best 100 after each step
        D_idx = theta.sort()
        theta_100 = theta[D_idx[0:100]]
        D = abs(theta_100[99] - theta_100[0])
        D_100 = theta_100[99]
        kcat_100 = kcat[:, D_idx[0:100]]
        tot = np.tile(tot_prot_weight, (1, len(theta_100)))
        kcat_genra = []
        for m in range(len(kcat_100)):
            for n in range(len(tot)):
                if m == n:
                    for k in range(len(theta_100)):
                        if n == k:
                            t = (kcat_100[m], tot[n], theta_100[k])
                            kcat_genra.append(t)
        np.savetxt('kcat_genra' + str(sampledgeneration) + '.txt', kcat_genra)

        # recalclate the sigma and mu
        sss = np.transpose(kcat_100[:,0])
        ss = [[x] for x in sss]

        [a, b] = updateprior(ss)
        enzymedata['kcat'] = np.transpose(a)
        enzymedata['kcat_var'] = np.transpose(b)
        sampledgeneration = sampledgeneration + 1

    else:
        D = rejectnum
        D_100 = D