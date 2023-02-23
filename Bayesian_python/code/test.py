import multiprocessing
import numpy as np
import scipy.stats as ss
from multiprocessing import Process,cpu_count,Manager
import pickle
import logging
import module
from typing import Callable, Dict, Iterable, List
import numpy.typing as npt
import pandas as pd
import json
import cobra
import math
import re
import random
import statistics
import os
import shutil
from cobra.core import Reaction
from cobra.io.dict import model_to_dict
from cobra.util.solver import set_objective
from xml.dom import minidom
from optlang.symbolics import Zero, add
import scipy.io as scio
import abc_python_max
import addCarbonNum
import anaerobicModel
import changeGAM
import changeMedia
import getcarbonnum
import getComponentIndexes
import getFraction
import getrSample
import rmsecal
#import solveModel
import sumBioMass
import updateprior

generation = 150
enzymedataFile = "F:/2022/SJTU/DLKcat/Saccharomyces_cerevisiae_dl.mat"
z = scio.loadmat(enzymedataFile)

enzymedata = z['enzymedata']
max_growth = z['max_growth']
growthdata = z['growthdata']
model = z['model']
strain = z['strain']
rxn2block = z['rxn2block']

#tot_prot_weight = sumBioMass(model)  # [~,tot_prot_weight,~,~,~,~,~,~] = sumBioMass(model);
#if strain == 'Kluyveromyces_marxianus':
#    tot_prot_weight = 0.325
#elif strain == 'Kluyveromyces_lactis':
#    tot_prot_weight = 0.245
#else:
#    tot_prot_weight = tot_prot_weight*0.5

tot_prot_weight = 0.5
proc = 18
numPerGeneration = 126
rejectnum = 0.2
generation = 100
max_growth = [strain,'D-glucose',0.2,0,0,0,0,0,0,0,0,0,0,'aerobic','Batch','MIN'];

if len(max_growth) != 0 and len(growthdata) != 0:  # isempty
    max_growth = [strain, 'D-glucose', 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'aerobic', 'Batch', 'MIN']
    simulated3 = abc_python_max(model, enzymedata, enzymedata['kcat'], tot_prot_weight, growthdata, max_growth, 1, 1, 1, rxn2block)[2]
    if simulated3[0, 0] > 0.2:
        max_growth = [strain, 'D-glucose', simulated3[0, 0], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'aerobic', 'Batch', 'MIN']

D = abc_python_max(model, enzymedata, enzymedata['kcat'], tot_prot_weight, growthdata, max_growth, 1, 1, 1, rxn2block)
D_100 = D
theta_100 = []
kcat_100 = []
sampledgeneration = 1
while D > rejectnum and D_100 > 0.5:

    nfound = len(os.listdir(r'F:\2022\SJTU\DLKcat\Saccharomyces_cerevisiae2'))
    if nfound > 0:
        tmp = np.loadtxt('kcat_genra' + str(nfound + 1) + '.txt')
        theta_100 = tmp[-1, :]
        kcat_100 = tmp[0:-3, :]
        tot_prot_weight = tmp[-2, 0]
        sampledgeneration = nfound + 1
        ss = np.transpose(kcat_100)  ## ss = num2cell(kcat_100',1);
        [a, b] = updateprior(ss)
        enzymedata['kcat'] = np.transpose(a)
        enzymedata['kcat_var'] = np.transpose(b)

    if sampledgeneration <= generation:
        print("No." + str(sampledgeneration) + 'generation')
        # generate a
        old = theta_100
        kcat_old_100 = kcat_100
        # repeat a generation

        if sampledgeneration == 1:
            sample_generation = 144
        else:
            sample_generation = numPerGeneration

        # generate one generation sample of kcats
        kcat_random_all = getrSample(enzymedata['kcat'], enzymedata['kcat_var'],
                                     enzymedata['enzyme_ec_kcat_range'][:, 0],
                                     enzymedata['enzyme_ec_kcat_range'][:, 1],
                                     np.matlib.repmat(sample_generation, len(enzymedata.kcat), 1),
                                     'Uniform')
        kcat_random_all = np.array(kcat_random_all)
        print("kcat random finish")

        # start sampling
        new_tmp = []
        for i in range(proc):
            i
            rmse_final = abc_python_max(model, enzymedata, kcat_random_all, tot_prot_weight, growthdata, max_growth,
                                        proc, sample_generation, i, rxn2block)
            new_tmp[i] = rmse_final

        print("RMSE calculation finish")

        # find the best 100 samples
        new = np.array(new_tmp)
        theta = [new, old]
        kcat = [kcat_random_all, kcat_old_100]

        # initialize an empty set to store the best 100 after each step
        D_idx = theta.sort()
        theta_100 = theta[D_idx[0:99]]
        D = abs(theta_100[99] - theta_100[0])
        D_100 = theta_100[99]
        kcat_100 = kcat[:, D_idx[0:100]]
        tot = np.matlib.repmat(tot_prot_weight, 1, len(theta_100))
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

