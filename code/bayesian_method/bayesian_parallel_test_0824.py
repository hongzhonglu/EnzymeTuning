def main():
    import numpy as np
    import os
    import cobra
    from pathlib import Path
    from sklearn.metrics import mean_squared_error
    import pandas as pd
    from cobra.io import load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
    import scipy.io as scio
    from scipy import sparse
    import cobra.util
    from cobra import Model,Reaction,Metabolite
    from multiprocessing import Pool,cpu_count

    from code0824 import abc_python_max
    from code0824 import parallel_task
    #from Bayesian_generation_cobra import anaerobicModel
    #from Bayesian_generation_cobra import changeGAM
    #from Bayesian_generation_cobra import changeMedia
    from code_generation_preparation import getrSample
    from code_generation_preparation import sumBiomass
    from code_generation_preparation import updateprior

    path = "F:\python\Bayesian_python"
    os.chdir(path)

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

    #load the model in cobra format
    scio.savemat('data/model.mat', {'model': model})
    model_cobra = load_matlab_model('data/model.mat')

    import scipy.io as io
    F = sumBiomass(model,model_cobra)
    #calculate the prot_weight
    #F = 0.46  #怎么从sumBiomass中导出结果，matlab版本为0.46
    tot_prot_weight = F[0]*0.5
    if strain == 'Kluyveromyces_marxianus':
        tot_prot_weight = 0.325
    elif strain == 'Kluyveromyces_lactis':
        tot_prot_weight = 0.245

    met = cobra.Metabolite ('cost' , name='cost' , compartment='c')
    model_cobra.add_metabolites (met)
    model_cobra.add_boundary (model_cobra.metabolites.get_by_id ("cost") , type="sink" , ub=tot_prot_weight * 1000 , lb=0)

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

    output = 0


    D = abc_python_max (model , model_cobra , enzymedata , enzymedata["kcat"] , tot_prot_weight , growthdata ,
                        max_growth , 1 , 1 , 1 , rxn2block,output)[0]
    D = D[0][0]
    D_100 = D
    theta_100 = []
    kcat_100 = []
    sampledgeneration = 1

    #j = 1

    #while j <= 1:
    while D > rejectnum and D_100 > 0.5:
        output = output + 1
        nfound = len(os.listdir(r'F:/python/Bayesian_python/result'))
        print(nfound)
        if nfound > 0:
            tmp = pd.read_csv('result/'+'kcat_genra' + str(output) + '.txt', delimiter=',',header= None)
            tmp = np.array(tmp)
            theta_100 = tmp[-1, :]
            kcat_100 = tmp[0:-2, :]
            tot_prot_weight = tmp[-2, 0]
            sampledgeneration = nfound + 1
            # recalculate the sigma and mu
            [a, b] = updateprior(kcat_100)
            enzymedata['kcat'] = np.transpose(a)
            enzymedata['kcat_var'] = np.transpose(b)
            print(len(enzymedata['kcat']))
        if sampledgeneration <= generation:
            print("No. " + str(output) + ' generation')

            # generate a
            old = theta_100
            kcat_old_100 = kcat_100

            # repeat a generation
            if sampledgeneration == 1:
                sample_generation = 144
            else:
                sample_generation = numPerGeneration

            # generate one generation sample of kcats
            #kcat_random_all = np.concatenate((kcat_100,kcat_100[:,0:26]),axis=1)
            kcat_random_all = getrSample(enzymedata['kcat'], enzymedata['kcat_var'], enzymedata['enzyme_ec_kcat_range'][:, 0],enzymedata['enzyme_ec_kcat_range'][:, 1], sample_generation,method = 'Uniform')
            print("kcat random finish !")
            #print(kcat_random_all)
            # start sampling
            new_tmp = np.zeros((18, 7))



            # Use a Pool to parallelize the work among 20 cores
            num_cpus = min (cpu_count () , 18)
            print(num_cpus)
            with Pool (processes=num_cpus) as pool:
                args = [(
                        i , output , model , model_cobra , enzymedata , kcat_random_all , tot_prot_weight , growthdata ,
                        max_growth , proc , sample_generation , rxn2block) for i in range (int (proc))]
                results = pool.starmap (parallel_task , args)

            new_tmp = np.vstack (results)
            print(new_tmp)
            print ("RMSE calculation finish !")

            # ... [继续其它代码]

            #for i in range(int(proc)):
            #    print(i+1)
                # Call start_memory_trace() before the code you want to profile
               # start_memory_trace ()

            #    rmse_final = abc_python_max(model, model_cobra, enzymedata, kcat_random_all, tot_prot_weight, growthdata, max_growth,int(proc), int(sample_generation), i+1, rxn2block,output)
            #    print(rmse_final)
            #    new_tmp[i,:] = rmse_final[0][0]
            #    #new_tmp = np.transpose(new_tmp)
            #    path = 'result/output_rmse'
            #    os.makedirs(path,exist_ok=True)
            #    filename = 'output_'+str(output)+'_rmse'+str(i+1)+'.txt'
            #    full_path = os.path.join(path,filename)
            #    np.savetxt (full_path , new_tmp)

                # Call stop_memory_trace() after the code you want to profile
              #  stop_memory_trace ()
            #print("RMSE calculation finish !")

            # find the best 100 samples
            new = new_tmp
            theta = np.append (old , new)
            kcat = np.append (kcat_old_100 , kcat_random_all , axis=1)

            # initialize an empty set to store the best 100 after each step
            theta_100 = np.sort (theta.flatten ())[:100]
            D = abs (theta_100[99] - theta_100[0])
            D_idx = np.argsort (theta.flatten ())[:100]
            D_100 = theta_100[99]
            kcat_100 = kcat[: , D_idx]
            tot = np.tile(tot_prot_weight, (1, len(theta_100)))
            kcat_genra = []
            for m in range(len(kcat_100[0,:])):
                for n in range(len(tot[0,:])):
                    if m == n:
                        for k in range(len(theta_100)):
                            if n == k:
                                t = np.append(kcat_100[:,m], tot[:,n])
                                t = np.append(t, theta_100[k])
                                kcat_genra.append(t)
            kcat_genra = np.transpose (kcat_genra)
            path = 'result'
            os.makedirs (path , exist_ok=True)
            filename = 'kcat_genra' + str (output+1) + '.txt'
            full_path = os.path.join (path , filename)
            np.savetxt (full_path , kcat_genra,delimiter=",", fmt="%f")

            # recalclate the sigma and mu
            sss = np.transpose(kcat_100[:,0])
            ss = [[x] for x in sss]

            [a, b] = updateprior(ss)
            enzymedata['kcat'] = np.transpose(a)
            enzymedata['kcat_var'] = np.transpose(b)
            sampledgeneration = sampledgeneration + 1
            print(D)
            print(D_100)


        else:
            D = rejectnum
            D_100 = D

if __name__ == '__main__':
    main()
