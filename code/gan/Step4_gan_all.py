from src.ec_N_lim_sen import *
from src.constrain_ecmodel import *
from src.mainFunction import *
from src.protein_process import *
from src.model_process import *
from src.code_generation_preparation import *
from multiprocessing import Pool , cpu_count
import numpy as np
import pandas as pd
import h5py
from io import StringIO
import os
import pickle
from src.GAN import CGAN
import helper as hp

def main():
    # load the enzyme_id
    #enzyme_preparation = pd.read_csv (f'gan_input/kcat_merge_sen_50.csv')
    enzyme_preparation = pd.read_csv('data/kcat_all.csv')
    enzyme_preparation_all = pd.read_csv('data/kcat_all.csv')
    #rxn_name = enzyme_preparation.drop (['kcat_value1'] , axis=1 , inplace=False)

    # parameter determination
    number_of_models = 1000  # No of kcat samples


    # Training hyperparameters
    latent_dim = 127      # Length of noise vector
    epochs = 50           # Total number of epochs
    n_sample = 1000        # No of parameter sets to generate at every sampling epoch
    repeats = 1           # number of training repeats
    batchsize = 50      # Batchsize
    sample_interval = 10  # Frequency of testing generator

    # load enzyme information
    enzymedataFile = "data/Saccharomyces_cerevisiae_dl.mat"
    z = scio.loadmat (enzymedataFile)
    model = z['model'][0 , 0]
    enzymedata = z['enzymedata'][0 , 0]
    rxn2block = z['rxn2block']
    max_growth_train = z['max_growth']
    growthdata_train = z['growthdata']
    model = z['model'][0, 0]
    strain = z['strain']
    rxnlist = enzyme_preparation_all.drop (['kcat_value1'], axis=1, inplace=False)
    savepath = 'data'
    file_path = os.path.join (savepath , 'rxnlist.pkl')
    with open (file_path , 'wb') as file:
        pickle.dump (rxnlist , file)
    print(rxnlist)

    #construct kcat_list
    kcat_initial= pd.read_csv("data/kcat_genra1.txt",delimiter=",",header=None)
    kcat_initial = np.array(kcat_initial)
    kcat_100 = kcat_initial[0:-2,:]
    [a,b] = updateprior(kcat_100)
    enzymedata['kcat'] = np.transpose(a)
    enzymedata['kcat_var'] = np.transpose(b)
    MWs = construct_MWs (enzymedata)
    kcat_dict = construct_kcat_dict (model , enzymedata)
    kcat_df = pd.DataFrame (list (kcat_dict.items ()) , columns=['gene_rxnID' , 'kcat_value'])
    kcat_df[['gene' , 'rxnID']] = pd.DataFrame (kcat_df['gene_rxnID'].tolist () , index=kcat_df.index)
    merged_df = pd.merge (enzyme_preparation , kcat_df , on=['gene' , 'rxnID'])
    kcat_final = merged_df[['gene' , 'rxnID' , 'kcat_value']]
    kcat_all = kcat_df.drop (['gene_rxnID'], axis=1, inplace=False)

    # load ecModel
    model_cobra = load_matlab_model ('data/model/model.mat')
    eModel = convertToEnzymeModel (model_cobra , kcat_dict)
    eModel = constrainPool (eModel , MWs , {} , eModel.enzymes , 230)


    # proteomics data (Clim & Nlim)
    eModel_genes = [gene.id for gene in eModel.genes]
    eModel_genes_series = pd.DataFrame (eModel_genes , columns=['emodel_genes'])
    proteomics_Clim = pd.read_csv ('data/Clim_all.csv')
    proteomicsdata_Clim = proteomics_Clim[proteomics_Clim["all_gene"].isin (eModel_genes_series["emodel_genes"])]
    proteomics_Nlim = pd.read_csv ('data/Nlim_train.csv')
    proteomicsdata_Nlim = proteomics_Nlim[proteomics_Nlim["all_gene"].isin (eModel_genes_series["emodel_genes"])]

    # growth data
    growthdata = pd.read_excel ('data/growthdata.xlsx' , sheet_name="Sheet1")
    growthdata = np.array (growthdata)  # 7N+9C

    # update the kcat
    ex_mets = ['biomass pseudoreaction' , 'D-glucose exchange' , 'acetate exchange' , 'ethanol exchange' ,
               'glycerol exchange' , 'pyruvate exchange' , 'ethyl acetate exchange' , 'carbon dioxide exchange' ,
               'oxygen exchange' , 'prot_pool_exchange']

    # find the related rxnID
    idx = []
    for name0 in ex_mets:
        s = getRxnByReactionName(model=eModel, name=name0)
        if len(s) > 1:
            print("need check")
        elif len(s) == 1:
            idx.append(s[0])

    # fetch training set eigenvalues
    path_stability = 'data/gan_input_all/pro_data_1000_both.csv'
    if not path_stability.endswith ('.csv'):
        raise ValueError ('Your data must be a .csv file')
    pro_data = pd.read_csv (path_stability)

    # load initial 1000 kcat samples

    path_kcat = 'gan_input/kcat_data_new_1000_50_score.h5'
    if "RMSE" in pro_data :
        path_kcat = 'data/gan_input_all/kcat_data_both_1000.h5'
    f = h5py.File (path_kcat , 'r')

    # load rxnlist(kcat_name)
    path_parameters = 'gan_input/parameter_sen_50.pkl'
    if "RMSE" in pro_data :
        path_parameters = 'data/gan_input_all/parameter.pkl'
    with open ( path_parameters, 'rb') as file:
        parameter_names = pickle.load (file)

    w1, w2 = 1, 1  # weight
    if "error" in pro_data :
        # rxnlist = rxn_name
        train_type = "protein"
        parameter_set_dim = len (enzyme_preparation)  # No of kcat parameters in model
        error = pro_data['error'].values
        corr = pro_data['corr'].values
        error_nor = (error - error.min ())/(error.max () - error.min ())
        corr_nor = 1 - (corr - corr.min ())/(corr.max () - corr.min ())
        stabilities_0 = error_nor + corr_nor
        if "RMSE" in pro_data :
            train_type = "protein_phenotype"
            parameter_set_dim = len (enzyme_preparation_all)  # No of kcat parameters in model
            RMSE = pro_data['RMSE'].values
            RMSE_nor = (RMSE - RMSE.min ())/(RMSE.max () - RMSE.min ())
            stabilities_0 = w1*(error_nor + corr_nor) + w2*RMSE_nor
    else :
        train_type = "phenotype"
        parameter_set_dim = len (enzyme_preparation_all)  # No of kcat parameters in model
        RMSE = pro_data['RMSE'].values
        stabilities_0 = RMSE

    start = time.time ()
    print(train_type)
    print(parameter_set_dim)
    print ('\nSTARTING PREPROCESSING')
    for j in range (3) : #10, 30, 50
        all_data = np.empty ([number_of_models, parameter_set_dim])
        all_stabilities = np.empty ([number_of_models])

        J_partition = np.median (stabilities_0)  # <--- Create class partition based on this eigenvalue
        count0, count1 = 0, 0

        for i in range (0, number_of_models) :

            if i%100==0 :
                print (f'current set processed: {i}')
            this_param_set = f'kcat{i}'
            param_values = np.array (f.get (this_param_set))

            mreal = stabilities_0[i]

            if mreal <= J_partition :
                stability = 1
                count0 += 1
            elif mreal > J_partition :
                stability = -1
                count1 += 1

            all_data[i] = param_values
            all_stabilities[i] = stability

        all_data = np.array (all_data)
        print (all_data)
        all_stabilities = np.array (all_stabilities)

        n_param = all_data.shape[0]
        print (f'% relevant models: {count1/n_param}')

        # take the log
        log_all_data = np.log (all_data)  # Log transform all parameters

        # train-val split
        ratio = float (0.9)  # Partition of training and test data
        n_data = log_all_data.shape[0]
        limit = int (ratio*n_data)
        all_idx = np.arange (n_data)
        np.random.shuffle (all_idx)

        idx_tr = all_idx[:limit]
        idx_val = all_idx[limit :]

        tr_data = log_all_data[idx_tr]
        val_data = log_all_data[idx_val]

        tr_stabi = all_stabilities[idx_tr]
        val_stabi = all_stabilities[idx_val]

        print (f'N data for training: {tr_data.shape[0]}')
        print (f'N data for validation: {val_data.shape[0]}')

        # save everything
        exp_id = 'test' + str (j)
        savepath = f'result/gan_all/gan_input/{exp_id}'
        os.makedirs (savepath, exist_ok=True)
        np.save (f'{savepath}/all_kcat_{exp_id}.npy', all_data)
        np.save (f'{savepath}/all_targets_{exp_id}.npy', all_stabilities)
        np.save (f'{savepath}/X_train_{exp_id}.npy', tr_data)
        np.save (f'{savepath}/X_val_{exp_id}.npy', val_data)

        np.save (f'{savepath}/y_train_{exp_id}.npy', tr_stabi)
        np.save (f'{savepath}/y_val_{exp_id}.npy', val_stabi)

        with open (f'{savepath}/parameter_names_{exp_id}.pkl', 'wb') as f :
            pickle.dump (parameter_names, f)

        path_generator = None  # <---if doing transfer learning put path to trained generator here else leave None
        #    if loading model using load_model gives an error upgrade tensorflow to v2.3.0
        #    > pip install tensorflow==2.3.0

        # load the data for appropriate experiment
        X_train = tr_data
        y_train = tr_stabi

        # Specify output folders
        savepath = f'result/gan_all/gan_output/{exp_id}/'
        os.makedirs (savepath, exist_ok=True)
        for k in range (0, repeats) :
            print (f'Current exp: {exp_id}: Samples used: {np.shape (X_train)[0]}, repeat {k}')
            # set save directory
            this_savepath = f'{savepath}repeat_{k}/'
            os.makedirs (this_savepath, exist_ok=True)

            cgan = CGAN (X_train, y_train, latent_dim, batchsize, path_generator, savepath=this_savepath)
            d_loss, g_loss, acc = cgan.train (epochs, sample_interval, n_sample)

            # store training summary

            this_train_savepath = f'{this_savepath}training_summary/'
            os.makedirs (this_train_savepath, exist_ok=True)

            with open (f'{this_train_savepath}d_loss.pkl', 'wb') as f :
                pickle.dump (d_loss, f)
            with open (f'{this_train_savepath}g_loss.pkl', 'wb') as f :
                pickle.dump (g_loss, f)
            with open (f'{this_train_savepath}acc.pkl', 'wb') as f :
                pickle.dump (acc, f)

        # plot metrics
        x_plot = np.arange (0, epochs, 1)
        fig = plt.figure (figsize=(10, 5))
        ax1 = fig.add_axes ([0.2, 0.2, 1, 1])
        ax1.plot (x_plot, d_loss, label='discriminator loss')
        ax1.plot (x_plot, g_loss, label='generator loss')
        ax1.set (ylabel='criterion_losses', xlabel='epochs')
        ax1.legend ()
        plt.savefig (f'{this_train_savepath}loss.svg', dpi=300,
                     transparent=False, bbox_inches='tight')

        fig = plt.figure (figsize=(10, 5))
        ax2 = fig.add_axes ([0.2, 0.2, 1, 1])
        ax2.plot (x_plot, d_loss, label='discriminator accuracy')
        ax2.set (ylabel='accuracy', xlabel='epochs')
        ax2.legend ()
        plt.savefig (f'{this_train_savepath}d_accuracy.svg', dpi=300,
                     transparent=False, bbox_inches='tight')

        # load new kcat
        kcat_new_name = f'result/gan_all/gan_output/test{j}/repeat_{repeats - 1}/490_r.npy'
        kcat_new = np.load (kcat_new_name)
        kcat_new = np.exp (kcat_new)
        kcat_new = np.transpose (kcat_new)
        kcat_df = pd.DataFrame (kcat_new)
        column_names = [f'kcat_value{i}' for i in range (kcat_df.shape[1])]
        kcat_df.columns = column_names
        kcat_final_tmp = pd.concat ([rxnlist, kcat_df], axis=1)
        print (kcat_final_tmp)

        # calculate new stabilities
        num_cpus = 8
        print (num_cpus)
        with Pool (processes=num_cpus) as pool :
            args = [(
                j, growthdata,growthdata_train, max_growth_train, eModel, model, kcat_final_tmp, proteomicsdata_Clim, proteomicsdata_Nlim, train_type, rxn2block) for
                j in range (number_of_models)]
            result = pool.starmap (parallel_task_preprocessing_all, args)
            print (result)
        # pool.close ()
        # pool.join ()

        results = np.vstack (result)
        toy_data = pd.DataFrame (results, columns=['error', 'corr', 'pvalue', 'number','RMSE'])

        # save the result and sort to new inputs
        RMSE_new = toy_data['RMSE']
        error_new = (toy_data['error'] - toy_data['error'].min ())/(toy_data['error'].max () - toy_data['error'].min ())
        corr_new = 1 - (toy_data['corr'] - toy_data['corr'].min ())/(toy_data['corr'].max () - toy_data['corr'].min ())
        if "phenotype" in train_type:
            RMSE_new = (toy_data['RMSE'] - toy_data['RMSE'].min ())/(toy_data['RMSE'].max () - toy_data['RMSE'].min ())
        w1, w2 = 1, 1  # weight
        stability_new = w1*(error_new + corr_new) + w2*RMSE_new
        if j > 0 :
            stability_old = pd.read_csv (file_path_sta).iloc[:, 0].values
        else :
            stability_old = pd.read_csv (path_stability).iloc[:, 0].values

        stability_all = np.append (stability_old, stability_new)
        stability_sort = sorted (range (len (stability_all)), key=lambda k : stability_all[k], reverse=False)
        stabilities_0 = stability_all[stability_sort[0 :1000]]
        savepath_sta = f'result/gan_all/gan_output/test{j}'
        os.makedirs (savepath_sta, exist_ok=True)
        file_path_sta = os.path.join (savepath_sta, f'pro_data_{j}.csv')
        file_path_toy = os.path.join (savepath_sta , f'result_{j}.csv')
        np.savetxt (file_path_sta, stabilities_0)
        np.savetxt (file_path_toy , toy_data)

        kcat_old = all_data
        kcat_new = np.transpose (kcat_new)
        kcat_total = np.concatenate ((kcat_old, kcat_new))
        kcat_next = kcat_total[stability_sort[0 :1000], :]

        path_parameters = f'result/gan_all/gan_input/test{j}/kcat_data_1000.h5'
        with h5py.File (path_parameters, 'w') as h5f :
            for q in range (len (kcat_next)) :
                # 获取第 q 个元素
                data_to_write = kcat_next[q, :]
                h5f.create_dataset (f'kcat{q}', data=data_to_write)

        f = h5py.File (path_parameters, 'r')
    end = time.time ()
    print (f'PROCESSING DONE in {end - start:.05} seconds')


if __name__=='__main__' :
    main ()