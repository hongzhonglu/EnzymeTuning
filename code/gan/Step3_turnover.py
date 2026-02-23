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
from src.ec_N_lim_sen import updateEcGEMkcat

def main():
    miu = 0.424
    N = 100
    turnover_data = pd.read_csv (f'data/turnover/test/top_{N}_{miu}.csv')
    number_of_models = 1000  # No of vsyn samples
    parameter_set_dim = len (turnover_data)  # No of vsyn parameters in model

    # Training hyperparameters
    latent_dim = 127      # Length of noise vector 127
    epochs = 500           # Total number of epochs
    n_sample = 1000        # No of parameter sets to generate at every sampling epoch
    repeats = 1           # number of training repeats 1
    batchsize = 50      # Batchsize 50
    sample_interval = 10  # Frequency of testing generator

    # load the dl_model
    enzymedataFile = "data/model/Saccharomyces_cerevisiae_dl.mat"
    z = scio.loadmat (enzymedataFile)
    model = z['model'][0, 0]
    enzymedata = z['enzymedata'][0, 0]
    # scio.savemat('data/model/model.mat', {'model': model})
    model_cobra = load_matlab_model ('data/model/model.mat')

    # construnct ecModel
    # kcat_tmp = pd.read_csv("data/kcat_genra4.txt",delimiter=",",header=None)
    enzymedata = z['enzymedata'][0, 0]
    MWs = construct_MWs (enzymedata)
    kcat_tmp = pd.read_csv ("data/kcat_genra100.txt", delimiter=",", header=None)
    kcat_tmp = np.array (kcat_tmp)
    kcat_100 = kcat_tmp[0 :-2, :]
    [a, b] = updateprior (kcat_100)
    enzymedata['kcat'] = np.transpose (a)
    enzymedata['kcat_var'] = np.transpose (b)
    MWs = construct_MWs (enzymedata)
    kcat_dict = construct_kcat_dict (model, enzymedata)
    eModel_bayesian = convertToEnzymeModel (model_cobra, kcat_dict)
    eModel_bayesian = constrainPool (eModel_bayesian, MWs, {}, eModel_bayesian.enzymes, 230)
    # s2 = get_r_max (eModel_bayesian , model_cobra)

    enzyme_preparation = pd.read_csv ('data/kcat_merge_sen_50.csv')
    # kcat_new = pd.read_csv ('/merged_data.csv')
    path_parameters = 'data/turnover/kcat_data_1000.h5'
    f = h5py.File (path_parameters, 'r')
    kcat_new = np.empty ([1000, len (enzyme_preparation)])
    for i in range (0, 1000) :
        this_param_set = f'kcat{i}'
        param_values = np.array (f.get (this_param_set))
        kcat_new[i] = param_values
    kcat_new = np.array (kcat_new)
    kcat_new_value = np.mean (kcat_new, axis=0)
    kcat_df_new = pd.DataFrame (kcat_new_value)
    column_names = ['kcat_value0']
    kcat_df_new.columns = column_names
    rxn_name = enzyme_preparation.drop (['kcat_value1'], axis=1, inplace=False)
    kcat_final_tmp = pd.concat ([rxn_name, kcat_df_new], axis=1)
    index = 'kcat_value0'
    with eModel_bayesian as model_tmp :
        # model_tmp = constrainPool(model_tmp, MWs, {}, model_tmp.enzymes, 230)
        # s2 = get_r_max(model_tmp,model_cobra)
        for k in range (len (kcat_final_tmp)) :
            target_gene0 = kcat_final_tmp['gene'][k]
            kcat_m0 = kcat_final_tmp[index][k]
            rxn0 = kcat_final_tmp['rxnID'][k]
            model_tmp = updateEcGEMkcat (ecModel=model_tmp, target_gene=target_gene0, rxnID=rxn0, kcat_m=kcat_m0)
        eModel_gan = model_tmp.copy ()

    # proteomics data (Clim & Nlim)
    eModel_genes = [gene.id for gene in eModel_gan.genes]
    eModel_genes_series = pd.DataFrame (eModel_genes , columns=['emodel_genes'])
    proteomics_Clim = pd.read_csv ('data/Clim_all.csv')
    proteomicsdata_Clim = proteomics_Clim[proteomics_Clim["all_gene"].isin (eModel_genes_series["emodel_genes"])]
    # proteomics_Nlim = pd.read_csv ('data/Nlim_train.csv')
    proteomics_Nlim = pd.read_csv ('data/data_PNAS_2021.csv')
    proteomicsdata_Nlim = proteomics_Nlim[proteomics_Nlim["all_gene"].isin (eModel_genes_series["emodel_genes"])]

    # growth data
    growthdata = pd.read_excel ('data/growthdata.xlsx' , sheet_name="Sheet3")
    growthdata = np.array (growthdata)  # 7N+9C

    # update the kcat
    ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
               'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
               'oxygen exchange', 'prot_pool_exchange']
    # find the related rxnID
    idx = []
    for name0 in ex_mets :
        s = getRxnByReactionName (model=eModel_gan, name=name0)
        if len (s) > 1 :
            print ("need check")
        elif len (s)==1 :
            idx.append (s[0])

    path_parameters = f'data/turnover/test/{N}/vsyn_data_{miu}_1000.h5'
    f = h5py.File (path_parameters, 'r')

    with open (f'data/turnover/test/{N}/parameter.pkl' , 'rb') as file:
        parameter_names = pickle.load (file)

    path_stability = f'data/turnover/test/{N}/pro_1000_{miu}.csv'
    if not path_stability.endswith ('.csv') :
        raise ValueError ('Your data must be a .csv file')
    pro_data = pd.read_csv (path_stability).iloc[:, 0 :4].values
    error = (pro_data[:, 0] - pro_data[:, 0].min ())/(pro_data[:, 0].max () - pro_data[:, 0].min ())
    corr = 1 - (pro_data[:, 1] - pro_data[:, 1].min ())/(pro_data[:, 1].max () - pro_data[:, 1].min ())
    number = 1 - (pro_data[ : , 3 ] - pro_data[ : , 3 ].min ()) / (pro_data[ : , 3 ].max () - pro_data[ : , 3 ].min ())
    w1, w2, w3 = 0, 1, 0  # weight
    stabilities_0 = w1*error + w2*corr + w3*number
    # stabilities_0 = 1-pro_data[:, 1]

    start = time.time ()
    print ('\nSTARTING PREPROCESSING')

    best_gen_loss = float ('inf')
    best_gan_output_path = None
    for j in range (3) : #10, 30, 50
        all_data = np.empty ([number_of_models, parameter_set_dim])
        all_stabilities = np.empty ([number_of_models])

        J_partition = np.median (stabilities_0)  # <--- Create class partition based on this eigenvalue
        count0, count1 = 0, 0

        for i in range (0, number_of_models) :

            if i%100==0 :
                print (f'current set processed: {i}')
            this_param_set = f'vsyn{i}'
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
        savepath = f'result/turnover_test/gan_turnover_{N}_{miu}/gan_input/{exp_id}'
        os.makedirs (savepath, exist_ok=True)
        np.save (f'{savepath}/all_vsyn_{exp_id}.npy', all_data)
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
        savepath = f'result/turnover_test/gan_turnover_{N}_{miu}/gan_output/{exp_id}/'
        os.makedirs (savepath, exist_ok=True)
        for k in range (0, repeats) :
            print (f'Current exp: {exp_id}: Samples used: {np.shape (X_train)[0]}, repeat {k}')
            # set save directory
            this_savepath = f'{savepath}repeat_{k}/'
            os.makedirs (this_savepath, exist_ok=True)

            cgan = CGAN (X_train, y_train, latent_dim, batchsize, path_generator, savepath=this_savepath)
            d_loss, g_loss, acc = cgan.train (epochs, sample_interval, n_sample)

            # store training summary
            if np.mean (g_loss) < best_gen_loss:
                best_gen_loss = np.mean (g_loss)
            best_gan_output_path = this_savepath

            this_train_savepath = f'{this_savepath}training_summary/'
            os.makedirs (this_train_savepath, exist_ok=True)

            with open (f'{this_train_savepath}d_loss.pkl', 'wb') as f :
                pickle.dump (d_loss, f)
            with open (f'{this_train_savepath}g_loss.pkl', 'wb') as f :
                pickle.dump (g_loss, f)
            with open (f'{this_train_savepath}acc.pkl', 'wb') as f :
                pickle.dump (acc, f)

        # plot metrics
        x_plot = np.arange (0 , epochs , 1)
        fig = plt.figure (figsize=(10 , 5))
        ax1 = fig.add_axes ([0.2 , 0.2 , 1 , 1])
        ax1.plot (x_plot , d_loss , label='discriminator loss')
        ax1.plot (x_plot , g_loss , label='generator loss')
        ax1.set (ylabel='criterion_losses' , xlabel='epochs')
        ax1.legend ()
        plt.savefig (f'{this_train_savepath}loss.svg' , dpi=300 ,
                     transparent=False , bbox_inches='tight')
        plt.close(fig)

        fig = plt.figure (figsize=(10 , 5))
        ax2 = fig.add_axes ([0.2 , 0.2 , 1 , 1])
        ax2.plot (x_plot , d_loss , label='discriminator accuracy')
        ax2.set (ylabel='accuracy' , xlabel='epochs')
        ax2.legend ()
        plt.savefig (f'{this_train_savepath}d_accuracy.svg' , dpi=300 ,
                     transparent=False , bbox_inches='tight')
        plt.close(fig)

        # load new kcat
        vsyn_new_name = f'{best_gan_output_path}490_r.npy'
        vsyn_new = np.load (vsyn_new_name)
        vsyn_new = np.exp (vsyn_new)
        vsyn_new = np.transpose (vsyn_new)
        vsyn_df = pd.DataFrame (vsyn_new)
        column_names = [f'vsyn_value{i}' for i in range (vsyn_df.shape[1])]
        vsyn_df.columns = column_names
        vsyn_final_tmp = pd.concat ([turnover_data, vsyn_df], axis=1)
        print (vsyn_final_tmp)

        # calculate new stabilities
        num_cpus = 16
        print (num_cpus)
        with Pool (processes=num_cpus) as pool :
            args = [(
                j, growthdata, eModel_gan, model, vsyn_final_tmp, proteomicsdata_Clim, proteomicsdata_Nlim, miu) for j
                in range (number_of_models)]
            result = pool.starmap (parallel_turnover, args)
            print (result)
        pool.close ()
        pool.join ()
        results = np.vstack (result)
        toy_data = pd.DataFrame (results, columns=['error', 'corr', 'pvalue', 'number'])
        toy_data[ 'corr' ].fillna (toy_data[ 'corr' ].min () , inplace=True)

        # save the result and sort to new inputs
        error_new = (toy_data['error'] - toy_data['error'].min ())/(toy_data['error'].max () - toy_data['error'].min ())
        corr_new = 1- (toy_data['corr'] - toy_data['corr'].min ())/(toy_data['corr'].max () - toy_data['corr'].min ())
        number_new = 1- (toy_data['number'] - toy_data['number'].min ())/(toy_data['number'].max () - toy_data['number'].min ())
        w1 , w2 , w3 = 0 , 1 , 0 # weight
        stability_new = w1*error_new + w2*corr_new + w3*number_new
        # stability_new = 1-toy_data['corr']

        if j > 0:
            stability_old = pd.read_csv (file_path_sta).iloc[: , 0].values
        else:
            stability_old = pd.read_csv (path_stability).iloc[: , 0].values

        stability_all = np.append (stability_old, stability_new)
        stability_sort = sorted (range (len (stability_all)), key=lambda k : stability_all[k], reverse=False)
        stabilities_0 = stability_all[stability_sort[0 :1000]]
        savepath_sta = f'result/turnover_test/gan_turnover_{N}_{miu}/gan_output/result_summary'
        os.makedirs (savepath_sta, exist_ok=True)
        file_path_sta = os.path.join (savepath_sta, f'pro_data_{j}.csv')
        file_toy_data = os.path.join (savepath_sta, f'result_{j}.csv')
        np.savetxt (file_path_sta, stabilities_0)
        np.savetxt (file_toy_data, toy_data)

        vsyn_old = all_data
        vsyn_new = np.transpose (vsyn_new)
        vsyn_total = np.concatenate ((vsyn_old, vsyn_new))
        vsyn_next = vsyn_total[stability_sort[0 :1000], :]
        path_parameters = f'result/turnover_test/gan_turnover_{N}_{miu}/gan_input/test{j}/vsyn_data_1000.h5'
        with h5py.File (path_parameters, 'w') as h5f :
            for q in range (len (vsyn_next)) :
                # 获取第 q 个元素
                data_to_write = vsyn_next[q, :]
                h5f.create_dataset (f'vsyn{q}', data=data_to_write)
        f = h5py.File (path_parameters, 'r')
    end = time.time ()
    print (f'PROCESSING DONE in {end - start:.05} seconds')


if __name__=='__main__' :
    main ()
