from cobra.io import load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import numpy as np
import pandas as pd
from scipy.io import loadmat
import scipy.io as scio
def mat_to_xml(species):
    '''
    :param species: the name of the dl_model saved in 'data'
    :return: a dl_model in xml format
    '''
    modelfile = 'data/'+str(species)+'_dl.mat'
    model = cobra.io.load_matlab_model(modelfile)
    cobra.io.write_sbml_model(model,str(species)+'_dl.xml')


def mat_ec_to_ec_xml(species):
    '''
    :param species: the name of the ec_post_model saved in 'data'
    :return: a ec_post_model in xml format
    '''