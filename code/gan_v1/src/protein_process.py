import math
import numpy as np
import statistics
import pandas as pd
from src.mainFunction import *

def getSurfaceRatio(volume_ratio = 0.005/100, Vcell = 82, Scell=91.27):
    """
    The function is used to calculate the surface area ratio between organelle and cell.
    :param volume_ratio:
    :param Vcell:
    :param Scell:
    :return:
    """
    Vany = Vcell * volume_ratio
    Dany = 2 * (3 * Vany / (4 * math.pi)) ** (1 / 3)
    Sanym = 4 * math.pi * (Dany / 2) ** 2
    ratio = Sanym / Scell
    return ratio


def setReferenceProCopy():
    """
    This function is used to set the reference protein molecular/cell! If we still could not find the values from the
    reference value, then we can use the values at 5 percentile or 10 percentile?
    :return:
    """
    pro_abundance = pd.read_excel("data/proteomics/yeast_proteomics_example_cell_system_2018.xlsx")
    pro_abundance = pro_abundance[["Systematic Name", "Mean molecules per cell", "Median molecules per cell"]]
    pro_abundance.columns = ["gene", "absolute_abundance", "median_absolute_abundance"]  # protein abundance per cell
    pro_abundance.columns = ['gene', 'molecular/cell', "median_absolute_abundance"]
    pro_abundance = pro_abundance[pro_abundance["molecular/cell"].notna()]
    reference_copy = pro_abundance
    # statistics_analysis = reference_copy.describe()
    v_five_percent = reference_copy['molecular/cell'].quantile(0.05)
    v_ten_percent = reference_copy['molecular/cell'].quantile(0.1)
    return reference_copy, v_five_percent, v_ten_percent


def getProAundance_old_version(genes_select0, pro_abundance0):
    """
    Note: this function need double check!!!

    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param genes_select0:
    :param pro_abundance0: the unite is moleculars per cell
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = pd.DataFrame({"gene": genes_select0}) # change it as a dataframe
    combine_df["molecular/cell"] = singleMapping(pro_abundance0["molecular/cell"], pro_abundance0['gene'],combine_df["gene"])

    # for the protein without abundance, use the median value from this group.
    # calculate the abundance median value
    combine_df.fillna(value=pd.np.nan, inplace=True) # change none into nan in the dataframe
    abundance0 = combine_df["molecular/cell"].tolist()
    abundance1 = [x for x in abundance0 if np.isnan(x) == False]
    if len(abundance1) < 1:
        return "no_abundance"
    else:
        # use the first choice: for gene with no measured abundance, use the median value for the gene from the same compartment
        abundance_median = statistics.median(abundance1)  # here for the protein without abundance, the median value from this group is used. But maybe not correct at some cases
        abundance_update = []
        for x in abundance0:
            if np.isnan(x) == False:
                x0 = x
            else:
                x0 = abundance_median
            abundance_update.append(x0)
        combine_df["molecular/cell_local"] = abundance_update


        # use the second choice: for gene with no measured abundance, use the value from the reference conditions??
        # load the reference molecular copies
        ref_abundance, v_5, v_10 = setReferenceProCopy()
        # set a dict
        gene_abundance = {}
        for i, x in ref_abundance.iterrows():
            print(i)
            gene_abundance[x['gene']] = x['molecular/cell']

        abundance_update2 = []
        for i, x in combine_df.iterrows():
            print(i)
            ss = x['molecular/cell']
            if np.isnan(ss) == False:
                x0 = ss
            elif x['gene'] in ref_abundance['gene'].tolist():
                x0 = gene_abundance[x['gene']]
            else:
                x0 = v_5
            abundance_update2.append(x0)

        combine_df["molecular/cell_global"] = abundance_update2

    return combine_df


def getProAundance(genes_select0, pro_abundance0):
    """
    Note: this function need double check!!!

    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param genes_select0:
    :param pro_abundance0: the unite is moleculars per cell
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = pd.DataFrame({"gene": genes_select0}) # change it as a dataframe
    #combine_df["molecular/cell"] = singleMapping(pro_abundance0["molecular/cell"], pro_abundance0['gene'],combine_df["gene"])
    combine_df = pd.merge(left=pro_abundance0, right=combine_df, left_on=['gene'], right_on=['gene'],how="right")


    # for the protein without abundance, use the median value from this group.
    # calculate the abundance median value
    combine_df.fillna(value=pd.np.nan, inplace=True) # change none into nan in the dataframe
    abundance0 = combine_df["molecular/cell"].tolist()
    abundance1 = [x for x in abundance0 if np.isnan(x) == False]
    if len(abundance1) < 1:
        return "no_abundance"
    else:

        """
        # use the first choice: for gene with no measured abundance, use the median value for the gene from the same compartment
        abundance_median = statistics.median(abundance1)  # here for the protein without abundance, the median value from this group is used. But maybe not correct at some cases
        abundance_update = []
        for x in abundance0:
            if np.isnan(x) == False:
                x0 = x
            else:
                x0 = abundance_median
            abundance_update.append(x0)
        combine_df["molecular/cell_local"] = abundance_update


        # use the second choice: for gene with no measured abundance, use the value from the reference conditions??
        # load the reference molecular copies
        ref_abundance, v_5, v_10 = setReferenceProCopy()
        # set a dict
        gene_abundance = {}
        for i, x in ref_abundance.iterrows():
            print(i)
            gene_abundance[x['gene']] = x['molecular/cell']

        abundance_update2 = []
        for i, x in combine_df.iterrows():
            print(i)
            ss = x['molecular/cell']
            if np.isnan(ss) == False:
                x0 = ss
            elif x['gene'] in ref_abundance['gene'].tolist():
                x0 = gene_abundance[x['gene']]
            else:
                x0 = v_5
            abundance_update2.append(x0)

        combine_df["molecular/cell_global"] = abundance_update2

    return combine_df"""
        return combine_df


def getStructureSize(pro_size0, abundance0, need_check="No"):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = abundance0
    combine_df["Volume"] = singleMapping(pro_size0['Total_Volume'], pro_size0['locus'], combine_df["gene"])
    combine_df["section_area"] = singleMapping(pro_size0['section_area_new'], pro_size0['locus'], combine_df["gene"])
    #combine_df["molecular/cell"] = singleMapping(abundance0["molecular/cell"], abundance0['gene'], combine_df["gene"])
    # it shows that some genes have no locus
    combine_df = combine_df[~combine_df["section_area"].isna()]

    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_volume = sum(combine_df["molecular/cell_global"] * combine_df["Volume"])
    # change nm^3 into um^3
    total_volume_um = total_volume / 1e9

    # 1 纳米(nm)=0.001 微米(um)
    total_area = sum(combine_df["molecular/cell_global"] * combine_df["section_area"])
    # change nm^2 into um^2
    total_area_um = total_area / 1e6
    if need_check=="No":
        return total_volume_um, total_area_um
    else:
        combine_df["total_volume"] = combine_df["molecular/cell_global"] * combine_df["Volume"]
        combine_df["total_area"] = combine_df["molecular/cell_global"] * combine_df["section_area"]
        combine_df = combine_df.sort_values(by=['total_area'], ascending=False)
        return total_volume_um, total_area_um, combine_df


def getStructureSize_MeasuredAbundances(pro_size0, abundance0, need_check="No"):
    """
    The function is used to calculate the total protein size and sectional area for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param pro_size0: the unit is nm^3 (volume) or nm^2 (area)
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = abundance0
    combine_df["Volume"] = singleMapping(pro_size0['Total_Volume'], pro_size0['locus'], combine_df["gene"])
    combine_df["section_area"] = singleMapping(pro_size0['section_area_new'], pro_size0['locus'], combine_df["gene"])
    #combine_df["molecular/cell"] = singleMapping(abundance0["molecular/cell"], abundance0['gene'], combine_df["gene"])
    # it shows that some genes have no locus
    combine_df = combine_df[~combine_df["section_area"].isna()]
    combine_df = combine_df[~combine_df["molecular/cell"].isna()] # newly added for this new function


    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_volume = sum(combine_df["molecular/cell"] * combine_df["Volume"])
    # change nm^3 into um^3
    total_volume_um = total_volume / 1e9

    # 1 纳米(nm)=0.001 微米(um)
    total_area = sum(combine_df["molecular/cell"] * combine_df["section_area"])
    # change nm^2 into um^2
    total_area_um = total_area / 1e6
    if need_check=="No":
        return total_volume_um, total_area_um
    else:
        combine_df["total_volume"] = combine_df["molecular/cell"] * combine_df["Volume"]
        combine_df["total_area"] = combine_df["molecular/cell"] * combine_df["section_area"]
        combine_df = combine_df.sort_values(by=['total_area'], ascending=False)
        return total_volume_um, total_area_um, combine_df


def getOrganelleAbundance(abundance0, need_check="No"):
    """
    The function is used to calculate the total protein abundance for a group of genes from specific location.
    It should be noted that the unit of pro_abundance is molecules per cell.
    :param pro_abundance0: the unite is moleculars per cell
    :param need_check:
    :return:
    """

    # should make sure no structure size data is nan
    combine_df = abundance0
    #combine_df["molecular/cell"] = singleMapping(abundance0["molecular/cell"], abundance0['gene'], combine_df["gene"])
    # calculate the size of all proteins for the selected gene list
    # 1 纳米(nm)=0.001 微米(um)
    total_abundance = sum(combine_df["molecular/cell_global"])

    if need_check=="No":
        return total_abundance
    else:
        return combine_df



def splitAbundance(pro_df):
    """
    The function is used to quality check of protein abundance in molecular/cell or other unit before entering next step.
    :param pro_df: A dataframe should columns-gene,molecular/cell.
    :return:
    """
    # sometimes it shows mutiple proteins together have one abundance value, here we need a function to do the quality check!
    len1 = pro_df.shape[0]
    colnames=pro_df.columns
    pro_df1 = pro_df[pro_df['gene'].str.contains(';')]
    if (len(pro_df1) > 0):
        print('Mutiple protein have one abudance value! Need quality check.')
    pro_df2 = pro_df[~pro_df['gene'].str.contains(';')]
    gene0 = []
    abundance0 = []
    for i, x in pro_df1.iterrows():
        print(i, x)
        s = x["gene"].split(";")
        len0 = len(s)
        gene0 = gene0 + s
        v = [x[colnames[1]] / len0] * len0
        abundance0 = abundance0 + v
    gene0 = [x.strip(" ") for x in gene0]
    pro_df1 = pd.DataFrame({"gene": gene0, colnames[1]: abundance0})

    pro_df = pd.concat([pro_df1, pro_df2], axis=0)
    len2 = pro_df.shape[0]

    if (len2 > len1):
        print('Complete the quality check!')

    return pro_df



def getGoTermGeneList(input1, input2):
    """
    This function is generate the GO term and its gene list
    :parameter1: a directory contain a excel file with gene GO term annotation
    :parameter2: a directory contain a excel file with gene id mapping
    :return: A dict. GO_term as key and gene list as value
    ------------
    Usage:
    GO_term_gene = getGoTermGeneList(input1="data/pnas.1921890117.sd01_GO_term.xlsx", input2="data/sce_protein_weight.tsv")

    Hongzhong Lu
    2022.01.07
    """
    # Input the datasets from paxDB
    GO_term = pd.read_excel(input1)
    GO_term['GO-slim mapper process term'] = GO_term['GO-slim mapper process term'].str.strip()

    gene_info = pd.read_csv(input2, sep="\t")
    gene_short_name = gene_info['gene_name'].tolist()
    gene_short_name = [str(x) for x in gene_short_name]
    gene_info['gene_name'] = gene_short_name
    gene_short_name2 = []
    gene_locus = gene_info['locus'].tolist()
    for x, y in zip(gene_short_name, gene_locus):
        print(x, y)
        if x == 'nan':
            gene_short_name2.append(y)
        else:
            gene_short_name2.append(x)
    # build a dict
    # it shows that some genes have no locus
    # be careful in this step
    gene_name_dict = {}
    for w, v in zip(gene_short_name2, gene_locus):
        gene_name_dict[w] = v

    # build GO_term dict
    GO_dict = {}
    for i, x in GO_term.iterrows():
        print(i)
        ss = list(x)
        name = ss[0]
        ss = ss[1:]
        mylist = [str(x) for x in ss]
        newlist = [v for v in mylist if v != 'nan']
        # get the gene OFR name base on short gene ID
        key_all = gene_name_dict.keys()
        newlist1 = []
        for x in newlist:
            if x in key_all:
                x0 = gene_name_dict[x]
            else:
                x0 = x
            newlist1.append(x0)
        GO_dict[name] = newlist1
    return GO_dict



def getGeneListFromLocation(gene_location_annotation, location):
    """
    The function is to extract gene list based on its compartment information
    :param gene_location_annotation: A dataframe contains the annotation of each protein
    :param location: A string represent the compartment name
    :return:
    """
    #location = 'mitochondrial envelope'
    gene_subset = gene_location_annotation[gene_location_annotation["GO_Name"]==location]
    gene_list = list(set(gene_subset["Systematic_name"].tolist()))
    return gene_list



def getCompartmentGeneList(filter="Yes"):
    """
    This function to build a compartment dict, with which we can get the gene list from the compartment name

    :param filter:
    :return:
    """

    # Input the datasets from paxDB
    compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')

    # extract compartment
    compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'GO_Qualifier',
                           'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']
    compartment1 = compartment[compartment["GO_Namespace"] == "cellular_component"]

    # filter out compartment with "complex" or "subunit"
    compartment2 = compartment1[~compartment1["GO_Name"].str.contains("complex")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("subunit")]
    # firstly remove some general cellular component
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("snRNP")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("spindle")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("actin")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("cellular_component")]

    # analyze the annotation type
    annotation_type = compartment2["Annot_Type"].tolist()
    annotation_type = list(set(annotation_type))
    # here if we remove "computational"
    compartment_with_evidence = compartment2[compartment2["Annot_Type"] != 'computational']
    compartment_with_computation = compartment2[compartment2["Annot_Type"] == 'computational']
    # in one procedure, if a protein has no compartment annotation from manual and high-throughput, then the computational is used!
    compartment_addition = compartment_with_computation[~compartment_with_computation["Systematic_name"].isin(compartment_with_evidence["Systematic_name"])]
    compartment_combine = pd.concat([compartment_with_evidence, compartment_addition])

    # build the dict
    compartment_dict_all = {}
    for i, x in compartment2.iterrows():
        print(i, x)
        if x['GO_Name'] in compartment_dict_all.keys():
            compartment_dict_all[x['GO_Name']] = list(set(compartment_dict_all[x['GO_Name']] + [x["Systematic_name"]]))
        else:
            compartment_dict_all[x['GO_Name']] = list(set([x["Systematic_name"]]))
    # filter
    compartment_dict_all0 = {}
    for key in compartment_dict_all.keys():
        print(key)
        value = compartment_dict_all[key]
        if len(value) >= 6:
            compartment_dict_all0[key] = value
        else:
            pass

    compartment_dict2 = {}
    for i, x in compartment_combine.iterrows():
        print(i, x)
        if x['GO_Name'] in compartment_dict2.keys():
            compartment_dict2[x['GO_Name']] = list(set(compartment_dict2[x['GO_Name']] + [x["Systematic_name"]]))
        else:
            compartment_dict2[x['GO_Name']] = list(set([x["Systematic_name"]]))
    # filter
    compartment_dict20 = {}
    for key in compartment_dict2.keys():
        print(key)
        value = compartment_dict2[key]
        if len(value) >= 6:
            compartment_dict20[key] = value
        else:
            pass
    if filter == "Yes":
        return compartment_dict20
    else:
        return compartment_dict_all0



def AllProteomicsAnalysis(pro_df):
    """
    The function is used to do the general statistical analysis of proteomics datasets across conditions.
    :param pro_df: A dataframe to store proteomics, with column "gene"
    :return:

    _____

    usage: AllProteomicsAnalysis(pro_df=protein_copy_all1)
    """
    protein_copy_all2 = pro_df.drop(columns=['gene'])
    ss = protein_copy_all2.describe()
    # get the top 1000 proteins based on their molecular copies
    new_df = pro_df[['gene']]
    # Get the top 1000 proteins
    column20 = protein_copy_all2.columns
    for x in column20:
        print(x)
        df = pro_df[['gene', x]]
        ss2 = df.sort_values(by=[x], ascending=False)
        ss2_top1000 = ss2.iloc[0:1000, ]
        df[x][~df['gene'].isin(ss2_top1000['gene'])] = None
        new_df[x] = df[x]
    # analysis
    new_df1 = new_df[column20]
    ss1 = new_df1.describe()
    ss1.to_excel("data/proteomics/protein_copy_statistical_top1000.xlsx")
    ss.to_excel("data/proteomics/protein_copy_statistical.xlsx")
    return ss


def calculateCoefficient(cell_volume0):
    cell_volume = cell_volume0  # 32.6 # fL/cell
    dry_content = 0.35  # https://onlinelibrary.wiley.com/doi/pdf/10.1002/j.2050-0416.1952.tb02660.x#:~:text=Yeast%20cakes%20produced%20by%20normal,the%20conditions%20of%20growth%2C%20to
    cell_density = 1.1126e-12  # g/fL yeast cell density under exponential growth, [g/fL] = 1e12  g/mL
    coefficent10 = 1000 / 6.022e+23 / cell_volume / dry_content / cell_density  # from molecular/cell into mmol/gDW
    coefficent20 = 1 / coefficent10  # from mmol/gDW into molecular/cell
    return coefficent20


def calculateCurationCoefficent():
    # curation of rosemary datasets based on the fitted cell volume under different growth rates
    # fitting formula to calculate the coefficients
    # when miu > 0.2, cell_volume = 77.32 miu + 15.771
    # when miu < 0.2, average volume is 28 um^3

    # Rosemary sample ID
    Sample_ID_select = ['prot.1','prot.2', 'prot.3','prot.7','prot.8','prot.9','prot.10','prot.11','prot.12','prot.13','prot.14','prot.15','prot.16','prot.17','prot.18','prot.19','prot.20','prot.21']
    coefficient2 = 7.8298e9  # this is the original coefficient used in cell system paper!
    growth_rate = [0.05, 0.1, 0.13, 0.18, 0.3, 0.35]
    #cell_volume_fit = [28, 28, 28, 28, 38.967, 42.833] # here assume there exist a minimum cell size
    cell_volume_fit = [x*47.458+22.742 for x in growth_rate] # here we assume there exist a linear increase of cell size when growth rate increased
    coefficent_list = [calculateCoefficient(cell_volume0=x) for x in cell_volume_fit]
    curation_coefficent = [x / coefficient2 for x in coefficent_list]
    new_coefficent = []
    for x in curation_coefficent:
        print(x)
        s = [x] * 3
        new_coefficent = new_coefficent + s
    cell_volume_all = []
    for x in cell_volume_fit:
        print(x)
        s = [x] * 3
        cell_volume_all = cell_volume_all + s
    growth_all = []
    for x in growth_rate:
        print(x)
        s = [x] * 3
        growth_all = growth_all + s

    curation_info_rosemary = pd.DataFrame({"ID":Sample_ID_select,"growth_rate":growth_all, "cell_size":cell_volume_all, "curation_coefficent": new_coefficent})

    return curation_info_rosemary


def linearFit(df, x_name, y_name):
    """
    This function is used to fitting a linear relation between two variables
    Linear fit for two colums in a dataframe
    :param df:
    :param x_name:
    :param y_name:
    :return:
    """
    from sklearn.metrics import r2_score
    import matplotlib.pyplot as plt
    df = df[[x_name, y_name]]
    df = df.dropna()
    x = df[x_name]
    y = df[y_name]
    x_name = x_name.split("(")[0]
    y_name = y_name.split("(")[0]
    x_name = x_name.replace("/", "_per_")
    y_name = y_name.replace("/", "_per_")
    coef = np.polyfit(x, y, 1)
    poly1d_fn = np.poly1d(coef)
    predict = np.poly1d(coef)
    R2 = r2_score(y, predict(x))
    print(R2)
    print(coef)
    R2 = "{:.3f}".format(R2)
    # poly1d_fn is now a function which takes in x and returns an estimate for y
    plt.figure()
    plt.plot(x, y, 'yo', x, poly1d_fn(x), '--k')  # '--k'=black dashed line, 'yo' = yellow circle marker
    plt.xlabel(x_name, fontsize=18)
    plt.ylabel(y_name, fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    x_max = max(x)
    y_max = max(y)
    plt.text(x_max/3, 2*y_max/3, "R2=" + str(R2), fontsize=18)
    plt.show()
    return coef[0], coef[1]


def collectOrganelleTerm(type):
    """
    Some compartment need manual check.
    The function is just to get the important organelle list for volume or membrane size calculation.
    :param type:
    :return:

    usage:
    organelle_v = collectOrganelleTerm(type="volume")
    organelle_m = collectOrganelleTerm(type="m")

    """
    volume_list = ['mitochondrion', 'nucleus', 'cytosol', 'endoplasmic reticulum', 'endosome', 'lipid droplet',
                   'fungal-type vacuole', 'peroxisome', 'ribosome', 'Golgi apparatus', 'cytosolic ribosome',
                   'mitochondrial ribosome', 'nucleolus']
    membrane_list = ['fungal-type vacuole membrane', 'plasma membrane', 'mitochondrial outer membrane',
                     'prospore membrane', 'endoplasmic reticulum membrane', 'mitochondrial inner membrane',
                     'Golgi membrane', 'cellular bud membrane',
                     'late endosome membrane', 'peroxisomal membrane', 'nuclear membrane', 'endosome membrane',
                     'nuclear inner membrane']
    if type == "volume":
        return volume_list
    else:
        return membrane_list


def FingGenesForOrganelle(gene_set, compartment_list, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein volume or sectional area as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
        # all_compartment = list(compartment.keys())

    # use some manually checked gene compartment definion
    gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
    all_compartment = compartment_list
    result_df = dict()
    for y in all_compartment:
            print(y)
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            else:
                genes_select = compartment[y]
            # here we need calculate the intersection
            result_df[y] = list(set(genes_select) & set(gene_set))
    return result_df


def Pro3DCal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein volume or sectional area as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
        all_compartment = list(compartment.keys())

    # input the protein structure information
    pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
    pro_size = pro_size[['DBID', 'locus', 'Total_Volume', 'section_area_new']]
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]

    # use some manually checked gene compartment definion
    gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")

    # creat two dataframe to save the result
    result1 = pd.DataFrame({"compartment": all_compartment})
    result2 = pd.DataFrame({"compartment": all_compartment})

    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        value2 = []
        for y in all_compartment:
            print(y)
            pro_abundance = protein_copy[['gene', col0]]
            pro_abundance.columns = ['gene', 'molecular/cell']
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            else:
                genes_select = compartment[y]
            pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
            if pro_abundance1 is "no_abundance":
                value1.append(None)
                value2.append(None)
            else:
                x, S = getStructureSize_MeasuredAbundances(pro_size0=pro_size, abundance0=pro_abundance1)
                value1.append(x)
                value2.append(S)
        result1[col0] = value1
        result2[col0] = value2
    return result1, result2


def ProAbsoluteCal(protein_copy, compartment_type="organelle"):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_copy:
    :param compartment_type:
    :return:
    """
    if compartment_type == "organelle":
        # compartment info
        compartment = getCompartmentGeneList(filter="Yes")  # based on the automatic way
        all_compartment = list(compartment.keys())

    # input the protein structure information
    pro_size = pd.read_excel("result/sce_protein_size_3D_structure.xlsx")
    pro_size = pro_size[['DBID', 'locus', 'Total_Volume', 'section_area_new']]
    # sample ID information
    Sample_ID_select = list(protein_copy.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]

    # use some manually checked gene compartment definion
    gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
    # creat a dataframe to save the result
    result1 = pd.DataFrame({"compartment": all_compartment})
    # run the cycle
    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        for y in all_compartment:
            print(y)
            pro_abundance = protein_copy[['gene', col0]]
            pro_abundance.columns = ['gene', 'molecular/cell']
            if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            else:
                genes_select = compartment[y]
            pro_abundance1 = getProAundance(genes_select0=genes_select, pro_abundance0=pro_abundance)
            if pro_abundance1 is "no_abundance":
                value1.append(None)
            else:
                pro_abundance1 = pro_abundance1.dropna()
                x = sum(pro_abundance1['molecular/cell'])
                value1.append(x)
        result1[col0] = value1
    return result1

