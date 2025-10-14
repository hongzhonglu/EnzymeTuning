import pandas as pd
import numpy as np

#data preparation
data_all = pd.read_csv("G:/my_github/EnzymeTuning/data/E.coli/uniprot_data_all_ecoli.tsv",delimiter='\t')
gene_list = pd.read_csv("G:/my_github/EnzymeTuning/data/E.coli/gene_list.csv",header=None)
proteomics_data = pd.read_csv("G:/my_github/EnzymeTuning/data/E.coli/proteomic_ecoli.csv")

#compare data_all and proteomics_data
proteomics_all = pd.merge(data_all,proteomics_data,how='inner',left_on="Entry",right_on="Uniprot Accession")
# proteomics_all = data_all
# merge gene_list and proteomics_all
proteomics_all['Gene_ID'] = proteomics_all['Gene Names (ordered locus)'].str.split()
gene_id_exploded = proteomics_all.explode('Gene_ID')
gene_id_map = gene_id_exploded[gene_id_exploded['Gene_ID'].isin(gene_list[0])]
gene_id_map_final = gene_id_map.drop_duplicates()
gene_id_map_final.fillna(0, inplace=True)
gene_id_map_final.to_csv("G:/my_github/EnzymeTuning/data/E.coli/proteomics_ecoli_final.csv",index=False)