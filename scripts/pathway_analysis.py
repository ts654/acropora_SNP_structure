# -*- coding: utf-8 -*-
"""
Description:
    count the occurrences of the GO terms and keywords for a condition
        filter to only genes with number of SNPs above threshold
        filter keywords to top frequency keywords
    uses the output from count_SNP_frequency.py and keyword_frequency_sp_with_GO.py

input:
    mapped_SNP_freq.tsv : # gene, count in protein non-syn, count in protein no change, count in intergenic
        output from count_SNP_frequency.py
    foldseek_hits_freq.tsv : # gene, keywords"(#), ", top model, GO terms (delimiter ", ")
        output from keyword_frequency_sp_with_GO.py
    
output:
    GO_freq.tsv
    keyword_freq.tsv
    
"""

# dependencies
import os
import numpy as np

# constants
CONDITION = 'bl'
THRESH_SNP_NUM = 3.2
THRESH_KEY_NUM = 8 #3.6*2
EXCLUDE_WORDS = ['model_v4', 'domain-containing', 'family', 'member', 'putative']

# paths
path_SNP_freq = fr'path\to\file\{CONDITION}_SNP_freq.tsv'
path_keyword_freq = r'path\to\file\foldseek_hits_freq.tsv'

dir_save = r'path\to\save\folder'
path_save_go = os.path.join(dir_save, f'aten-{CONDITION}_freq_GO.tsv')
path_save_keyword = os.path.join(dir_save, f'aten-{CONDITION}_freq_{THRESH_KEY_NUM}_keyword.tsv')

# Function
def SaveDict(_dict, _path):
    _list = np.array(list(_dict.items()))
    _list = np.array(sorted(_list, key=lambda x: x[0], reverse=False)) # sort names
    _list = np.array(sorted(_list, key=lambda x: int(x[1]), reverse=True)) # sort most common at top
    
    np.savetxt(_path, _list, fmt="%s", delimiter="\t")
    return


# Main

# load data
list_SNP = np.genfromtxt(path_SNP_freq, delimiter='\t', dtype='str')
list_keyword = np.genfromtxt(path_keyword_freq, delimiter='\t', dtype='str')

# make list of genes to use based on number of SNP
list_genes_filtered = []
for i in range(1, list_SNP.shape[0]): # has header
    tempsum = int(list_SNP[i, 1]) + int(list_SNP[i, 2]) + int(list_SNP[i, 3])
    if tempsum > THRESH_SNP_NUM:
        list_genes_filtered.append(str(list_SNP[i, 0]))

# count annotation if gene is on list
dict_GO = {}
dict_keyword = {}
for gene in list_genes_filtered: # for each gene
    for j in range(0, list_keyword.shape[0]): # find the annotation entry
        if gene in list_keyword[j, 0]:
            # add to GO term list
            terms = list_keyword[j, 3].split(', ')
            for term in terms:
                if term not in dict_GO:
                    dict_GO[term] = 1
                else:
                    dict_GO[term] += 1
            
            # add to keyword list
            keywords_filtered = []
            keywords = list_keyword[j, 1].split(', ')
            for keyword in keywords:    # remove generic terms
                keyword_final = keyword.split('(')[0]
                if not 'model_v4' in keyword_final:
                    if not keyword_final in EXCLUDE_WORDS:
                        keywords_filtered.append(keyword_final)
            if len(keywords_filtered) < THRESH_KEY_NUM:  #if list size is less than THRESH_KEY_NUM, take all
                loopend = len(keywords_filtered)
            else:
                loopend = THRESH_KEY_NUM

            for k in range(0, loopend):
                keyword = keywords_filtered[k]
                if keyword not in dict_keyword:
                    dict_keyword[keyword] = 1
                else:
                    dict_keyword[keyword] += 1
            
            break

# save
SaveDict(dict_GO, path_save_go)
SaveDict(dict_keyword, path_save_keyword)
print('done')