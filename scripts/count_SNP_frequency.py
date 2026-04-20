# -*- coding: utf-8 -*-
"""
Description:
    This script counts the number of SNPs that are associated with each protein,
    divided into categories based on whether it causes an amino acid change.
    It uses the outputs from map_SNP_to_proteins_amil.py or map_SNP_to_proteins_aten.py

input:
    analyzed_SNP_in_proteins.tsv
    analyzed_SNP_in_intergenic.tsv

output:
    .tsv with list of proteins and number of times it was in the list
"""

import os
import numpy as np

# paths
dir_parent = r'folder\with\SNP_mapping'

file_base = 'bl_SNP'
file_prot = f'{file_base}_in_proteins.tsv'
file_intergen = f'{file_base}_intergenic.tsv'
file_save = f'{file_base}_freq.tsv'

path_prot = os.path.join(dir_parent, file_prot)
path_intergen = os.path.join(dir_parent, file_intergen)
path_save = os.path.join(dir_parent, file_save)

header = ['gene', 'count in protein non-syn', 'count in protein no change', 'count in intergenic']
num_list = len(header) - 1

# Functions
def AddList(_path, _pos, _col_check, _col_save):
    global dict_count
    
    _list_input = np.genfromtxt(_path, delimiter='\t', dtype='str')
    
    #create dictionary of genes with count
    for i in range(1, _list_input.shape[0]): # file has header
        if '_' in _list_input[i, _col_check]: # '_' will indentify SNPs in protein and gene name in intergenic
            _temp_protein = _list_input[i, _col_save]
            
            if not _temp_protein in dict_count:
                dict_count[_temp_protein] = np.zeros(num_list)
                dict_count[_temp_protein][_pos] = 1
            else:
                dict_count[_temp_protein][_pos] += 1
        elif 'intron' in _list_input[i, _col_check] or 'synonymous' in _list_input[i, _col_check]: # in protein mutation that is not non-synonymous
            _temp_protein = _list_input[i, _col_save]
            
            if not _temp_protein in dict_count:
                dict_count[_temp_protein] = np.zeros(num_list)
                dict_count[_temp_protein][_pos + 1] = 1
            else:
                dict_count[_temp_protein][_pos + 1] += 1
    return


# Main
dict_count = {}

AddList(path_prot, 0, 10, 3)
AddList(path_prot, 0, 13, 3)
AddList(path_intergen, 2, 7, 7)
AddList(path_intergen, 2, 13, 13)

#convert to sorted np
np_output = np.empty([len(dict_count), len(header)], dtype=object)
for i, key in enumerate(dict_count):
    np_output[i, 0] = str(key)
    for j in range(0, num_list):
        np_output[i, j + 1] = int(dict_count[key][j])
np_output_sorted = np.asarray(sorted(np_output, key=lambda x: x[1], reverse=True))

#add header
np_output_sorted_head = np.insert(np_output_sorted, 0, header, axis=0)

# save
np.savetxt(path_save, np_output_sorted_head, fmt="%5s", delimiter="\t")
print('done')