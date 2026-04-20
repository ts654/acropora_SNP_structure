# -*- coding: utf-8 -*-
"""
Description:
    export a list of the distances between SNPs
    use both outputs from map_SNP_to_proteins.py
        merge lists and sort
        if sequential SNPs are on the same scaffold/chromosome but different proteins, then return distance
        
input
    Tim_Aten_Ning_bl_SNP_in_proteins.tsv
    Tim_Aten_Ning_bl_SNP_intergenic.tsv
    
output
    .tsv list of distances and the associated SNPs

"""

# Dependencies
import os
import numpy as np

# Constants
CONDITIONS = ['bl', 'fvfm', 'larvsurv', 'surv']
NUMINPUT = len(CONDITIONS)

header_merge = ['SNP_num', 'scaffold_num', 'nt_position', 'in_protein/intergenic', 'in_protein', 'intergenic upstream', 'intergenic downsteam']
header_distance = ['distance', 'scaffold_num', 'SNP_1', 'SNP_1_position', 'SNP_1_category', 'SNP_2', 'SNP_2_position', 'SNP_2_category']

# Paths
dir_parent = r'path\to\folder\of\SNP\files'

path_files_protein = np.empty(NUMINPUT, dtype=object)
path_files_intergenic = np.empty(NUMINPUT, dtype=object)
for i in range(0, NUMINPUT):
    path_files_protein[i] = os.path.join(dir_parent, f'Tim_Aten_Ning_{CONDITIONS[i]}_SNP_in_proteins.tsv')
    path_files_intergenic[i] = os.path.join(dir_parent, f'Tim_Aten_Ning_{CONDITIONS[i]}_SNP_intergenic.tsv')


# Main
for i in range(0, NUMINPUT):
    # import list of mapped SNPs in proteins and intergenic
    list_protein = np.genfromtxt(path_files_protein[i], delimiter='\t', dtype='str')
    list_intergenic = np.genfromtxt(path_files_intergenic[i], delimiter='\t', dtype='str')
    
    # merge lists with needed information
    list_merge = np.empty([(list_protein.shape[0] - 1 + list_intergenic.shape[0] - 1), len(header_merge)], dtype=object) # both lists have headers
    count_merge = 0
    for j in range(1, list_protein.shape[0]): # has header
        list_merge[count_merge, 0] = list_protein[j, 0]
        list_merge[count_merge, 1] = list_protein[j, 1]
        list_merge[count_merge, 2] = int(list_protein[j, 2])
        list_merge[count_merge, 3] = 'in_protein'
        list_merge[count_merge, 4] = list_protein[j, 3]
        list_merge[count_merge, 5] = ''
        list_merge[count_merge, 6] = ''
        count_merge += 1

    for j in range(1, list_intergenic.shape[0]): # has header
        list_merge[count_merge, 0] = list_intergenic[j, 0]
        list_merge[count_merge, 1] = list_intergenic[j, 1]
        list_merge[count_merge, 2] = int(list_intergenic[j, 2])
        list_merge[count_merge, 3] = 'intergenic'
        list_merge[count_merge, 4] = ''
        list_merge[count_merge, 5] = list_intergenic[j, 7]
        list_merge[count_merge, 6] = list_intergenic[j, 13]
        count_merge += 1
    
    # sort by scaffold then position
    list_merge = np.asarray(sorted(list_merge, key=lambda x: x[2], reverse = False))
    list_merge = np.asarray(sorted(list_merge, key=lambda x: x[1], reverse = False))
    list_merge = np.insert(list_merge, 0, header_merge, axis=0)
    
    # measure distance between SNPs
    list_measure = []
    list_measure.append(header_distance)
    for j in range(2, list_merge.shape[0]): # has header, and compare to previous row
        measure = False
        if list_merge[j, 1] == list_merge[j - 1, 1]: # same scaffold
            if not list_merge[j, 3] == list_merge[j - 1, 3]: # if they're different categories, measure
                measure = True
            else: # if same category
                if not list_merge[j, 4] == list_merge[j - 1, 4]: # check if in different proteins
                    measure = True
                elif (not list_merge[j, 5] == list_merge[j - 1, 5]) and (not list_merge[j, 6] == list_merge[j - 1, 6]): # check if between different proteins
                    measure = True
        if measure:
            list_measure.append([list_merge[j, 2] - list_merge[j - 1, 2], list_merge[j - 1, 1], \
                                 list_merge[j - 1, 0], list_merge[j - 1, 2], list_merge[j - 1, 3],\
                                 list_merge[j, 0], list_merge[j, 2], list_merge[j, 3]])
    list_measure_save = np.asarray(list_measure)

    # save    
    np.savetxt(os.path.join(dir_parent, f'Tim_Aten_Ning_{CONDITIONS[i]}_SNP_all_mapped.tsv'), list_merge, fmt="%5s", delimiter="\t")
    np.savetxt(os.path.join(dir_parent, f'Tim_Aten_Ning_{CONDITIONS[i]}_SNP_distance.tsv'), list_measure_save, fmt="%5s", delimiter="\t")
      
print('done')