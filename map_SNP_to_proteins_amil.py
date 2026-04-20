# -*- coding: utf-8 -*-
"""
Description:
    This script maps a list of SNPs with genome positions to the A. millepora annotated genome.
    This genome should be the same one used in the SNP analysis.

input:
    The chromosome files should be in a single folder. They should be named with the format: ch_01.gb
    The SNP file is .xlsx containing the following columns: SNP# [12], Chromosome# [1], nt-Position [2], Major allel [3], Minor allel [4]
    
output:
    if SNP is in a protein:
        _SNP_in_proteins.tsv: SNP#, Ch#, nt-position, reference DNA, reference protein, minor DNA, minor protein, minor mutation, major DNA (if different), major protein, major mutation
    if SNP is NOT in a protein:
        _SNP_intergenic.tsv: SNP#, Ch#, nt-position, gene upstream distance, gene upstream DNA, gene upstream protein, gene downstream distance, gene downstream DNA, gene downstream protein
"""

# import packages
import os
import pandas as pd
import numpy as np
from Bio import SeqIO


# Define paths and constants
dir_ch = r'path\to\chromosome\folder'

dir_parent = r'path\to\snp\folder'
file_snp = 'snp.xlsx'
file_save_base = file_snp[:-5]
path_snp = os.path.join(dir_parent, file_snp)

# indices in SNP input file
ISNP = 12
ICH = 1
IPOS = 2
IMAJOR = 3
IMINOR = 4

HEADER_list_snp_inter = ['SNP_num', 'Ch_num', 'nt-position', \
                     'upstream gene distance', 'upstream gene DNA', 'upstream gene protein', 'upstream gene protein genbank','upstream protein id', 'upstream protein description', \
                     'downstream gene distance', 'downstream gene DNA', 'downstream gene protein', 'downstream gene protein genbank', 'downstream protein id', 'downstream protein description']
HEADER_list_snp_prot = ['SNP_num', 'Ch_num', 'nt-position', 'protein id', 'protein description', \
                    'reference DNA', 'reference protein', 'reference protein genbank', \
                    'minor DNA', 'minor protein', 'minor mutation', \
                    'major DNA (if different)', 'major protein', 'major mutation']

OUTPUTINTER = len(HEADER_list_snp_inter)
OUTPUTPROT = len(HEADER_list_snp_prot)

# Functions
def FindChNum(_file):
    if '_' in _file and '.' in _file:
        _delim1 = _file.split('.')[0]
        _num = int(_delim1.split('_')[1])
    else:
        _num = -1    
    return _num

def SplitSearch(_value, _start, _end):
    # seach in the global variable ch (was set to the correct chromosome)
    # returns the index of the protein in front or equal to the snp position

    global ch

    #print(_value, ch[_start][0], ch[_end][0], _start, _end)

    if _end - _start == 0 or _end - _start == 1:
        if _value == ch[_start][0]:
            _index = _start
        elif _value == ch[_end][0]:
            _index = _end
        elif _value < ch[_end][0]:
            _index = _start
        else:
            _index = -1
    else:
        _middle = int(np.round(np.mean([_start, _end])))
        if _value > ch[_middle][0]:
            _index = SplitSearch(_value, _middle, _end)
        else:
            _index = SplitSearch(_value, _start, _middle)  
    return _index

def FindDNAProteinSeq(_rec, _feature):
    _cds_locations= _feature.location
    _cds_DNA = _cds_locations.extract(_rec.seq)
    _cds_protein = _cds_DNA.translate(to_stop=False)
    return [str(_cds_DNA), str(_cds_protein[:len(_cds_protein) - 1])]

def ChangeNT(_rec, _pos, _nt):
    _rec_snp = _rec
    _rec_snp.seq = _rec.seq[:_pos] + _nt + _rec.seq[_pos + 1:]
    return _rec_snp

def FindMutation(_ref, _snp):
    if _ref[0] == _snp[0]:
        _change = 'SNP is in intron'
    else:
        if _ref[1] == _snp[1]:
            _change = 'synonymous mutation'
        elif not len(_ref) == len(_snp):
            _change = 'truncation'
        else:
            for _i in range(0, len(_ref[1])):
                if not _ref[1][_i] == _snp[1][_i]:
                    _change = f'{_ref[1][_i]}_{_i + 1}_{_snp[1][_i]}' #residue number is index + 1
                    break
    return _change


## Main

# load all chromosome data into an array sorted by position [start, end, orientation (1/-1), genbank feature]
list_dir_ch = os.listdir(dir_ch)
list_dir_ch.sort()

numch = FindChNum(list_dir_ch[len(list_dir_ch)-1])
list_ch = np.empty(numch + 1, dtype='object')
list_ch_rec = np.empty(numch + 1, dtype='object')

for i in range(0, len(list_dir_ch)):
    tmp_ch_num = FindChNum(list_dir_ch[i])
    if not tmp_ch_num == -1:
        tmp_ch = []
        list_ch_rec[tmp_ch_num] = SeqIO.read(open(os.path.join(dir_ch, list_dir_ch[i]), 'r'), 'genbank')
        
        for feature in list_ch_rec[tmp_ch_num].features:
            if feature.type == "CDS":
                side5 = int(feature.location.start)
                side3 = int(feature.location.end)
    
                if side5 < side3:
                    start = side5
                    end = side3
                    orientation = 1 #forward, positive sense
                else:
                    start = side3
                    end = side5
                    orientation = -1 #reverse, negative sense
                tmp_ch.append([start, end, orientation, feature])
        tmp_ch.sort(key=lambda x: x[0]) #may not need this, gb should be sorted in order
        list_ch[tmp_ch_num] = tmp_ch
del tmp_ch
del tmp_ch_num


# load SNP data
df_snp = pd.read_excel(path_snp)
list_snp = df_snp.to_numpy()

# adjust nt position because python is base 0
for i in range(0, list_snp.shape[0]):
    list_snp[i, IPOS] = list_snp[i, IPOS] - 1

# define empty output arrays with header
list_snp_inter = np.empty([len(list_snp) + 1, OUTPUTINTER], dtype='object')
list_snp_prot = np.empty([len(list_snp) + 1, OUTPUTPROT], dtype='object')

# add header
list_snp_inter[0] = HEADER_list_snp_inter
list_snp_prot[0] = HEADER_list_snp_prot

# have a header, so start saving at index 1
count_inter = 1
count_prot = 1

for i in range(0, list_snp.shape[0]):
    if (list_snp[i, ICH] <= numch) and (not list_ch[list_snp[i, ICH]] == None):
        ch = list_ch[list_snp[i, ICH]]
        infront_index = SplitSearch(list_snp[i, IPOS], 0, len(ch) - 1) #index of the the protein in front of the nucleotide
        #print('final:', infront_index)
        
        # check if intergenic
        if list_snp[i, IPOS] > ch[infront_index][1]: # if snp position is after the end of the CDS
            intergenic = True
        else:
            intergenic = False
        
        # save to intergenic list
        if intergenic:
            seqs_up = FindDNAProteinSeq(list_ch_rec[list_snp[i, ICH]], ch[infront_index][3])
            seqs_down = FindDNAProteinSeq(list_ch_rec[list_snp[i, ICH]], ch[infront_index + 1][3])
            
            list_snp_inter[count_inter, 0] = list_snp[i, ISNP]
            list_snp_inter[count_inter, 1] = list_snp[i, ICH]
            list_snp_inter[count_inter, 2] = list_snp[i, IPOS]
            list_snp_inter[count_inter, 3] = list_snp[i, IPOS] - ch[infront_index][1]
            list_snp_inter[count_inter, 4] = seqs_up[0]
            list_snp_inter[count_inter, 5] = seqs_up[1]
            list_snp_inter[count_inter, 6] = ch[infront_index][3].qualifiers['translation'][0]
            list_snp_inter[count_inter, 7] = ch[infront_index][3].qualifiers['protein_id'][0]
            list_snp_inter[count_inter, 8] = ch[infront_index][3].qualifiers['product'][0]
            list_snp_inter[count_inter, 9] = ch[infront_index + 1][0] - list_snp[i, IPOS] 
            list_snp_inter[count_inter, 10] = seqs_down[0]
            list_snp_inter[count_inter, 11] = seqs_down[1]
            list_snp_inter[count_inter, 12] = ch[infront_index + 1][3].qualifiers['translation'][0]
            list_snp_inter[count_inter, 13] = ch[infront_index + 1][3].qualifiers['protein_id'][0]
            list_snp_inter[count_inter, 14] = ch[infront_index + 1][3].qualifiers['product'][0]
            
            count_inter += 1
            
        # save to in_protein list
        else:
            tmp_rec_ref = list_ch_rec[list_snp[i, ICH]]
            seqs_ref = FindDNAProteinSeq(tmp_rec_ref, ch[infront_index][3])
            
            if list_snp[i, IMAJOR] == tmp_rec_ref[list_snp[i, IPOS]]: # if the major allel is the same as the reference
                exportmajor = False
                seqs_major = ['', '']
                change_major = ''
            else:
                exportmajor = True
                tmp_rec_major = ChangeNT(tmp_rec_ref, list_snp[i, IPOS], list_snp[i, IMAJOR])
                seqs_major = FindDNAProteinSeq(tmp_rec_major, ch[infront_index][3])
                change_major = FindMutation(seqs_ref, seqs_major)
    
            tmp_rec_minor = ChangeNT(tmp_rec_ref, list_snp[i, IPOS], list_snp[i, IMINOR])
            seqs_minor = FindDNAProteinSeq(tmp_rec_minor, ch[infront_index][3])
            change_minor = FindMutation(seqs_ref, seqs_minor)
    
            list_snp_prot[count_prot, 0] = list_snp[i, ISNP]
            list_snp_prot[count_prot, 1] = list_snp[i, ICH]
            list_snp_prot[count_prot, 2] = list_snp[i, IPOS]
            list_snp_prot[count_prot, 3] = ch[infront_index][3].qualifiers['protein_id'][0]
            list_snp_prot[count_prot, 4] = ch[infront_index][3].qualifiers['product'][0]
            list_snp_prot[count_prot, 5] = seqs_ref[0]
            list_snp_prot[count_prot, 6] = seqs_ref[1]
            list_snp_prot[count_prot, 7] = ch[infront_index][3].qualifiers['translation'][0]
            list_snp_prot[count_prot, 8] = seqs_minor[0]
            list_snp_prot[count_prot, 9] = seqs_minor[1]
            list_snp_prot[count_prot, 10] = change_minor
            list_snp_prot[count_prot, 11] = seqs_major[0]
            list_snp_prot[count_prot, 12] = seqs_major[1]
            list_snp_prot[count_prot, 13] = change_major
    
            count_prot += 1
    else:
        print(f'Chromosome {list_snp[i, ICH]} not present in folder')
list_snp_inter = np.delete(list_snp_inter, range(count_inter, np.shape(list_snp_inter)[0]), 0)
list_snp_prot = np.delete(list_snp_prot, range(count_prot, np.shape(list_snp_prot)[0]), 0)
print(f'#SNP intergenic: {count_inter - 1}')
print(f'#SNP in proteins: {count_prot - 1}')

# save
np.savetxt(os.path.join(dir_parent, f'{file_save_base}_SNP_intergenic.tsv'), list_snp_inter, delimiter="\t", fmt="%5s")
np.savetxt(os.path.join(dir_parent, f'{file_save_base}_SNP_in_proteins.tsv'), list_snp_prot, delimiter="\t", fmt="%5s")

print('Done')
