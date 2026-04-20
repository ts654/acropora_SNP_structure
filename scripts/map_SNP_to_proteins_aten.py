# -*- coding: utf-8 -*-
"""
Description:
    This script maps a list of SNPs with genome positions to the A. tenuis annotated genome.
    This genome should be the same one used in the SNP analysis.

input:
    aten_0.11.genbank file - multi fasta of all chromosomes
        note: scaffold 375 is numbered 154
    The SNP file is .xlsx containing the following columns: SNP# [3], Scaffold# [2], nt-Position [1], Major allel [5], Minor allel [6]
    
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
from Bio.SeqFeature import SeqFeature, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# Define paths and constants
path_ch = r'path\to\aten_0.11.genbank'

dir_parent = r'path\to\snp\folder'
file_snp = 'snp.xlsx'
file_save_base = os.path.splitext(file_snp)[0]
path_snp = os.path.join(dir_parent, file_snp)

    
# indices in SNP input file
ISNP = 3 # identifier
ICH = 2 # scaffold
IPOS = 1 # position
IMAJOR = 5
IMINOR = 6

HEADER_list_snp_inter = ['SNP num', 'scaffold num', 'nt-position', \
                     'upstream gene distance', 'upstream gene DNA', 'upstream gene protein', 'upstream gene protein genbank','upstream protein id', 'upstream protein description', \
                     'downstream gene distance', 'downstream gene DNA', 'downstream gene protein', 'downstream gene protein genbank', 'downstream protein id', 'downstream protein description']
HEADER_list_snp_prot = ['SNP num', 'scaffold num', 'nt-position', 'protein id', 'protein description', \
                    'reference DNA', 'reference protein', 'reference protein genbank', \
                    'minor DNA', 'minor protein', 'minor mutation', \
                    'major DNA (if different)', 'major protein', 'major mutation']

OUTPUTINTER = len(HEADER_list_snp_inter)
OUTPUTPROT = len(HEADER_list_snp_prot)

# Functions
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
            _index = -99
    else:
        _middle = int(np.round(np.mean([_start, _end])))
        if _value > ch[_middle][0]:
            _index = SplitSearch(_value, _middle, _end)
        else:
            _index = SplitSearch(_value, _start, _middle)  
    return _index

def CreatedJoinedCDSFeature(_idx_feature, _idx_chrom):
    global chromosomes
    
    #join CDS exons
    _list_CDS_locations = []
    _ID = chromosomes[_idx_chrom].features[_idx_feature].qualifiers['ID']
    _protein_ID = str(_ID[0]).split('.CDS')[0]
    _strand = chromosomes[_idx_chrom].features[_idx_feature].location.strand
    _i = _idx_feature
    while  chromosomes[_idx_chrom].features[_i].type == 'CDS':
        _list_CDS_locations.append(chromosomes[_idx_chrom].features[_i].location)
        _i += 1
        if _i == len(chromosomes[_idx_chrom].features):
            break
        
    if len(_list_CDS_locations) > 1:
        if _strand == 1:
            _list_CDS_locations.sort(key=lambda x: int(x.start), reverse=False)
        else:
            _list_CDS_locations.sort(key=lambda x: int(x.start), reverse=True)
        _locations_joined = CompoundLocation(_list_CDS_locations)
        
        #find translation
        _cds_DNA = _locations_joined.extract(chromosomes[_idx_chrom].seq)
        _cds_protein = _cds_DNA.translate(to_stop=True)
        
        _qualifiers = {'ID': _ID,
                       'protein_id': [_protein_ID],
                       'translation': [_cds_protein]}
        _feature_CDS_joined = SeqFeature(_locations_joined, type="CDS", qualifiers=_qualifiers)
    else:
        _cds_locations = chromosomes[_idx_chrom].features[_idx_feature].location
        _cds_DNA = _cds_locations.extract(chromosomes[_idx_chrom].seq)
        _cds_protein = _cds_DNA.translate(to_stop=True)
        
        _feature_CDS_joined = chromosomes[_idx_chrom].features[_idx_feature]
        _feature_CDS_joined.qualifiers['protein_id'] = [_protein_ID]
        _feature_CDS_joined.qualifiers['translation'] = [_cds_protein]
    return _feature_CDS_joined

def FindDNAProteinSeq(_rec, _feature):
    _cds_locations= _feature.location
    _cds_DNA = _cds_locations.extract(_rec.seq)
    _cds_protein = _cds_DNA.translate(to_stop=True)
    return [str(_cds_DNA), str(_cds_protein)]

def ChangeNT(_rec, _pos, _nt):
    _rec_snp = SeqRecord(Seq(""), id="", description="")
    _rec_snp.seq = _rec.seq[:_pos] + _nt + _rec.seq[_pos + 1:]
    return _rec_snp

def FindMutation(_ref, _snp):
    _same = True
    if _ref[0] == _snp[0]:
        _change = 'SNP is in intron'
    else:
        if _ref[1] == _snp[1]:
            _change = 'synonymous mutation'
        elif not len(_ref[1]) == len(_snp[1]):
            _change = 'truncation'
            _same = False
        else:
            for _i in range(0, len(_ref[1])):
                if not _ref[1][_i] == _snp[1][_i]:
                    _change = f'{_ref[1][_i]}_{_i + 1}_{_snp[1][_i]}' #residue number is index + 1
                    _same = False
                    break
    return _change, _same


## Main
print('script started')

# load chromosome data and extract CDS sorted by position [start, end, orientation (1/-1), genbank feature]
chromosomes = list(SeqIO.parse(path_ch, 'genbank'))
print('chromosomes loaded')
dict_ch = {}
dict_ch_rec = {}
for i_ch in range(0, len(chromosomes)):
    tmp_ch_num = chromosomes[i_ch].id
    dict_ch_rec[tmp_ch_num] = chromosomes[i_ch]
    tmp_ch = []
    for i_feat in range(0, len(chromosomes[i_ch].features)):
        if i_feat % 10 == 0:
            print(f'checking feature {i_feat} on scaffold {i_ch}') #i_ch = 375 not written

        if chromosomes[i_ch].features[i_feat].type == 'CDS' and not chromosomes[i_ch].features[i_feat - 1].type == 'CDS':
            feature_CDS_joined = CreatedJoinedCDSFeature(i_feat, i_ch)
            side5 = int(feature_CDS_joined.location.start)
            side3 = int(feature_CDS_joined.location.end)
            if side5 < side3: #forward, positive sense
                start = side5
                end = side3
                orientation = 1
            else: #reverse, negative sense
                start = side3
                end = side5
                orientation = -1
            tmp_ch.append([start, end, orientation, feature_CDS_joined])

        tmp_ch.sort(key=lambda x: x[0]) #may not need this, gb should be sorted in order
        dict_ch[tmp_ch_num] = tmp_ch


del tmp_ch
del tmp_ch_num


# load SNP data
df_snp = pd.read_excel(path_snp)
list_snp = df_snp.to_numpy()
print('SNP file loaded')

# adjust nt position because python is base 0
for i in range(0, list_snp.shape[0]):
    list_snp[i, IPOS] = int(list_snp[i, IPOS]) - 1

# define empty output arrays with header
list_snp_inter = np.empty([len(list_snp) + 1, OUTPUTINTER], dtype='object')
list_snp_prot = np.empty([len(list_snp) + 1, OUTPUTPROT], dtype='object')

# add header
list_snp_inter[0] = HEADER_list_snp_inter
list_snp_prot[0] = HEADER_list_snp_prot

# have a header, so start saving at index 1
count_inter = 1
count_prot = 1

count_mut = 0

for i in range(0, list_snp.shape[0]):
    if list_snp[i, ICH] in dict_ch:
        if i % 10 == 0:
            print(f'checking SNP {i} with ID {list_snp[i, ISNP]}')

        ch = dict_ch[list_snp[i, ICH]]
                
        # find index
        # check if SNP is in front of first CDS or last CDS, otherwise search
        if list_snp[i, IPOS] < ch[0][0]:
            infront_index = -1
        elif list_snp[i, IPOS] > ch[len(ch) - 1][1]:
            infront_index = -2
        elif list_snp[i, IPOS] > ch[len(ch) - 1][0]:
            infront_index = len(ch) - 1
        else:
            infront_index = SplitSearch(list_snp[i, IPOS], 0, len(ch) - 1) #index of the the protein in front of the nucleotide
        
        # check if intergenic
        if infront_index == -1 or infront_index == -2: # if snp position is after the end of the CDS
            intergenic = True
        elif list_snp[i, IPOS] > ch[infront_index][1]:
            intergenic = True
        else:
            intergenic = False
        
        # save to intergenic list
        if intergenic:
            list_snp_inter[count_inter, 0] = list_snp[i, ISNP]
            list_snp_inter[count_inter, 1] = list_snp[i, ICH]
            list_snp_inter[count_inter, 2] = list_snp[i, IPOS]

            if infront_index == -1: # before first gene
                seqs_down = FindDNAProteinSeq(dict_ch_rec[list_snp[i, ICH]], ch[infront_index + 1][3])

                list_snp_inter[count_inter, 9] = ch[infront_index + 1][0] - list_snp[i, IPOS] 
                list_snp_inter[count_inter, 10] = seqs_down[0]
                list_snp_inter[count_inter, 11] = seqs_down[1]
                list_snp_inter[count_inter, 12] = ch[infront_index + 1][3].qualifiers['translation'][0]
                list_snp_inter[count_inter, 13] = ch[infront_index + 1][3].qualifiers['protein_id'][0]
                list_snp_inter[count_inter, 14] = ch[infront_index + 1][3].qualifiers['ID'][0]
                
            elif infront_index == -2: # after last gene
                infront_index = len(ch) - 1
                seqs_up = FindDNAProteinSeq(dict_ch_rec[list_snp[i, ICH]], ch[infront_index][3])

                list_snp_inter[count_inter, 3] = list_snp[i, IPOS] - ch[infront_index][1]
                list_snp_inter[count_inter, 4] = seqs_up[0]
                list_snp_inter[count_inter, 5] = seqs_up[1]
                list_snp_inter[count_inter, 6] = ch[infront_index][3].qualifiers['translation'][0]
                list_snp_inter[count_inter, 7] = ch[infront_index][3].qualifiers['protein_id'][0]
                list_snp_inter[count_inter, 8] = ch[infront_index][3].qualifiers['ID'][0]
                
            else:
                seqs_up = FindDNAProteinSeq(dict_ch_rec[list_snp[i, ICH]], ch[infront_index][3])

                list_snp_inter[count_inter, 3] = list_snp[i, IPOS] - ch[infront_index][1]
                list_snp_inter[count_inter, 4] = seqs_up[0]
                list_snp_inter[count_inter, 5] = seqs_up[1]
                list_snp_inter[count_inter, 6] = ch[infront_index][3].qualifiers['translation'][0]
                list_snp_inter[count_inter, 7] = ch[infront_index][3].qualifiers['protein_id'][0]
                list_snp_inter[count_inter, 8] = ch[infront_index][3].qualifiers['ID'][0]

                seqs_down = FindDNAProteinSeq(dict_ch_rec[list_snp[i, ICH]], ch[infront_index + 1][3])

                list_snp_inter[count_inter, 9] = ch[infront_index + 1][0] - list_snp[i, IPOS] 
                list_snp_inter[count_inter, 10] = seqs_down[0]
                list_snp_inter[count_inter, 11] = seqs_down[1]
                list_snp_inter[count_inter, 12] = ch[infront_index + 1][3].qualifiers['translation'][0]
                list_snp_inter[count_inter, 13] = ch[infront_index + 1][3].qualifiers['protein_id'][0]
                list_snp_inter[count_inter, 14] = ch[infront_index + 1][3].qualifiers['ID'][0]
                
            count_inter += 1
            
        # save to in_protein list
        else:
            # extract sequences for ref, major, and minor
            tmp_rec_ref = dict_ch_rec[list_snp[i, ICH]]
            seqs_ref = FindDNAProteinSeq(tmp_rec_ref, ch[infront_index][3])

            tmp_rec_major = ChangeNT(tmp_rec_ref, list_snp[i, IPOS], list_snp[i, IMAJOR])
            seqs_major = FindDNAProteinSeq(tmp_rec_major, ch[infront_index][3])

            tmp_rec_minor = ChangeNT(tmp_rec_ref, list_snp[i, IPOS], list_snp[i, IMINOR])
            seqs_minor = FindDNAProteinSeq(tmp_rec_minor, ch[infront_index][3])

            # compare major and minor; if dif, export for mut prediction; if same, pick ref or SNP for annotation prediction
            change_min_maj, same_min_maj = FindMutation(seqs_minor, seqs_major)
            if not same_min_maj:
                count_mut += 1

            # compare to major and minor to ref for exporting
            if list_snp[i, IMAJOR] == tmp_rec_ref[list_snp[i, IPOS]]: # if the major allel is the same as the reference
                seqs_major = ['', '']
                change_major = 'same as ref'
            else:
                change_major, tmp = FindMutation(seqs_ref, seqs_major)
    
            if list_snp[i, IMINOR] == tmp_rec_ref[list_snp[i, IPOS]]: # if the minor allel is the same as the reference
                seqs_minor = ['', '']
                change_minor = 'same as ref'
            else:
                change_minor, tmp = FindMutation(seqs_ref, seqs_minor)

            list_snp_prot[count_prot, 0] = list_snp[i, ISNP]
            list_snp_prot[count_prot, 1] = list_snp[i, ICH]
            list_snp_prot[count_prot, 2] = list_snp[i, IPOS]
            list_snp_prot[count_prot, 3] = ch[infront_index][3].qualifiers['protein_id'][0]
            list_snp_prot[count_prot, 4] = ch[infront_index][3].qualifiers['ID'][0]
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
print(f'#point mutations: {count_mut}')

# save
np.savetxt(os.path.join(dir_parent, f'{file_save_base}_SNP_intergenic.tsv'), list_snp_inter, delimiter="\t", fmt="%5s")
np.savetxt(os.path.join(dir_parent, f'{file_save_base}_SNP_in_proteins.tsv'), list_snp_prot, delimiter="\t", fmt="%5s")

print('Done')
