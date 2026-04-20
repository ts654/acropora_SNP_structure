#-*- coding: utf-8 -*-
"""
Description:
    Analyzes the AlphaFold 3 outputs and summaries the scores into a single file.

input:
    output folder (folder of folder_prediction)

output:
    file with scores
    folder of images of PAE, pLDDT, and ipTMpair and model.cif
   
"""

# install/import dependencies and define functions
import numpy as np
import pandas as pd
import json
import os
from pathlib import Path
import shutil
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

NUMSEED = 1
NUMMODEL = 5

# file paths
dirsourcestr = r'path\to\folder\of\folder_predictions'
dirsource  = Path(dirsourcestr)
folder = dirsource.name
dirparent = dirsource.parent

folderanalyze = f'{folder}_analyze'
diranalyze = os.path.join(dirparent, folderanalyze)
if os.path.exists(diranalyze)==False:
    os.mkdir(diranalyze)

folderfiles = f'{folder}_files'
dirfilessave = os.path.join(diranalyze, folderfiles)
if os.path.exists(dirfilessave)==False:
    os.mkdir(dirfilessave)

prefix = 'TS5Scores'
filepath = os.path.join(diranalyze, f'{prefix}_{folder}.xlsx')

numoutput = 15
"""
 cif old name
 cif new name
 "fraction_disordered"
 "has_clash"
 "iptm"
 "chain_pair_iptm"
 "ptm"
 "chain_ptm"
 "ranking_score"
 mean pLDDT
 stdev pLDDT
 "30-PAEmin"
 "PAEstd/PAEmean"
 "chain_30-PAEmin"
 "chain_PAEstd/PAEmean"
"""


# Main
listdir = os.listdir(dirsource)
listdir.sort()
scores = np.zeros(((len(listdir)*NUMMODEL*NUMSEED), numoutput), dtype=object)
print('\nsave in folder:', diranalyze)
print("proteins found: ", len(listdir))

index = 0
for folder_num in range (0, len(listdir)): # for every prediction
    if folder_num%1==0:
        print(f'currently at model: {folder_num + 1} of {len(listdir)}, model: {listdir[folder_num]}')

    dirtemp = os.path.join(dirsource, listdir[folder_num])
    if os.path.isfile(os.path.join(dirtemp, 'ranking_scores.csv')):
        listrank = np.genfromtxt(os.path.join(dirtemp, 'ranking_scores.csv'), delimiter=",", dtype='str')
        listranksort = np.asarray(sorted(listrank[1:, :], key=lambda x: x[2], reverse = True)) # remove header
            
        for model_num in range(0, listranksort.shape[0]):
            rank_num = model_num + 1
            
            # pad rank_num
            rank_num_str = str(rank_num)
            if len(rank_num_str) == 1:
                rank_num_pad = f'00{rank_num_str}'
            elif len(rank_num_str) == 2:
                rank_num_pad = f'0{rank_num_str}'
            else:
                rank_num_pad = rank_num_str
    
            # extract scores
            scores[index, 0] = f'{listdir[folder_num]}_rank_{rank_num_pad}'
            scores[index, 1] = listdir[folder_num]
    
            # copy .cif with rank number
            if not os.path.isfile(os.path.join(dirfilessave, f'{scores[index, 1]}_rank_{rank_num_pad}.cif')):
                shutil.copy2(os.path.join(dirtemp, f'seed-{listranksort[model_num, 0]}_sample-{listranksort[model_num, 1]}', 'model.cif'), os.path.join(dirfilessave, f'{scores[index, 1]}_rank_{rank_num_pad}.cif'))
    
    
            # read json files
            f = open(os.path.join(dirtemp, f'{listdir[folder_num]}_data.json'))
            metadata_request = json.load(f)
            
            f = open(os.path.join(dirtemp, f'seed-{listranksort[model_num, 0]}_sample-{listranksort[model_num, 1]}', 'confidences.json'))
            metadata_data = json.load(f)
            
            f = open(os.path.join(dirtemp, f'seed-{listranksort[model_num, 0]}_sample-{listranksort[model_num, 1]}', 'summary_confidences.json'))
            metadata_conf = json.load(f)
            
    
            # extract scores directly
            scores[index, 2] = metadata_conf['fraction_disordered']
            scores[index, 3] = metadata_conf['has_clash']
            scores[index, 4] = metadata_conf['iptm']
            scores[index, 5] = str(metadata_conf['chain_pair_iptm'])
            scores[index, 6] = metadata_conf['ptm']
            scores[index, 7] = str(metadata_conf['chain_ptm'])
            scores[index, 8] = metadata_conf['ranking_score']
    
            iptmpair = np.asarray(metadata_conf['chain_pair_iptm'], dtype=float)
            
            # find total length of proteins
            specieslength = [] # array of each species length
            totalproteinchains = 0
            totalproteinlength = 0
            for species in metadata_request['sequences']:
                if 'protein' in species:
                    molecule = species['protein']['sequence']
                    specieslength.append(len(molecule))
                    
                    protein = species['protein']['sequence']
                    totalproteinchains = totalproteinchains + 1
                    totalproteinlength = totalproteinlength + len(protein)
                elif 'dna' in species:
                    molecule = species['dna']['sequence']
                    specieslength.append(len(molecule))
                elif 'rna' in species:
                    molecule = species['rna']['sequence']
                    specieslength.append(len(molecule))
                elif 'ligand' in species:
                    if 'smiles' in species['ligand']:
                        molecule = species['ligand']['smiles']
                    elif 'ccdCodes' in species['ligand']:
                        molecule = species['ligand']['ccdCodes']
                    else:
                        print('problem', species['ligand'])
                    moleculeatomletters = ''.join(char for char in molecule if char.isalpha())
                    moleculeatomall = ''.join(char for char in moleculeatomletters if char.isupper()) #remove lower case letters (Br is only 1 atom)
                    moleculeatom = ''.join(char for char in moleculeatomall if not char == 'H') #remove hydrogens
                    specieslength.append(len(moleculeatom))
            numspecies = len(specieslength)
        
            #extract plddt plot for total proteins
            plddtall = metadata_data['atom_plddts']
            plddtid = metadata_data['atom_chain_ids']
            
            chaincount = 1
            atomlength = -1
            for j in range(0, len(plddtid)):
                if not j == len(plddtid)-1:
                    if not plddtid[j] == plddtid[j+1]:
                        chaincount += 1
                if chaincount > totalproteinchains:
                    atomlength = j
                    break
            if atomlength == -1:
                atomlength = len(plddtall)
            plddt = plddtall[0:atomlength]
             
        
            #plddt calculations
            scores[index, 9] = np.mean(plddt)
            scores[index, 10] = np.std(plddt)
    
    
            #extract pae plot for total proteins
            paeall = np.array(metadata_data['pae'])
            pae = paeall[0:totalproteinlength, 0:totalproteinlength]        
            if pae.shape == (0, 0): # if pae is empty, which means no proteins in prediction
                pae = 1
                   
            # PAE calculations
            PAEminprot = np.amin(pae)
            PAEmeanprot = np.mean(pae)
            PAEstdprot = np.std(pae)
              
            scores[index, 11] = 30-PAEminprot
            scores[index, 12] = PAEstdprot/PAEmeanprot
            
            #extract pairwise pae
            paechainall = np.empty([numspecies, numspecies], dtype=object)
            paechainmin = np.empty([numspecies, numspecies], dtype=object)
            paechainstd = np.empty([numspecies, numspecies], dtype=object)
            xstart=0
            for position1 in range(0, numspecies):
                ystart=0
                for position2 in range(0, numspecies):
                    paequad = paeall[xstart:xstart + specieslength[position1], ystart:ystart + specieslength[position2]]
                    paechainall[position1, position2] = paequad.flatten()
                    
                    PAEmin = np.amin(paequad)
                    PAEmean = np.mean(paequad)
                    PAEstd = np.std(paequad)
                    paechainmin[position1, position2] = 30-PAEmin
                    paechainstd[position1, position2] = PAEstd/PAEmean
                    
                    ystart = ystart + specieslength[position2]
                xstart = xstart + specieslength[position1]
                
            listmin = []
            liststd = []
            for species in range(0, numspecies):
                listmin.append(paechainmin[species,:].tolist())
                liststd.append(paechainstd[species,:].tolist())
            scores[index, 13] = str(listmin)
            scores[index, 14] = str(liststd)
        
        
            #initialize empty plot
            plt.clf()
            
            # Plot a pLDDT line graph
            plt.plot(range(0, len(plddtall)), plddtall, label=f'rank_{rank_num_pad}')
        
            # Add labels and title
            plt.xlabel('Atom')
            plt.ylabel('Predicted LDDT')
            plt.title('Predicted LDDT per atom')
            
            plt.legend()    # Add a legend   
            plt.savefig(os.path.join(dirfilessave, f'{scores[index, 1]}_rank_{rank_num_pad}_plddt.png'))   # Save the graph
            plt.close() #remove this to display graph
            
        
            # Plot pae graph
            cmap = plt.colormaps['coolwarm']
            fig, ax = plt.subplots(1, 1, figsize=(4, 4)) # Create a figure and axis objects
            
            # Plot each heat map
            im = ax.imshow(paeall, cmap=cmap, vmin=0, vmax=30)
            ax.set_title(f'rank_{rank_num_pad}', fontsize=18)
            ax.axis('on')  # Hide axis labels and ticks
        
            # Add a colorbar
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
            # Save    
            plt.savefig(os.path.join(dirfilessave, f'{scores[index, 1]}_rank_{rank_num_pad}_pae.png'))   # Save the graph
            plt.close() #remove this to display graph
        
            
            # plot pairwise ipTM
            cmap = cmap = plt.colormaps['Greys_r']
            fig, ax = plt.subplots(1, 1, figsize=(4, 4)) # Create a figure and axis objects

            # Plot each heat map
            ipTMgraphsize = iptmpair.shape[0]
            im = ax.imshow(iptmpair[0:ipTMgraphsize, 0:ipTMgraphsize], cmap=cmap, vmin=0, vmax=1)
            ax.set_title(f'rank_{rank_num_pad}', fontsize=18)
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.axis('on')  # Hide axis labels and ticks

            # Annotate each cell with the numeric value
            for i in range(0, ipTMgraphsize):
                for j in range(0, ipTMgraphsize):
                    value = iptmpair[i, j]
                    if value <= 0.5:
                        ax.text(j, i, f'{value:.2f}', ha='center', va='center', fontsize=8, color="white")
                    else:
                        ax.text(j, i, f'{value:.2f}', ha='center', va='center', fontsize=8, color="black")
            
            # Add a colorbar
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
       
            # Save
            plt.savefig(os.path.join(dirfilessave, f'{scores[index, 1]}_rank_{rank_num_pad}_ipTMpair.png'))   # Save the graph
            plt.savefig(os.path.join(dirfilessave, f'{scores[index, 1]}_rank_{rank_num_pad}_ipTMpair.svg'), format='svg')   # Save the graph
            plt.close() #remove this to display graphs

            index += 1
    else:
        print(f'   no predictions in folder: {listdir[folder_num]}')     
scoressave = np.delete(scores, range(index, scores.shape[0]), 0)


# save scoressave file
scoresDF = pd.DataFrame(scoressave)
with pd.ExcelWriter(filepath) as writer:
    scoresDF.to_excel(writer, header=["file original", "file new", "fraction_disordered", "has_clash", "global ipTM", "chain_pair_ipTM", "global pTM", "chain_pair_pTM", "ranking_score",\
                                      "global protein mean pLDDT", "global protein stdev pLDDT", \
                                      "global protein 30-PAEmin", "global protein PAEstd/PAEmean", \
                                      "chain_pair 30-PAEmin", "chain_pair PAEstd/PAEmean"])

# end
print("\nproteins found:  ", len(listdir))
print('done calculating and saved in: \n'+filepath)
