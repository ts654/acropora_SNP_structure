# -*- coding: utf-8 -*-
"""
Description:
    Take pathway lists for all conditions and merge into 1 array and generate a heatmap.
    Uses the output from pathway_analysis.py

input:
    pathway lists (keyword or GO) for multiple conditions
    
output:
    .tsv of merged list (absolute counts and relative total)
    .png and .svg of heatmap

"""

# Dependencies
import os
import numpy as np
from matplotlib import pyplot as plt


# Constants
CONDITIONS = ['bl', 'fvfm', 'larvsurv', 'surv']
NUMINPUT = len(CONDITIONS)
NUMSHOW = 30

dir_parent = r'path\to\folder\of\conditions'
path_files_keywords = np.empty(NUMINPUT, dtype=object)
path_files_GO = np.empty(NUMINPUT, dtype=object)
for i in range(0, NUMINPUT):
    path_files_keywords[i] = os.path.join(dir_parent, f'aten-{CONDITIONS[i]}_freq_8_keyword.tsv')
    path_files_GO[i] = os.path.join(dir_parent, f'aten-{CONDITIONS[i]}_freq_GO.tsv')


# Functions
def MergeLists(list_input_words, category):
    # merge lists into a dictionary
    dict_words_count = {}
    count_words = np.zeros(NUMINPUT)
    for i in range(0, NUMINPUT):
        for j in range(0, list_input_words[i].shape[0]):
            count_words[i] = count_words[i] + int(list_input_words[i][j, 1])
            if not list_input_words[i][j, 0] in dict_words_count: #if not in list, create empty entry
                dict_words_count[list_input_words[i][j, 0]] = np.zeros(NUMINPUT)
            dict_words_count[list_input_words[i][j, 0]][i] = list_input_words[i][j, 1]        

    # convert dictionary to np of counts and save
    list_words_count = np.empty([len(dict_words_count.keys()), 1 + NUMINPUT + 1], dtype='object')
    for i, word in enumerate(dict_words_count):
        list_words_count[i, 0] = word
        for j in range(0, NUMINPUT):
            list_words_count[i, j+1] = int(dict_words_count[word][j])
        list_words_count[i, NUMINPUT + 1] = max(list_words_count[i, 1:NUMINPUT + 1])    
    list_words_count_sorted = np.asarray(sorted(list_words_count, key=lambda x: x[NUMINPUT + 1], reverse = True))
    list_words_count_sorted_header = np.insert(list_words_count_sorted, 0, [category] + CONDITIONS + ['max in row'], axis=0)
    np.savetxt(os.path.join(dir_parent, f'merge_list_{category}_count.tsv'), list_words_count_sorted_header, fmt="%5s", delimiter="\t")

    # convert dictionary to np of percentages and save 
    list_words_percent = np.empty([len(dict_words_count.keys()), 1 + NUMINPUT + 1], dtype='object')
    for i, word in enumerate(dict_words_count):
        list_words_percent[i, 0] = word
        for j in range(0, NUMINPUT):
            list_words_percent[i, j+1] = int(dict_words_count[word][j])/count_words[j]
        list_words_percent[i, NUMINPUT + 1] = max(list_words_percent[i, 1:NUMINPUT + 1])    
    list_words_percent_sorted = np.asarray(sorted(list_words_percent, key=lambda x: x[NUMINPUT + 1], reverse = True))
    list_words_percent_sorted_header = np.insert(list_words_percent_sorted, 0, [category] + CONDITIONS + ['max in row'], axis=0)
    np.savetxt(os.path.join(dir_parent, f'merge_list_{category}_percent.tsv'), list_words_percent_sorted_header, fmt="%5s", delimiter="\t")
        
    # if GO_term, remove GO:### from front of label
    if category == 'GO_terms':
        for i in range(0, list_words_percent_sorted.shape[0]):
            if not list_words_percent_sorted[i, 0] == '':
                list_words_percent_sorted[i, 0] = list_words_percent_sorted[i, 0].split(' ', 1)[1]
    
    return list_words_percent_sorted

def MakeHeatmap(list_merge_words_percent, category):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4)) # Create a figure and axis objects
    data = list_merge_words_percent[0:NUMSHOW, 1:NUMINPUT + 1].astype(float)
    cmap = plt.colormaps['viridis']
    
    im = ax.imshow(data, cmap=cmap, vmin=0, vmax=np.amax(data))
    ax.set_title(f'heatmap: {category}', fontsize=18)
    plt.xticks(np.arange(data.shape[1]), labels=CONDITIONS, rotation=90, ha='center')
    plt.yticks(np.arange(data.shape[0]), labels=list_merge_words_percent[0:NUMSHOW, 0])
    ax.axis('on')  # Hide axis labels and ticks
    
    # Add a colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
       
    # Save
    plt.savefig(os.path.join(dir_parent, f'heatmap_{category}.png'))                 # Save the graph
    plt.savefig(os.path.join(dir_parent, f'heatmap_{category}.svg'), format='svg')   # Save the graph in vector format
    

# Main
#load files
list_input_keywords = []
list_input_GO = []
for i in range(0, NUMINPUT):
    list_input_keywords.append(np.genfromtxt(path_files_keywords[i], delimiter='\t', dtype='str'))
    list_input_GO.append(np.genfromtxt(path_files_GO[i], delimiter='\t', dtype='str'))

#merge lists and save
list_merge_keywords_percent = MergeLists(list_input_keywords, 'keywords')
list_merge_GO_percent = MergeLists(list_input_GO, 'GO_terms')

#make heatmaps and save
MakeHeatmap(list_merge_keywords_percent, 'keywords')
MakeHeatmap(list_merge_GO_percent, 'GO_terms')

print('done')