# -*- coding: utf-8 -*-
"""
adapted from Saskia Sanders' file from publication
    Keyword_Frequency_SP_PDB.ipynb

    Soh TK, Ognibene S, Sanders S, Schäper R, Kaufer BB, Bosse JB. A proteome-wide structural systems approach reveals insights into protein families of all human herpesviruses.
    Nat Commun. 2024 Nov 26;15(1):10230. doi: 10.1038/s41467-024-54668-2. PMID: 39592652; PMCID: PMC11599850.

Description:
    The Word Frequency Analysis Tool is a Python script to analyze text data and generate frequency statistics for words within specific domains. The tool is primarily intended for processing TSV (Tab-Separated Values) files containing domain and word information.
    This Python script processes a TSV file containing information about virus proteins and their associated words. It extracts word frequency data for each virus protein domain, excluding specified common words. The resulting frequency data is then saved to a new TSV file for further analysis.
    Exclusion of Common and Irrelevant Words: The script allows for the exclusion of a predefined list of common English words, Greek alphabet characters, and specific user-defined words. These exclusions help focus the analysis on more meaningful content.
    Case Insensitivity and Punctuation Handling: The script processes text in a case-insensitive manner, ensuring that words are counted regardless of their capitalization. Additionally, it handles punctuation marks, such as commas and double quotes, ensuring they do not interfere with word recognition.
    Exclusion of Single Letters and Numbers: The tool automatically excludes words consisting of only a single letter or a single number, ensuring that the analysis focuses on meaningful terms.
    Excluded Words: The list of words to be excluded can be customized by modifying the exclude_words list within the script
    Example input: HCMV_UL18_domain_0.pdb Crystal structure of HLA DR52c 3c5j_B 7.74E-10 HCMV_UL18_domain_0.pdb Structure of chicken CD1-2 with bound fatty acid 3dbx_A 4.10E-11 HCMV_UL18_domain_0.pdb TCR complex 2wbj_F 7.32E-10 HCMV_UL18_domain_0.pdb NEONATAL FC RECEPTOR, PH 6.5" 3fru_C 6.58E-12 HCMV_UL18_domain_0.pdb HLA-DR1 with GMF Influenza PB1 Peptide 6qza_BBB 1.14E-09 HCMV_UL18_domain_0.pdb Crystal structure of FcRn bound to UCB-84 6c98_C 6.23E-12 HCMV_UL18_domain_0.pdb immune receptor complex 6v15_B 8.18E-10 HCMV_UL18_domain_0.pdb Immune Receptor 4mdi_B 3.02E-10 [...]

input:
    foldseek_sp_hits.tsv
        file that contains only the significant hits from the foldseek output

output:
    table containing: protein, top hit, foldseek keyword frequency, GO annotation of top hit
    
    foldseek_sp_freq.tsv   # list of description words excluding common words
    foldseek_sp_freq-1.tsv # list of description words excluding common words and words that appear only once
"""

# Dependencies

import os
import requests


# Define the file path manually
dir_parent = r'path\to\dir'
path_sp = os.path.join(dir_parent, 'foldseek_sp_hits.tsv')

# Add your list of excluded words here
exclude_words = ["and", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota",
                 "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon",
                 "phi", "chi", "psi", "omega", "a", "an", "the", "in", "on", "to", "with", "for", "of",
                 "as", "by", "at", "ii", "cell", "protein", "probable", "motif", "chain", "type", "molecule",
                 "uncharacterized", "c-x-c", "c-c", "basic", "structure", "complex", "crystal", "bound", "fab", "domain",
                 "antibody", "mutant", "peptide", "cryo-em", "fragment", "complexed", "the", "compound", "atomic", "structural",
                 "by", "coli", "cryoem", "resolution", "its" , "form", "e.", "x-ray", "class", "structures", "-", "class",
                 "therapeutic", "deletion", "variant", "site", "functional", "reveal", "c1", "implications", "monoclonal",
                 "different", "allosteric", "ligand", "c-terminal", "low-ph", "escherichia", "e.coli", "between", "complexes.",
                 "(crystal", "engineered", "full", "xfel", "natural", "angstroms", "allosteric", "design", "discovery", "open"
                 "edition", "mode", "opposite", "novel", "region", "ligand-binding", "apo", "holoenzyme", "high", "full-length",
                 "loop", "angstrom", "template-primer", "based", "using", "from", "open", "conserved", "significance", "unit",
                 "loader", "neutralizing", "specificity", "substrate", "one", "inhibitor", "inhibitors", "after", "dehydration",
                 "de", "novo", "symmetry", "element", "iii", "ion", "", "herpes", "virus", "herpesvirus", 
                 "subunit", "factor", "homolog"]


# Functions

def is_valid_word(word):
    return not (len(word) == 1 and (word.isalpha() or word.isdigit())) and not (len(word) == 2 and word.isdigit())


def get_go_terms(uniprot_id):
    # Get UniProt entry with GO annotations
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    data = requests.get(uniprot_url).json()

    go_terms = []
    for xref in data.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "GO":
            go_terms.append({
                "id": xref["id"],
                "term": xref.get("properties", [{}])[0].get("value", "")
            })
            
    # Sanitize and format each GO term
    clean_terms = [f"{term['id']} {term['term'].replace(',', '')}" for term in go_terms]
    # Join into a single comma-separated string
    go_string = ', '.join(clean_terms)

    return go_string


def extract_word_frequency(filename, exclude_words, exclude_single_count=True):
    domain_word_freq = {}
    domain_top_hit = {}
    domain_go = {}
    protein_previous = ''
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            virus_protein_domain_parts = parts[0].split('_')
            virus_protein_domain = '_'.join(virus_protein_domain_parts).replace('.pdb', '')  # Take all parts and exclude ".pdb" extension
            words = parts[2].split()

            if virus_protein_domain not in domain_word_freq:
                domain_word_freq[virus_protein_domain] = {}
                domain_top_hit[virus_protein_domain] = {}
                domain_go[virus_protein_domain] = {}

            if not parts[0] == protein_previous: # if at new protein, i.e. most significant hit since it is at the top
                # save top hit
                domain_top_hit[virus_protein_domain] = parts[2]
                
                # once have top hit, get GO terms for that protein based on uniprot
                uniprot_id = parts[2].split('-')[1]
                domain_go[virus_protein_domain] = get_go_terms(uniprot_id)
            protein_previous = parts[0]

            for word in words:
                # Convert both the word and exclude_words to lowercase for case insensitivity
                word = word.lower().replace(',', '').replace('"', '').replace('(',"").replace(')',"").replace(':', "")  # Remove commas, double quotesparentheses and quotes from word
                if word not in exclude_words and is_valid_word(word):
                    domain_word_freq[virus_protein_domain][word] = domain_word_freq[virus_protein_domain].get(word, 0) + 1

    # Exclude words with only one count if exclude_single_count is True
    if exclude_single_count:
        for domain in domain_word_freq:
            domain_word_freq[domain] = {word: count for word, count in domain_word_freq[domain].items() if count > 1}

    return domain_word_freq, domain_top_hit, domain_go


def save_frequency_tsv(input_path, domain_word_freq, domain_top_hit, domain_go, suffix="_freq"):
    output_folder, input_filename = os.path.split(input_path)
    output_filename = os.path.splitext(input_filename)[0] + f"{suffix}.tsv"
    output_path = os.path.join(output_folder, output_filename)

    with open(output_path, 'w') as output_file:
        sorted_domains = sorted(domain_word_freq.keys())  # Sort domains alphabetically

        for domain in sorted_domains:
            word_freq = domain_word_freq[domain]
            sorted_freq = sorted(word_freq.items(), key=lambda x: x[1], reverse=True)
            result = [f'{key}({value})' for key, value in sorted_freq]
            output_file.write(f'{domain}\t{", ".join(result)}\t{domain_top_hit[domain]}\t{domain_go[domain]}\n')
            
    print(f'Frequency data saved in file: {output_filename}')


def analyze_file(_file_path):
    # Save with "_freq-1" suffix (exclude single count)
    domain_word_freq, domain_top_hit, domain_go = extract_word_frequency(_file_path, exclude_words)
    save_frequency_tsv(_file_path, domain_word_freq, domain_top_hit, domain_go, suffix="_freq-1")

    # Save with "_freq" suffix (do not exclude single count)
    domain_word_freq_no_exclusion, domain_top_hit, domain_go = extract_word_frequency(_file_path, exclude_words, exclude_single_count=False)
    save_frequency_tsv(_file_path, domain_word_freq_no_exclusion, domain_top_hit, domain_go, suffix="_freq")


# Main
if __name__ == "__main__":
    print(f'analyzing:\n{path_sp}\n')
    analyze_file(path_sp)
