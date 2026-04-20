# acropora_SNP_structure

Custom Python scripts for the structure-guided analysis of GWAS-derived SNPs in Acropora sp. heat tolerance, as described in:

**Protein structure prediction identifies stability-modulating mutations alongside signaling and immunity pathways as drivers of heat tolerance in reef-building Acropora corals
Timothy K. Soh, David Edwards, Jens B. Bosse, Kate M. Quigley**

The scripts map heat-tolerance-associated SNPs from ANGSD GWAS output onto annotated coral genomes, classify variants (intergenic / intronic / synonymous / non-synonymous), prepare reference and mutant sequences for AlphaFold structure prediction, extract AlphaFold 3 confidence scores, and perform keyword/GO-based pathway analysis of Foldseek structure-similarity hits.

---

## Overview

```
  ANGSD GWAS output (significant SNPs per trait)
                    │
                    ▼
      map_SNP_to_proteins_{aten,amil}.py
                    │
      ┌─────────────┼──────────────┐
      ▼             ▼              ▼
  intergenic     intron/          CDS
                 synonymous       (non-synonymous)
      │             │              │
      ▼             ▼              ▼
  reference     reference       reference + mutant
  sequences     sequences       sequences
      │             │              │
      └──────┬──────┘              │
             ▼                     ▼
     AlphaFold 3 / LocalColabFold  AlphaFold 3 (100 seeds)
             │                     │
             ▼                     ▼
         score_AF3.py          score_AF3.py
             │                     │
             ▼                     ▼
         Foldseek              RaSP (ΔΔG)
             │                     │
             ▼                     ▼
  keyword_frequency_sp_with_GO.py  (stability analysis)
             │
             ▼
      pathway_analysis.py
             │
             ▼
  pathway_analysis_merged_heatmap.py
```

## Dependencies:
	python		3.9
	bio			1.6.2
	matplotlib	3.5.3
	numpy		1.21.6
	pandas		1.3.5
	requests	2.32.2

