# ochem-data
This repo compiles published organic chemistry data, including raw data and calculated descriptor sets for molecules and reactions. 

## File structure: 

Each folder contains one reaction (see publication list for details). 

Every reaction folder was divided into two sub-folders: 
- `mols`: a list of molecules categorized by reaction roles)
- `rxns`: a list of reaction entries with multiple reaction components and yield). 

## Reading data directly from this repo

  ```Py
  import pandas as pd
  
  REPO_PATH = 'https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/'
  FP = 'deoxyF/paper-dft/train.csv'  # change this  
  df = pd.read_csv(REPO_PATH + FP)
  
  # do things with df...
  ```
  
Each subfolder (`ohe`, `mol2vec`, `mordred`...) includes different descriptor encodings for the list of molecules or reaction entries

## Publications for reaction data: 

- deoxyF: [Deoxyfluorination with Sulfonyl Fluorides: Navigating Reaction Space with Machine Learning](https://pubs.acs.org/doi/10.1021/jacs.8b01523)

<p align="center">
  <img src=https://github.com/beef-broccoli/misc-files/blob/f370f89ab3370839df18bd79fa1f89a091ada630/ochem-data/deoxyF_scheme.png alt=deoxyF width=600/>
</p>

- Vbur: [Linking Mechanistic Analysis of Catalytic Reactivity Cliffs to Ligand Classification](https://chemrxiv.org/engage/chemrxiv/article-details/60c758aabdbb89d828a3ade9)

<p align="center">
  <img src=https://github.com/beef-broccoli/misc-files/blob/cde79f8752feddfe23e346591c10aae826d403ac/ochem-data/Vbur_scheme.png alt=Vbur width=600/>
</p>

- CN: [Predicting reaction performance in C–N cross-coupling using machine learning](https://www.science.org/doi/full/10.1126/science.aar5169?versioned=true)

<p align="center">
  <img src=https://github.com/beef-broccoli/misc-files/blob/cde79f8752feddfe23e346591c10aae826d403ac/ochem-data/CN_scheme.png alt=CN width=600/>
</p>

- asym_epox: [Ni/Photoredox-Catalyzed Enantioselective Cross-Electrophile Coupling of Styrene Oxides with Aryl Iodides](https://chemrxiv.org/engage/chemrxiv/article-details/60c7574abb8c1a7a333dc7d0)

<p align="center">
  <img src=https://github.com/beef-broccoli/misc-files/blob/9b3b3ef4a612c0aa232c0f5b138fa492823aa555/ochem-data/asym_epox_scheme.png alt=asym_epox width=600/>
</p>

- CH-aryl-1: [Bayesian reaction optimization as a tool for chemical synthesis](https://www.nature.com/articles/s41586-021-03213-y)

<p align="center">
  <img src=https://github.com/beef-broccoli/misc-files/blob/cde79f8752feddfe23e346591c10aae826d403ac/ochem-data/CH-aryl-1_scheme.png alt=CHaryl-1 width=600/>
</p>

- Kraken: [A Comprehensive Discovery Platform for Organophosphorus Ligands for Catalysis](https://chemrxiv.org/engage/chemrxiv/article-details/60c757f9702a9bdb7018cbd4)

<p align="center">
  <img src=https://github.com/beef-broccoli/misc-files/blob/60475f1bcbec6d5c7d400e5f5dab6f06eb993aea/ochem-data/kraken_scheme.png alt=kraken width=600/>
</p>

Supports a search function for kraken ID with cas number, name or smiles. Check out the demo here: 

## Publications/resources for descriptor sets: 

- [Mordred: a molecular descriptor calculator](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0258-y)
- [Mol2vec: Unsupervised Machine Learning Approach with Chemical Intuition](https://pubs.acs.org/doi/full/10.1021/acs.jcim.7b00616)
- [alvaDesc](https://www.alvascience.com/alvadesc/)
- Kraken: [A Comprehensive Discovery Platform for Organophosphorus Ligands for Catalysis](https://chemrxiv.org/engage/chemrxiv/article-details/60c757f9702a9bdb7018cbd4)
- [auto-qchem](https://github.com/PrincetonUniversity/auto-qchem)
