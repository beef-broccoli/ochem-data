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
  
  # change file path accordingly
  PATH = 'https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/deoxyF/paper-dft/train.csv'
  df = pd.read_csv(PATH)
  
  # do things with df...
  ```
  
Each subfolder (`ohe`, `mol2vec`, `mordred`...) includes different descriptor encodings for the list of molecules or reaction entries

## Publications for reaction data: 

- deoxyF: [Deoxyfluorination with Sulfonyl Fluorides: Navigating Reaction Space with Machine Learning](https://pubs.acs.org/doi/10.1021/jacs.8b01523)
- Vbur: [Linking Mechanistic Analysis of Catalytic Reactivity Cliffs to Ligand Classification](https://chemrxiv.org/engage/chemrxiv/article-details/60c758aabdbb89d828a3ade9)
- CN: [Predicting reaction performance in C–N cross-coupling using machine learning](https://www.science.org/doi/full/10.1126/science.aar5169?versioned=true)

## Publications for descriptor sets: 

- [Mordred: a molecular descriptor calculator](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0258-y)
- [Mol2vec: Unsupervised Machine Learning Approach with Chemical Intuition](https://pubs.acs.org/doi/full/10.1021/acs.jcim.7b00616)
- [alvaDesc](https://www.alvascience.com/alvadesc/)
- Kraken: [A Comprehensive Discovery Platform for Organophosphorus Ligands for Catalysis](https://chemrxiv.org/engage/chemrxiv/article-details/60c757f9702a9bdb7018cbd4)


