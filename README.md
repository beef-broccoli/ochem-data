# ochem-data
This repo compiles published organic chemistry data, including raw data and calculated descriptor sets for molecules and reactions. 

## Publication list: 

- deoxyF: [Deoxyfluorination with Sulfonyl Fluorides: Navigating Reaction Space with Machine Learning](https://pubs.acs.org/doi/10.1021/jacs.8b01523)

- Vbur: [Linking Mechanistic Analysis of Catalytic Reactivity Cliffs to Ligand Classification](https://chemrxiv.org/engage/chemrxiv/article-details/60c758aabdbb89d828a3ade9)

## Reading data directly from this repo

  `import pandas as pd`
  
  `PATH = 'https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/deoxyF/paper-dft/train.csv'`
  
  `df = pd.read_csv(PATH)`
  
