# for all mols, generate mol from SMILES, then generate 2048 bit morganFP using rdkit implementation
# write a csv file, each row is one experiment, each column is one bit of the morgan fingerprint
# df size: (# of molecules, n_bits)

import numpy as np
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Chem import DataStructs


def fingerprint_csv(mols, n_bits=2048, radius=2, output_path=''):

    mol_df = pd.DataFrame(mols, columns=['SMILES'])
    PandasTools.AddMoleculeColumnToFrame(mol_df, smilesCol='SMILES', molCol='ROMol')
    assert mol_df.isnull().sum().sum() == 0, 'some rdkit mol files fail to generate'

    # featurize with morgan FP
    mol_df['morganFP'] = mol_df.apply(lambda x: GetMorganFingerprintAsBitVect(x['ROMol'], radius=radius, nBits=n_bits), axis=1)
    mol_df = mol_df.drop(['ROMol'], axis=1)  # faster lookup
    mol_df = mol_df.set_index('SMILES')  # use SMILES as df index

    cols = ['']*n_bits
    df = pd.DataFrame(columns=cols, index=mol_df.index)
    for index, row in mol_df.iterrows():  # not ideal, but only run it once to create full set, okay
        fp = np.zeros((0,))
        DataStructs.ConvertToNumpyArray(row['morganFP'], fp)
        df.loc[index] = list(fp)
    assert df.isnull().sum().sum() == 0

    # save to csv
    if output_path is not None:
        df.to_csv(output_path)  # with index (name)


if __name__ == '__main__':
    df = pd.read_csv('https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/deebo/ami.csv')
    fingerprint_csv(df['nucleophile_smiles'].unique(), output_path='test.csv')
