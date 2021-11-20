# lookup molecules in kraken

import requests

import yaml
import pandas as pd
import numpy as np
from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors


# search priority: cas -> name -> canonical smiles -> atom count and formula -> inchi -> keywords
# This search method displays result as a dataframe. User selection is still required
# Note 1: remove salt before querying with smiles

def lookup(cas=None, name=None, smiles=None, inchi=None, keywords=None):

    identifiers = pd.read_csv('identifiers.csv')
    identifiers['id'] = identifiers['id'].astype('int')
    identifiers['atomcount'] = identifiers['atomcount'].astype('int')
    results_df = pd.DataFrame()

    # inner function: serach by smiles
    def search_by_smiles(smi):

        can_smiles_list = list(identifiers['can_smiles'].astype('str'))
        can_smiles = Chem.CanonSmiles(smi)

        results = [c for c in can_smiles_list if can_smiles == c]

        if results:  # find smiles match, get results
            temp_df = identifiers[identifiers['can_smiles'].isin(results)]
        else:  # no smiles match, try to match with mw, atom count, formula
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            atom_count = mol.GetNumAtoms()
            formula = rdMolDescriptors.CalcMolFormula(mol)

            temp_df = identifiers.loc[
                (identifiers['atomcount'] == int(atom_count)) &
                (identifiers['formula'] == str(formula))
                ]

        return temp_df

    if cas is not None:

        if type(cas) != str:
            raise TypeError('Search value must be a string!')

        cas_list1 = list(identifiers['cas_pr3'].astype('str'))
        results_1 = [c for c in cas_list1 if cas in c]
        if results_1:
            df_1 = identifiers[identifiers['cas_pr3'].isin(results_1)]
        else:
            df_1 = pd.DataFrame()

        cas_list2 = list(identifiers['cas_pr3_hx'].astype('str'))
        results_2 = [c for c in cas_list2 if cas in c]
        if results_2:
            df_2 = identifiers[identifiers['cas_pr3_hx'].isin(results_2)]
        else:
            df_2 = pd.DataFrame()

        cas_list3 = list(identifiers['cas_conj_acid'].astype('str'))
        results_3 = [c for c in cas_list3 if cas in c]
        if results_3:
            df_3 = identifiers[identifiers['cas_conj_acid'].isin(results_3)]
        else:
            df_3 = pd.DataFrame()

        results_df = pd.concat([df_1, df_2, df_3], axis=1)

        if len(results_df) == 1:  # only find one entry, return search result
            print('Found result by CAS number')
            return results_df[['ligand', 'id', 'can_smiles']]

    if name is not None:

        if type(name) != str:
            raise TypeError('Search value must be a string!')

        name_list = list(identifiers['ligand'].astype('str'))

        # exact match for name
        results = [n for n in name_list if name.lower() == n.lower()]

        # no exact match, search for partial match
        if len(results) == 0:
            results = [n for n in name_list if name.lower() in n.lower()]

        if results:
            temp_df = identifiers[identifiers['ligand'].isin(results)]
            results_df = pd.concat([results_df, temp_df])
            if len(results_df) == 1:  # only find one entry, return search result
                print('Found result by ligand name')
                return results_df[['ligand', 'id', 'can_smiles']]

    if smiles is not None:

        if type(smiles) != str:
            raise TypeError('Search value must be a string!')

        temp_df = search_by_smiles(smiles)
        results_df = pd.concat([results_df, temp_df])

        if len(results_df) == 1:  # only find one entry, return search result
            print('Found result by smiles')
            return results_df[['ligand', 'id', 'can_smiles']]

    if inchi is not None:

        if type(inchi) != str:
            raise TypeError('Search value must be a string!')

        mol = Chem.MolFromInchi(inchi)
        smi = Chem.MolToSmiles(mol)
        temp_df = search_by_smiles(smi)
        results_df = pd.concat([results_df, temp_df])

        if len(results_df) == 1:  # only find one entry, return search result
            print('Found result by inchi')
            return results_df[['ligand', 'id', 'can_smiles']]

    if keywords is not None:

        if type(keywords) != list:
            message = 'Search value for keywords must be a list of strings! ' \
                      'Example: lookup_kraken(keywords=[\'CF3\', \'Ph\'])'
            raise TypeError(message)

        boo = np.zeros(shape=(len(identifiers), len(keywords)))
        for i in range(len(keywords)):
            boo[:, i] = identifiers['ligand'].str.contains(keywords[i], na=False, case=False)
        boo = np.all(boo.astype(bool), axis=1)  # name must contain all keywords

        if np.sum(boo):
            temp_df = identifiers.loc[boo]
            results_df = pd.concat([results_df, temp_df])
            if len(results_df) == 1:  # only find one entry, return search result
                print('Found result by keywords')
                return results_df[['ligand', 'id', 'can_smiles']]

    if len(results_df) == 0:
        print('No matches found. This ligand is possibly not in kraken.')
        return None

    print('No exact match found, here are some possibilities')
    return results_df[['ligand', 'id', 'can_smiles']]


def access(id, mode=None):
    """
    access raw kraken data with kraken id and mode
    returns loaded data and mode

    :param id: kraken id as integer
    :type id: int
    :param mode: confdata, data, energy
    :type mode: str
    :return: loaded data, mode
    :rtype:
    """

    github_url = 'https://raw.githubusercontent.com/doyle-lab-ucla/krkn-raw/main/raw/'

    # add leading zeros (kraken currently uses 8 digits)
    k_id = str(id).zfill(8)

    if mode == 'confdata':
        file_name = k_id + '_confdata.yml'
    elif mode == 'data':
        file_name = k_id + '_data.yml'
    elif mode == 'energy':
        file_name = k_id + '_relative_energies.csv'
    else:
        raise ValueError('unknown mode; choose between confdata, data, energy')

    url = str(github_url + file_name)

    if mode == 'energy':
        data = pd.read_csv(url, delimiter=';')
    else:
        headers = {'Accept': 'text/yaml'}
        r = requests.get(url, headers=headers)
        if r.status_code == 200:
            data = yaml.load(r.text, Loader=yaml.Loader)
        else:
            raise ConnectionError('cannot retrieve data from Github (likely internet issue), try again')

    # return mode to differentiate between formats
    return data, mode


def find():
    """lookup and return only the id for one ligand (the only match, the best match, no match)"""
    return None


def access_vburmin_conf():
    return None


# for testing only
if __name__ == '__main__':
    data = access(1583, mode='confdata')
    for d in data:
        print(d)
        exit()
    # print(lookup(keywords=['Ph', 'tBu', 'OMe', 'jason']))