import requests
from io import StringIO
import multiprocessing
import itertools
import timeit
import yaml
import json
import pandas as pd
import numpy as np
from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
from urllib.error import HTTPError

pd.set_option("display.max_rows", None, "display.max_columns", None, 'display.max_colwidth', None)


def lookup(cas=None, name=None, smiles=None, inchi=None, keywords=None, verbose=1):
    """
    lookup a single ligand in kraken
    Input: cas, name, smiles, inchi or keywords
    Output: phosphines in kraken database with its ID, name and smiles
    Notes:
        - remove salt before querying with smiles
        - search priority: cas -> name -> canonical smiles -> atom count and formula -> inchi -> keywords
        - This search method displays result as a dataframe. User selection is still required

    :param cas:
    :type cas: str
    :param name:
    :type name: str
    :param smiles:
    :type smiles: str
    :param inchi:
    :type inchi: str
    :param keywords:
    :type keywords:
    :param verbose:
    :type verbose:
    :return: search results
    :rtype: pd.DataFrame
    """

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
            if verbose:
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
                if verbose:
                    print('Found result by ligand name')
                return results_df[['ligand', 'id', 'can_smiles']]

    if smiles is not None:

        if type(smiles) != str:
            raise TypeError('Search value must be a string!')

        temp_df = search_by_smiles(smiles)
        results_df = pd.concat([results_df, temp_df])

        if len(results_df) == 1:  # only find one entry, return search result
            if verbose:
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
            if verbose:
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
                if verbose:
                    print('Found result by keywords')
                return results_df[['ligand', 'id', 'can_smiles']]

    if len(results_df) == 0:
        if verbose:
            print('No matches found. This ligand is possibly not in kraken.')
        return None

    if verbose:
        print('No exact match found, here are some possibilities')
    return results_df[['ligand', 'id', 'can_smiles']]


def lookup_multi(values, identifier):
    """
    Builds on top of lookup()
    Input: a list of identifier values, type of identifiers
    Returns: if found one result, return kraken ID. If multiple or no results are found, outputs possibilities and prompt
    user to enter the correct kraken ID for ligand

    :param values:
    :param identifier:
    :return:
    """

    if identifier not in ['cas', 'name', 'smiles', 'inchi']:
        raise ValueError('identifier needs to be cas, name, smiles or inchi. (keyword is not supported)')

    results = []
    for val in values:
        d = {identifier: val,
             'verbose': 0}
        lookup_result = lookup(**d)
        if len(lookup_result) == 1:
            results.append(lookup_result['id'].values[0])
        else:  # no single match, return original identifier value
            print('For ligand with \"{0}\" \"{1}\", no match is found. Here are some possibilities:'.format(identifier, val))
            print(lookup_result)
            input_id = int(input('Please input a kraken ID for this ligand (0 if none matches): '))
            if input_id:
                results.append(input_id)
            else:
                results.append(val)

    return results


def access(k_id, mode=None, verbose=1):
    """
    For a single ligand, access raw kraken data with kraken id (k_id) and mode (data, confdata, energy).
    Note: Use access_multi() for batch operation. For loop is slow

    :param k_id: kraken id
    :type k_id: int
    :param mode: 'confdata', 'data', 'energy'
    :type mode: str
    :param verbose: verbose
    :type verbose: int
    :return: loaded data
    :rtype: Union[pd.DataFrame, dict]
    """

    github_url = 'https://raw.githubusercontent.com/doyle-lab-ucla/krkn/main/raw/'

    if mode == 'confdata':
        suffix = '_confdata.yml'
    elif mode == 'data':
        suffix = '_data.yml'
    elif mode == 'energy':
        suffix = '_relative_energies.csv'
    else:
        raise ValueError('unknown mode; choose between confdata, data, energy')

    # add leading zeros (kraken currently uses 8 digits)
    full_id = str(k_id).zfill(8)

    r = requests.get(github_url + full_id + suffix)
    if r.status_code == 404:
        if verbose:
            print('ligand {0} not found'.format(full_id))
        return None
    elif r.status_code == 200:
        if mode == 'energy':
            return pd.read_csv(StringIO(r.text), delimiter=';')
        else:
            return yaml.load(r.text, Loader=yaml.Loader)
    else:
        raise ConnectionError('unknown networking issue, status code {0}'.format(r.status_code))


def access_multi(ids, mode=None, howmanycores=8):
    """
    Access the data for a list ligands. Use multiprocessing to speed up access()

    :param ids: list of kraken ids
    :type ids: list of int
    :param mode: 'confdata', 'data', 'energy'
    :type mode: str
    :param howmanycores: number of cores of this computer (decides how many parallel processes can be initiated)
    :type howmanycores: int
    :return: retrieved data
    :rtype: dict
    """

    if mode not in ['data', 'confdata', 'energy']:
        raise ValueError('unknown mode; choose between confdata, data, energy')

    args = zip(ids, itertools.repeat(mode), itertools.repeat(0))  # zipped args for multiprocessing

    with multiprocessing.Pool(howmanycores) as p:  # num_workers goes with number of cores of computer
        datas = p.starmap(access, args)

    data = dict(zip(ids, datas))

    return data


def access_vburmin_conf(id):
    # TODO: this is slow. Working on access_multi()

    data = access(id, mode='data', verbose=0)
    if data is None:
        return None
    else:
        vburmin_conf = data['vburminconf']
        del data
        confdata = access(id, mode='confdata', verbose=0)
        vburmin_confdata = confdata[vburmin_conf]

    return vburmin_confdata


def _access_with_persistent_http(k_ids, mode=None, verbose=1):
    """
    Use a persistent http connection to keep querying github
    Not really faster, at least when first accessing the data
    Weirdly there is some speed up if query the second time (sometimes)
    """

    github_url = 'https://raw.githubusercontent.com/doyle-lab-ucla/krkn/main/raw/'

    if mode == 'confdata':
        suffix = '_confdata.yml'
    elif mode == 'data':
        suffix = '_data.yml'
    elif mode == 'energy':
        suffix = '_relative_energies.csv'
    else:
        raise ValueError('unknown mode; choose between confdata, data, energy')

    # add leading zeros, add suffix to complete url
    full_ids = [str(k_id).zfill(8) for k_id in k_ids]
    headers = [full_id + suffix for full_id in full_ids]

    s = requests.Session()

    data = {}
    for ii in range(len(k_ids)):
        r = s.get(github_url+headers[ii])
        if r.status_code == 404:
            if verbose:
                print('ligand {0} not found'.format(full_ids[ii]))
            data[full_ids[ii]] = None
        elif r.status_code == 200:
            if mode == 'energy':
                data[full_ids[ii]] = pd.read_csv(StringIO(r.text), delimiter=';')
            else:
                data[full_ids[ii]] = yaml.load(r.text, Loader=yaml.Loader)
        else:
            raise ConnectionError('unknown networking issue, status code {0}'.format(r.status_code))

    return data


def _access_speed_test():
    # 50 ligands
    # Access_multi_v1: 60.582969332999994s
    # For loop with access(): 265.57729820400004s

    with open('buchwald/buchwald_found.json', 'r') as f:
        bids = json.load(f)

    # # for loop access() vs. access_multi()
    #
    # start = timeit.default_timer()
    # access_multi(test_list, mode='confdata')
    # print('Access_multi_v1: {0}s'.format(timeit.default_timer() - start))
    #
    # start = timeit.default_timer()
    # for l in test_list:
    #     access(l, mode='confdata')
    # print('For loop with access(): {0}s'.format(timeit.default_timer() - start))

    return None


# for testing only
if __name__ == '__main__':

    with open('buchwald/buchwald_found.json', 'r') as f:
        bids = json.load(f)

    output = access_multi(bids, mode='data')
    print(output[1])

    #access_multi_v1(bids, mode='energy')
    
    # with open('buchwald/buchwald_found.json', 'r') as f:
    #     bids = json.load(f)
    
    # values = ['cyjohnphos', 'brett', 'pcy3', 'pph3']
    # identifier = 'name'
    # print(
    #     lookup_multi(values=values, identifier=identifier)
    # )

    # print(lookup(keywords=['Ph', 'tBu', 'OMe', 'jason']))