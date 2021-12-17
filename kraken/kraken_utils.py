import requests
from io import StringIO
import multiprocessing
import itertools
import yaml
import json
import pandas as pd
import numpy as np
from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors


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

    identifiers = pd.read_csv('https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/kraken/identifiers.csv')
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
            results = [n for n in name_list if (name.lower() in n.lower()) or (n.lower() in name.lower())]

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

        # catch invalid smiles
        # Note: rdkit will throw an error for invalid smiles. Since warnings are done at C++ level, there is no easy
        # way to catch it in python. The function still runs.
        # Try:
        #   import rdkit.RDLogger
        #   RDLogger.DisableLog('rdapp.*')

        m = Chem.MolFromSmiles(smiles)
        if m is None:
            temp_df = pd.DataFrame()
        else:
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


def lookup_multi(values, identifier, dataspell=False):
    """
    Builds on top of lookup()
    Input: a list of identifier values, type of identifiers
    Returns: if found one result, return kraken ID. If multiple or no results are found, outputs possibilities and prompt
    user to enter the correct kraken ID for ligand

    :param values:
    :param identifier:
    :param dataspell: if this notebook is run in dataspell (better pandas table rendering, no need to display all)
    :return:
    """
    if dataspell is False:
        pd.set_option("display.max_rows", None, "display.max_columns", None, 'display.max_colwidth', None)

    if identifier not in ['cas', 'name', 'smiles', 'inchi']:
        raise ValueError('identifier needs to be cas, name, smiles or inchi. (keyword is not supported)')

    results = []
    for val in values:
        d = {identifier: val,
             'verbose': 0}
        lookup_result = lookup(**d)
        if lookup_result is None:  # no results, return original value
            results.append(val)
        elif len(lookup_result) == 1:
            results.append(lookup_result['id'].values[0])
        else:  # no single match, return original identifier value
            print('For ligand with \"{0}\" \"{1}\", no match is found. Here are some possibilities:'.format(identifier, val))
            print(lookup_result)
            input_id = int(input('Please input a kraken ID for this ligand (0 if none matches): '))
            if input_id:
                results.append(input_id)
            else:
                results.append(val)

    pd.reset_option("display.max_rows")
    pd.reset_option("display.max_columns")
    pd.reset_option('display.max_colwidth')
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
    :rtype: Union[pandas.core.frame.DataFrame, dict]
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
    if r.status_code == 200:
        if mode == 'energy':
            return pd.read_csv(StringIO(r.text), delimiter=';')
        else:
            return yaml.load(r.text, Loader=yaml.Loader)
    elif r.status_code == 404:
        if verbose:
            print('ligand {0} not found'.format(full_id))
        return None
    else:
        raise ConnectionError('unknown networking issue, status code {0}'.format(r.status_code))


def access_one_conf(id_and_conf, mode='confdata', verbose=1):
    """
    Access data for one conformer
    Hacks it, combines id and conf into one argument.
    Note: Need to parse the string (use : as delimiter when combining id and conf_name)

    :param id_and_conf:
    :param mode:
    :param verbose:
    :return:
    """

    github_url = 'https://raw.githubusercontent.com/doyle-lab-ucla/krkn/main/raw/'

    mode = 'confdata'
    suffix = '_confdata.yml'

    id, conf = id_and_conf.split(':')

    full_id = str(id).zfill(8)

    r = requests.get(github_url + full_id + suffix)
    if r.status_code == 200:
        data = yaml.load(r.text, Loader=yaml.Loader)
        return data[conf]
    elif r.status_code == 404:
        if verbose:
            print('ligand {0} not found'.format(full_id))
        return None
    else:
        raise ConnectionError('unknown networking issue, status code {0}'.format(r.status_code))


def access_multi(ids, access_func=access, mode=None, howmanycores=8):
    """
    Access the data for a list ligands. Use multiprocessing to speed up access()
    Note: if use mode='confdata', it will take a while to fetch (~1 min for 50 ligands)

    :param ids: list of kraken ids
    :type ids: list of int
    :param access_func: access function to fetch data for individual ligand (default: access() fetches all data
    specified only by mode)
    :type access_func: function
    :param mode: 'confdata', 'data', 'energy'
    :type mode: str
    :param howmanycores: number of cores on this computer (decides how many parallel processes can be initiated)
    :type howmanycores: int
    :return:
        - data_filtered - a dict of retrieved data
        - ids_with_no_data - a list of ids(int) with no data available
    """

    if mode not in ['confdata', 'data', 'energy']:
        raise ValueError('unknown mode; choose between confdata, data, energy')

    args = zip(ids, itertools.repeat(mode), itertools.repeat(0))  # zipped args for multiprocessing

    with multiprocessing.Pool(howmanycores) as p:  # num_workers goes with number of cores of computer
        datas = p.starmap(access_func, args)

    data = dict(zip(ids, datas))  # outputs a dict with None
    ids_no_data = [k for k, v in data.items() if v is None]
    data_filtered = {k: v for k, v in data.items() if v is not None}  # seems wasteful, but probably okay for now...

    # maybe get ids_yes_data from data_filtered.keys(), then do a set difference with ids to get ids_no_data

    return data_filtered, ids_no_data


def featurize_ligand(id):
    """
    for each ligand, create a dataframe with conformer name as rows, properties as column

    :param id: kraken id for ligand
    :type id: int
    :return: dataframe with conformer name as rows and properties as columns
    :rtype: pandas.core.frame.DataFrame
    """

    data = access(id, mode='confdata')

    energy_param_list = ['e_dz', 'e_tz_gas', 'e_tz_solv', 'g', 'g_tz_gas', 'g_tz_solv']

    # flatten dictionaries. Output dict: first keys conf_name, second keys feature name
    filtered_data = {}
    for conf, confdata in data.items():
        energies = {l: confdata[l] for l in energy_param_list}
        properties = confdata['properties']
        filtered_data[conf] = {**energies, **properties}

    df = pd.DataFrame.from_dict(filtered_data, orient='index')

    return df


def featurize(ids):

    """
    For a list of ligands, create a dataframe with ligand id as rows, properties as column

    :param ids:
    :return:
    """

    # TODO: select features: sterics, electronics, interactions...
    # TODO: scale, corr analysis

    features = pd.read_csv('https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/kraken/kraken_features_only.csv')

    feature_ids = set(features['id'])
    query_ids = set(ids)
    no_data = query_ids.difference(feature_ids)

    features = features.loc[features['id'].isin(ids)]

    # add a name for features, since resulting features are sorted by ids
    identifiers = pd.read_csv('https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/kraken/identifiers.csv')
    features.insert(0, 'name', features['id'].apply(lambda x: identifiers.loc[identifiers['id'] == x]['ligand'].values[0]))

    return features, no_data


def fetch_xyz(id, conf_name, file_path=None, metal='Ni'):
    """
    for a conformer, write .xyz file

    :param id: kraken id
    :type id: int
    :param conf_name: name of the conformer
    :type conf_name: str
    :param file_path: file path to save .xyz file
    :type file_path: str
    :param metal: metal center used in calculation
    :type metal: str
    :return: string block that goes into the xyz file (in case you don't want to write to file)
    :rtype: str
    """

    data = access_one_conf(str(id) + ':' + conf_name)
    if metal == 'Ni':
        coords = data['coords']
        elements = data['elements']
    elif metal == 'Pd':
        coords = data['coords_pd']
        elements = data['elements_pd']
    else:
        raise ValueError('Ni or Pd for metal center')

    xyz = str(len(coords)) + '\n'
    xyz += str(id) + ' ' + conf_name + '\n'
    for ii in range(len(elements)):
        str_ints = [str(t) for t in coords[ii]]
        s = ' '.join(str_ints)
        xyz += elements[ii] + ' ' + s + '\n'

    if file_path is not None:
        with open(file_path, 'w') as f:
            f.write(xyz)

    return xyz


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
    # Access_multi: 60.582969332999994s
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
    # identifiers = pd.read_csv('https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/kraken/identifiers.csv')
    # print(identifiers.loc[identifiers['id'] == 3]['ligand'].values[0])

    print(featurize([1,2,3]))

    #fetch_xyz(1, '00000001_Ni_00014', './scratch/test.xyz', metal='Pd')

    # data = access(4, mode='confdata')
    # print(data.keys())

    # with open('buchwald/buchwald.json', 'r') as f:
    #     bids = json.load(f)
    #
    # data, ids = access_multi(bids, mode='energy')
    # print(data[1])

    # output, l = access_multi([1, 2, 359, 360], mode='data')
    # print(l)

    #access_multi_v1(bids, mode='energy')

    # with open('buchwald/buchwald_found.json', 'r') as f:
    #     bids = json.load(f)

    # values = ['cyjohnphos', 'brett', 'pcy3', 'pph3']
    # identifier = 'name'
    # print(
    #     lookup_multi(values=values, identifier=identifier)
    # )

    # print(lookup(keywords=['Ph', 'tBu', 'OMe', 'jason']))