import json
import yaml
import requests

from kraken_utils import access_multi, access_one_conf


"""
Demo of functionality
For a list of ligands, find the vbur_min conformer, and retrieve the vbur 
"""
def get_vbur_for_vburmin_conf(ids):

    data, ids_no_data = access_multi(ids, mode='data')
    ids_yes_data = list(data.keys())

    vburmin_conf_names = {}
    for id in ids_yes_data:
        vburmin_conf_names[id] = data[id]['vburminconf']

    ids_and_confs = [str(id) + ':' + vburmin_conf_names[id] for id in ids_yes_data]
    data, _ = access_multi(ids_and_confs, access_func=_get_vbur, mode='confdata')

    print(data)

    # data, _ = access_multi(ids_yes_data, mode='confdata')
    # temp = data[1][vburmin_conf_names[1]]


    return None


def _get_vbur(id_and_conf, mode, verbose=1):

    github_url = 'https://raw.githubusercontent.com/doyle-lab-ucla/krkn/main/raw/'

    suffix = '_confdata.yml'

    id, conf = id_and_conf.split(':')

    full_id = str(id).zfill(8)

    r = requests.get(github_url + full_id + suffix)
    if r.status_code == 200:
        data = yaml.load(r.text, Loader=yaml.Loader)
        return data[conf]['properties']['vbur_vbur']
    elif r.status_code == 404:
        if verbose:
            print('ligand {0} not found'.format(full_id))
        return None
    else:
        raise ConnectionError('unknown networking issue, status code {0}'.format(r.status_code))


if __name__ == '__main__':

    with open('buchwald/buchwald.json', 'r') as f:
        bids = json.load(f)

    bids = bids[0:10]

    get_vbur_for_vburmin_conf(bids)
