import json

from kraken_utils import access_multi, access_one_conf


def access_vburmin_conf(ids):

    data, ids_no_data = access_multi(ids, mode='data')
    ids_yes_data = list(data.keys())

    vburmin_conf_names = {}
    for id in ids_yes_data:
        vburmin_conf_names[id] = data[id]['vburminconf']

    ids_and_confs = [str(id) + ':' + vburmin_conf_names[id] for id in ids_yes_data]

    data, ids_no_data = access_multi(ids_and_confs, access_func=access_one_conf, mode='confdata')

    print(data)

    # data, _ = access_multi(ids_yes_data, mode='confdata')
    # temp = data[1][vburmin_conf_names[1]]


    return None


if __name__ == '__main__':

    with open('buchwald/buchwald.json', 'r') as f:
        bids = json.load(f)

    bids = bids[0:2]

    access_vburmin_conf(bids)
