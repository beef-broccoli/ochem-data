import json
import yaml
import requests
import pandas as pd
import numpy as np
import itertools
import matplotlib.pyplot as plt

from kraken_utils import access_multi, access_one_conf


"""
Demo of functionality
For a list of ligands, find the vbur_min conformer, and retrieve the vbur 
"""
def get_percentvbur_for_vburmin_conf(ids):

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


def plot_vbur_min_boltz():
    """
    plot %Vbur(min) + %Vbur(boltz) for buchwald ligands

    :param ids:
    :return:
    """

    with open('id_list/buchwald_found.json', 'r') as f:
        ids = json.load(f)
    with open('id_list/commercial.json', 'r') as f:
        comm_ids = json.load(f)

    ids = list(set(ids) & set(comm_ids))

    data, ids_no_data = access_multi(ids, mode='data')
    ids_yes_data = list(data.keys())

    mins = []
    boltzs = []
    for id in ids_yes_data:
        mins.append(data[id]['min_data']['vbur_vbur']/180*100)  # percent vbur, divide by the volume of a 3.5A-radius ball
        boltzs.append(data[id]['boltzmann_averaged_data']['vbur_vbur']/180*100)

    # fetch names
    identifiers = pd.read_csv('identifiers.csv')
    labels = []
    for id in ids_yes_data:
        name = identifiers.loc[identifiers['id'] == id]['ligand'].values[0]
        labels.append(str(id) + '.' + str(name))

    df = pd.DataFrame(data={'labels': labels,
                            '%Vbur(min)': mins,
                            '%Vbur(boltz)': boltzs,
                            'sum': [m+b for m,b in zip(mins, boltzs)]})

    df = df.loc[(df['sum']>80) & (df['sum']<85)]
    df = df.sort_values(by=['sum'])

    # draw a horizontal bar chart for discrete categorical values
    def categorical_bar(labels, data, category_names):
        """
            Parameters
            ----------
            labels : list
                A list of entries (individual reactions, ligands, chemicals...)
            data : numpy array
                Data as numpy array; for each label count each category
            category_names : list of str
                The category labels
        """
        # labels = list(results.keys())
        # data = np.array(list(results.values()))
        data_cum = data.cumsum(axis=1)
        category_colors = plt.get_cmap('Spectral')(
            np.linspace(0.15, 0.85, data.shape[1]))

        fig, ax = plt.subplots(figsize=(9.2, 5))
        ax.invert_yaxis()
        ax.xaxis.set_visible(True)
        ax.set_xlim(0, np.sum(data, axis=1).max())

        for i, (colname, color) in enumerate(zip(category_names, category_colors)):
            widths = data[:, i]
            starts = data_cum[:, i] - widths
            ax.barh(labels, widths, left=starts, height=0.5,
                    label=colname, color=color)
            xcenters = starts + widths / 2

            r, g, b, _ = color
            # text_color = 'white' if r * g * b < 0.5 else 'black'  # auto adjust text color based on color
            text_color = 'black'
            for y, (x, c) in enumerate(zip(xcenters, widths)):
                if int(c) != 0:
                    ax.text(x, y, str(round(c, 1)), ha='center', va='center', color=text_color)

        ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
                  loc='lower left', fontsize='medium')

        return fig, ax

    # plot
    fig, ax = categorical_bar(df['labels'], df[['%Vbur(min)', '%Vbur(boltz)']].to_numpy(), ['%Vbur(min)', '%Vbur(boltz)'])
    ax.axvline(80, ls='--', color='gray')
    ax.axvline(85, ls='--', color='gray')
    ax.set_xticks(np.arange(0,90,5))
    plt.show()
    return plt


def _get_vbur(id_and_conf, mode, verbose=1):

    github_url = 'https://raw.githubusercontent.com/doyle-lab-ucla/krkn/main/raw/'

    suffix = '_confdata.yml'

    id, conf = id_and_conf.split(':')

    full_id = str(id).zfill(8)

    r = requests.get(github_url + full_id + suffix)
    if r.status_code == 200:
        data = yaml.load(r.text, Loader=yaml.Loader)
        return data[conf]['properties']['vbur_vbur']/180
    elif r.status_code == 404:
        if verbose:
            print('ligand {0} not found'.format(full_id))
        return None
    else:
        raise ConnectionError('unknown networking issue, status code {0}'.format(r.status_code))


if __name__ == '__main__':

    plot_vbur_min_boltz()
