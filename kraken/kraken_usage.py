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


def merck_hte():

    """
    Written for notebook, some parts might look weird for a function

    :return:
    """

    import json
    from sklearn import linear_model
    from sklearn import preprocessing
    from kraken_utils import lookup_multi, featurize

    # HTE data from Merck
    data = pd.read_csv('https://raw.githubusercontent.com/beef-broccoli/ochem-data/main/Vbur/rxns/smiles.csv')

    # Select one set of experiments
    chloride_smiles = 'ClC1=CC=C(OC)C=C1'
    boronic_acid_smiles = 'OB(O)C1=CC=C(C(F)(F)F)C=C1'
    data = data.loc[(data['chloride_smiles'] == chloride_smiles) & (data['boronic_acid_smiles'] == boronic_acid_smiles)]
    data = data[['ligand_smiles', 'ligand_name', 'yield']]

    # lookup kraken ids for these ligands by name (smiles also work, not as well)
    names = data['ligand_name'].to_list()
    ids = lookup_multi(names, 'name')
    # smi = data['ligand_smiles'].to_list()
    # ids = lookup_multi(smi, 'smiles')
    data['id'] = ids

    # featurize
    features, _ = featurize(ids)
    Xs = features.drop(['id'], axis=1)
    feature_names = list(Xs.columns)

    # plot the reactivity threshold
    with open('id_list/buchwald_found.json') as f:
        b = json.load(f)
    non_buchwald_ids = []
    buchwald_ids = []
    for i in ids:
        if i in b:
            buchwald_ids.append(i)
        else:
            non_buchwald_ids.append(i)

    v_nb = [v/180*100 for v in list(features.loc[features['id'].isin(non_buchwald_ids)]['vbur_vbur_min'])]
    y_nb = list(data.loc[data['id'].isin(non_buchwald_ids)]['yield'])
    v_b = [v/180*100 for v in list(features.loc[features['id'].isin(buchwald_ids)]['vbur_vbur_min'])]
    y_b = list(data.loc[data['id'].isin(buchwald_ids)]['yield'])
    plt.scatter(v_nb, y_nb, zorder=3, label='non-buchwald ligand')
    plt.scatter(v_b, y_b, zorder=4, label='buchwald ligand')
    plt.axhline(10, c='k', linestyle='dashed', zorder=1)
    plt.axvline(33, c='k', linestyle='dashed', zorder=2)
    plt.xlabel('%Vbur(min)')
    plt.ylabel('yield (%)')
    plt.legend()

    # Logistic regression with L1 penalty to select features
    X = Xs.to_numpy()
    X = preprocessing.StandardScaler().fit_transform(X)
    y = np.array([i > 10 for i in list(data['yield'])]).astype(int)  # 10% yield as cutoff
    clf = linear_model.LogisticRegressionCV(cv=5, solver='liblinear', penalty='l1').fit(X, y)
    coeff = abs(clf.coef_[0, :])
    n = sum(clf.coef_[0, :] != 0)  # how many non-zero features
    selected_features = [x for _, x in sorted(zip(coeff, feature_names), reverse=True)]  # sorted by abs(coeff) for all coeff
    selected_features = selected_features[:n]

    # do LASSO regression again with the selected features
    X = Xs[selected_features].to_numpy()
    X = preprocessing.StandardScaler().fit_transform(X)
    y = data['yield']
    reg = linear_model.LassoCV(cv=5).fit(X, y)

    # screen all kraken ligands
    all = pd.read_csv('kraken_features_only.csv')
    ids = all['id']
    all = all.drop(['id'], axis=1)
    all = all[selected_features]
    all_X = preprocessing.StandardScaler().fit_transform(all.to_numpy())  # TODO: this is wrong. Should scale with training samples together
    all_y_pred = reg.predict(all_X)
    sorted_id_list = sorted(zip(all_y_pred, ids), reverse=True)

    # select top-10 ligands
    suggest = []
    for i in range(10):
        suggest.append(sorted_id_list[i][1])
    identifiers = pd.read_csv('identifiers.csv')
    suggest_names = identifiers.loc[identifiers['id'].isin(suggest)]['ligand']

    # plot non-zero features
    n = sum(np.array(reg.coef_) != 0)  # how many non-zero features
    coeff = abs(reg.coef_)
    ind = (-coeff).argsort()  # argsort() descending order
    selected2_features = np.array(selected_features)[ind][:n]
    weights = reg.coef_[ind][:n]

    ys = np.arange(len(weights))
    fig, ax = plt.subplots()
    ax.barh(ys, weights)
    ax.set_yticks(ys, labels=selected2_features)
    ax.invert_yaxis()
    ax.set_xlabel('coefficient')
    ax.set_xticks(np.arange(int(min(weights))-1, int(max(weights))+1, 1))



if __name__ == '__main__':

    merck_hte()

