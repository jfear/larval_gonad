"""Common elements used for the scRNA-Seq data.

This is a place to store functions and variables that are used over and over
again with the scRNA-seq data.

"""
import itertools

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


CLUSTER_ANNOT = {
    0: 'Late 1º Spermatocytes (0)',
    1: 'Mid Cyst Cells (1)',
    2: 'Mid 1º Spermatocytes (2)',
    3: 'Early 1º Spermatocytes (3)',
    4: 'Late Cyst Cells (4)',
    5: 'Early Cyst Cells (5)',
    6: 'Spermatogonia (6)',
    7: 'Terminal Epithelium (7)',
    8: 'Pigment Cells (8)',
    9: 'Unknown (9)',
    10: 'Unknown (10)',
    11: 'Unknown (11)',
}

CLUSTER_ORDER = [
    'Spermatogonia (6)',
    'Early 1º Spermatocytes (3)',
    'Mid 1º Spermatocytes (2)',
    'Late 1º Spermatocytes (0)',
    'Early Cyst Cells (5)',
    'Mid Cyst Cells (1)',
    'Late Cyst Cells (4)',
    'Terminal Epithelium (7)',
    'Pigment Cells (8)',
    'Unknown (9)',
    'Unknown (10)',
    'Unknown (11)',
]


# CLUSTER_ANNOT = {
#     0: 'Late Primary Spermatocytes (0)',
#     1: 'Early Somatic Cyst Cells (1)',
#     2: 'Late Somatic Cyst Cells (2)',
#     3: 'Late Somatic Cyst Cells (3)',
#     4: 'Spermatogonia (4)',
#     5: 'Terminal Epithelium (5)',
#     6: 'Mid Primary Spermatocytes (6)',
#     7: 'Late Somatic Cyst Cells (7)',
#     8: 'Early Primary Spermatocytes (8)',
#     9: 'Pigment Cells (9)',
#     10: 'Early Somatic Cyst Cells (10)',
#     11: 'Unknown (11)',
# }
#
# CLUSTER_ORDER = [
#     'Spermatogonia (4)',
#     'Early Primary Spermatocytes (8)',
#     'Mid Primary Spermatocytes (6)',
#     'Late Primary Spermatocytes (0)',
#     'Early Somatic Cyst Cells (1)',
#     'Early Somatic Cyst Cells (10)',
#     'Late Somatic Cyst Cells (2)',
#     'Late Somatic Cyst Cells (3)',
#     'Late Somatic Cyst Cells (7)',
#     'Terminal Epithelium (5)',
#     'Pigment Cells (9)',
#     'Unknown (11)',
# ]

# List of T2A Gal4 lines we have looked at.
T2A = [
    'vari',
    'CG42458',
    'CG34383',
    'CG34394',
    'CG7255',
    'Wnt4',
    'cv-2',
    'Eaf',
    'tun',
    'Zasp52',
    'Su(var)2-10',
    'GEFmeso',
    'ths',
    'pk',
    'rdo',
    'bin3',
    'CG2082',
    'Khc-73',
    'FER',
    'CG31075',
    'QC',
    'dally',
    'CG11658',
    'PyK',
    'PH4alphaEFB',
    'RasGAP1',
    'Ino80',
]

# List of GAL4 lines we have looked at
GAL4 = [
    'CG11658',
    'Notum',
    'CadN ',
    'Tsp74F',
    'sano',
    'rau',
    'CG32354',
    'trn',
    'Irk1',
    'bnl',
    'hng3',
    'tkv',
    'CG34383',
    'CG34323',
    'qjt',
    'qin',
    'qin',
    'klu',
    'CR46216',
    'mael',
    'CG9272',
    'CR33319',
    'Grx1t',
    'CG43059',
    'CG6888',
    'CG4218',
    'CR44852',
    'CG31644',
    'bol',
    'svp',
    'svp',
    'CG4822',
    'Meltrin',
    'Fili',
    'Papss',
    'AdamTS-A',
    'cdi',
    'eya',
    'CG42673',
]

# List of Ben White's T2A gal4 that we have looked at.
WHITE = [
    'AkhR',
    'burs',
    'CG2187',
    'ChaT',
    'CrzR',
    'Dh31R',
    'EH',
    'EHR',
    'ETHR',
    'ETHR',
    'ETHRB',
    'GluRIIB',
    'HR46',
    'HR46',
    'LGR1',
    'pburs',
    'Proc',
    'ProcR',
    'rk',
]

# List of protein traps that Erika has looked at.
PTRAP = [
    'ADD1', 'Ance', 'bbg', 'bol', 'CadN', 'CG11044',
    'CG17349', 'CG17646', 'CG3277', 'CG42321', 'CG43088', 'CG43373',
    'CG6356', 'CG8036', 'CG8851', 'CG9747', 'cindr', 'Dek',
    'DIP-delta', 'dlp', 'dpr17', 'dpr3', 'dsx', 'Efa6',
    'egr', 'eIF5B', 'e(y)3', 'Fas3', 'foxo', 'fru',
    'Fs(2)Ket', 'hid', 'ImpL2', 'klu', 'LRP1', 'Mapmodulin',
    'mbl', 'Mi-2', 'NetA', 'NFAT', 'nkd', 'nord',
    'Nrg', 'osa', 'p130CAS', 'p53', 'Pdcd4', 'Pde11',
    'Piezo', 'ppan', 'Ptp4E', 'Pvr', 'Rab3', 'Reck',
    'Sap-r', 'sdk', 'Sox102F', 'spir', 'SRPK', 'stai',
    'stg', 'svp', 'Ten-a', 'Ten-m', 'Tep2', 'tok',
    'tutl', 'twin', 'Vha68-3', 'Wdr62', 'wnd',
]


def TSNEPlot(x='tSNE_1', y='tSNE_2', data=None, hue=None, cmap=None,
             palette=None, ax=None, class_names=None,
             legend_kws=None, cbar=True, **kwargs):
    """ Make a TSNE plot using either continuous or discrete data.

    Parameters
    ----------
    x : str
        tSNE to plot on x-axis
    y : str
        tSNE to plot on y-axis
    data : pd.DataFrame
        DataFrame containg tSNEs and any additional meta data.
    hue : str
        Column in the data frame to color by. If column is strings, colors will
        be assined using current color palette. If column are numbers, colors
        will be a heatmap. If colors are bool then a grey/black color
        scheme is used.
    cmap : dict or ListedColorMap
        A ditionary mapping of the colors to use with which labels. Or a
        matplotlib colormap.
    ax : matplotlib.axes
        Axes which to plot.
    kwargs :
        Additional kwargs are passed to pd.DataFrame.plot.scatter

    """
    defaults = {
        'vmin': 0,
        'vmax': 5,
        'edgecolor': 'k',
    }
    defaults.update(kwargs)

    df = data.copy()
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    if palette is None:
        palette = sns.color_palette()

    legend_defaults = {
        'loc': 'center left',
        'bbox_to_anchor': (1, 0.5),
    }
    if legend_kws is not None:
        legend_defaults.update(legend_kws)

    if isinstance(hue, pd.Series):
        df['on'] = hue.astype(int).apply(lambda x: str(x))
        hue = 'on'

    values = sorted(df[hue].unique())
    if (isinstance(values[0], (int, np.integer)) & (len(values) > 20)) | isinstance(values[0], (float, np.float64)):
        if cmap is None:
            cmap = mpl.colors.ListedColormap(palette)
        zeros = df[df[hue] == 0]
        if len(zeros) > 0:
            zeros.plot.scatter(x, y, c=zeros[hue], cmap=cmap, ax=ax,
                               colorbar=False, zorder=1, **defaults)

        expressed = df[df[hue] > 0]
        if len(expressed) > 0:
            expressed.plot.scatter(x, y, c=expressed[hue], cmap=cmap, ax=ax,
                                   colorbar=cbar, zorder=3, **defaults)
    else:
        if cmap is None:
            cmap = {k: v for k, v in zip(values, palette)}

        for l, dd in df.groupby(hue):
            if class_names is None:
                _class = str(l).title()
            else:
                _class = class_names[l]

            try:
                dd.plot.scatter(x, y, c=cmap[l], label=_class, ax=ax,
                                **defaults)
            except KeyError as e:
                print('Try Setting Palette with the correct number of colors.')
                print(len(cmap), len(values))
                print(type(values[0]))
                raise e

        ax.legend(**legend_defaults)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_yticks([])


def plot_confusion_matrix(cm, classes, normalize=False, title='Confusion matrix',
                     cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')


