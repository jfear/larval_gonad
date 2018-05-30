import itertools
from functools import wraps

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


def add_styles(dirname):
    mpl.style.core.USER_LIBRARY_PATHS.append(dirname)
    mpl.style.core.update_user_library(mpl.style.library)


def make_ax(*args, **kwargs):
    fig, ax = plt.subplots(*args, **kwargs)
    return ax


def make_figs(fname=None, styles=None, formats=None, layout=True,
              kws_layout=None):
    if isinstance(formats, str):
        formats = [formats]
    elif formats is None:
        formats = ['png', 'eps']

    if isinstance(styles, str):
        styles = [styles]
    elif styles is None:
        styles = ['notebook']

    if kws_layout is None:
        kws_layout = {}

    def _plot_all(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            def plot_style(style, formats):
                with plt.style.context(['common', style]):
                    func(*args, **kwargs)
                    if layout:
                        plt.tight_layout(**kws_layout)
                    if ('notebook' not in style) & (fname is not None):
                        fn = fname + '_' + style
                        for f in formats:
                            plt.savefig('{}.{}'.format(fn, f))
                        plt.close()

            for style in styles:
                plot_style(style, formats)

        return wrapper
    return _plot_all


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
