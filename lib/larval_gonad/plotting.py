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


def make_figs(fname=None, styles=None, formats=None, kws_layout=None):
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
                with plt.style.context([style, 'common']):
                    func(*args, **kwargs)
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


def TSNEPlot(x, y, data=None, hue=None, cmap=None, palette=None, ax=None,
             class_names=None, legend_kws=None, **kwargs):
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
        'vmax': 5
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
    if isinstance(values[0], (int, float)) & (len(values) > 2):
        if cmap is None:
            cmap = mpl.colors.ListedColormap(palette)
        df.plot.scatter(x, y, c=df[hue], cmap=cmap, ax=ax, **defaults)
    else:
        if cmap is None:
            cmap = {k: v for k, v in zip(values, palette)}

        for l, dd in df.groupby(hue):
            if class_names is None:
                _class = str(l).title()
            else:
                _class = class_names[l]

            dd.plot.scatter(x, y, c=cmap[l], label=_class, ax=ax, **defaults)

        ax.legend(**legend_defaults)


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
