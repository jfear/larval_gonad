import os
from functools import wraps

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
                with plt.style.context([style, 'default']):
                    func(*args, **kwargs)
                    plt.tight_layout(**kws_layout)
                    if (not 'notebook' in style) & (fname is not None):
                        fn = fname + '_' + style
                        for f in formats:
                            plt.savefig('{}.{}'.format(fn, f))
                        plt.close()

            for style in styles:
                plot_style(style, formats)

        return wrapper
    return _plot_all


def TSNEPlot(x, y, data=None, hue=None, cmap=None, palette=None, ax=None, **kwargs):
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

    if isinstance(hue, pd.Series):
        df['on'] = hue.astype(int).apply(lambda x: str(x))
        hue = 'on'

    values = sorted(df[hue].unique())
    if isinstance(values[0], (int, float)) & (len(values) > 2):
        if cmap is None:
            cmap = mpl.colors.ListedColormap(sns.color_palette('Reds'))
        df.plot.scatter(x, y, c=df[hue], cmap=cmap, ax=ax, **defaults)
    else:
        if cmap is None:
            cmap = {k: v for k, v in zip(values, palette)}

        for l, dd in df.groupby(hue):
            dd.plot.scatter(x, y, c=cmap[l], label=str(l).title(), ax=ax, **defaults)

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
