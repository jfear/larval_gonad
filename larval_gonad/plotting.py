from functools import wraps

from numpy import arange
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from .config import config

# colormaps
cluster_cmap = dict(zip(config['cluster_order'], config['colors']['clusters']))
chrom_cmap = dict(zip(config['chrom_order'], config['colors']['chrom']))

# I have a separate color scheme for boxplots, that does not contain Y.
chrom_boxplot_cmap = dict(zip(config['chrom_order'][:-1], config['colors']['chrom_boxplot']))


def add_styles(dirname):
    mpl.style.core.USER_LIBRARY_PATHS.append(dirname)
    mpl.style.core.update_user_library(mpl.style.library)
    mpl.style.reload_library()


def make_ax(*args, **kwargs):
    fig, ax = plt.subplots(*args, **kwargs)
    return ax


def figure_element(func):
    """Decorate figure elements to create an Axes if none is given.

    I define a figure element as a function that imports, munges, and plots. A
    figure panel can be built by combining multiple figure elements using
    matplotlib.Axes. However, some times I will only want to plot the element,
    so this wrapper will generate a new Axes and pass it to the function for
    me.

    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        ax = kwargs.get('ax', None)
        if ax is None:
            fig, ax = plt.subplots()
            kwargs.update({'ax': ax})
        return func(*args, **kwargs)
    return wrapper


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


def add_color_labels(ax, s=5, germ=False):
    clusters = config['cluster_order']

    if germ:
        clusters = clusters[:5]

    lclus = len(clusters)

    ax.set_xticks(arange(0, lclus + 1, 0.5))
    ax.set_xlim(0, lclus)

    for i, clus in enumerate(clusters):
        ax.plot(i + 0.5, 1, 'bo', markersize=s, color=cluster_cmap[clus])
        sns.despine(ax=ax, left=True, bottom=True)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)


def flip_ticks(ax, pos='left'):
    ax.yaxis.set_ticks_position(pos)
    ax.yaxis.set_label_position(pos)


def add_triangle(ax, add_text=True, **kwargs):
    points = [[0, 0], [1, 0], [1, 1]]
    polygon = plt.Polygon(points, alpha=.6)
    ax.add_artist(polygon)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    if add_text:
        ax.text(0.5, .1, 'Pseudotime', ha='center', )
