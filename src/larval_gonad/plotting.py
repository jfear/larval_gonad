from functools import wraps

from numpy import arange
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from .config import config

# colormaps
cluster_cmap = dict(zip(config["cluster_order"], config["colors"]["clusters"]))
chrom_cmap = dict(zip(config["chrom_order"], config["colors"]["chrom"]))

# I have a separate color scheme for boxplots, that does not contain Y.
chrom_boxplot_cmap = dict(zip(config["chrom_order"][:-1], config["colors"]["chrom_boxplot"]))


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
        ax = kwargs.get("ax", None)
        if ax is None:
            fig, ax = plt.subplots()
            kwargs.update({"ax": ax})
        return func(*args, **kwargs)

    return wrapper


def make_figs(fname=None, styles=None, formats=None, layout=True, kws_layout=None):
    if isinstance(formats, str):
        formats = [formats]
    elif formats is None:
        formats = ["png", "eps"]

    if isinstance(styles, str):
        styles = [styles]
    elif styles is None:
        styles = ["notebook"]

    if kws_layout is None:
        kws_layout = {}

    def _plot_all(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            def plot_style(style, formats):
                with plt.style.context(["common", style]):
                    func(*args, **kwargs)
                    if layout:
                        plt.tight_layout(**kws_layout)
                    if ("notebook" not in style) & (fname is not None):
                        fn = fname + "_" + style
                        for f in formats:
                            plt.savefig("{}.{}".format(fn, f))
                        plt.close()

            for style in styles:
                plot_style(style, formats)

        return wrapper

    return _plot_all


def add_color_labels(ax, s=5, germ=False):
    clusters = config["sel_cluster_order"]

    if germ:
        clusters = clusters[:5]

    lclus = len(clusters)

    ax.set_xticks(arange(0, lclus + 1, 0.5))
    ax.set_xlim(0, lclus)

    for i, clus in enumerate(clusters):
        ax.plot(i + 0.5, 1, "bo", markersize=s, color=cluster_cmap[clus])
        sns.despine(ax=ax, left=True, bottom=True)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)


def add_color_labels_w_rep(ax, s=5, germ=False):
    clusters = config["sel_cluster_order_w_rep"]

    if germ:
        clusters = clusters[:12]

    lclus = len(clusters)

    ax.set_xticks(arange(0, lclus + 1, 0.5))
    ax.set_xlim(0, lclus)

    for i, clus in enumerate(clusters):
        ax.plot(i + 0.5, 1, "bo", markersize=s, color=cluster_cmap_w_rep[clus])
        sns.despine(ax=ax, left=True, bottom=True)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)


def flip_ticks(ax, pos="left"):
    ax.yaxis.set_ticks_position(pos)
    ax.yaxis.set_label_position(pos)


def add_triangle(ax, add_text=True, **kwargs):
    points = [[0, 0], [1, 0], [1, 1]]
    polygon = plt.Polygon(points, alpha=0.6)
    ax.add_artist(polygon)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    if add_text:
        ax.text(0.5, 0.1, "Pseudotime", ha="center")


def centerify(text, width=-1):
    """Center multiline text."""
    lines = text.split(" ")
    width = max(map(len, lines)) if width == -1 else width
    return "\n".join(line.center(width) for line in lines)


def dechr(ax, axis=0):
    """Remove chr prefix from axis labels."""
    labels = []
    if axis == 0:
        for label in ax.get_xticklabels():
            labels.append(label.get_text().replace("chr", ""))
        ax.set_xticklabels(labels)
    elif axis == 1:
        for lable in ax.get_yticklabels():
            labels.append(label.get_text().replace("chr", ""))
        ax.set_yticklabels(labels)


def format_pval(ax, x, y, pvalue, **kwargs):
    if pvalue <= 0.001:
        annotation = "***"
    elif pvalue <= 0.01:
        annotation = "**"
    elif pvalue <= 0.05:
        annotation = "*"
    else:
        return ax

    ax.text(x, y, annotation, ha="center", **kwargs)
    return ax


def add_pval(x, y, pval, ax, **kwargs):
    """Loops over and adds formatted p-value to plot.
    
    Parameters
    ----------
    x : array-like
        A list of x locations.
    y : array-like
        A list of y locations.
    pval : array-like
        A list of p-values or q-values for plotting.
    ax : plt.Axes
        The axes to add formatted p-values.
    """
    for x_i, y_i, pval_i in zip(x, y, pval):
        format_pval(ax, x_i, y_i, pval_i, **kwargs)


def plot_statsmodels_results(file: str, results: str):
    plt.text(0.01, 0.05, results, {'fontsize': 10}, fontproperties = 'monospace') # approach improved by OP -> monospace!
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(file, dpi=200)


def demasculinization(data, ax=None, title=None, legend=False, **kwargs):
    """Stacked barplot common for display of demasculinization.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame organized with chomosome as index, columns as Up/NS/Down,
        and values as proportions.
    ax : plt.Axes, optional
        Alternative axes to draw the plot, by default None
    title : str, optional
        Title to add to the plot, by default None
    legend : bool, optional
        If to keep the legend, by default False

    Returns
    -------
    plt.Axes
        Matplotlib axes with the stacked bar plot.

    Example
    -------
    >>> data = pd.DataFrame({"Up": [.1, .3], "NS": [.7, .4], "Down": [.2, .3]}, index=["X", "A"])
    >>> demasculinization(data)
    """
    if ax is None:
        fig, ax = plt.subplots()

    plot_defaults = dict(
        stacked=True,
        color=["red", "lightgrey", "blue"],
        width=0.9,
        edgecolor="k",
        linewidth=0.2
    )
    plot_defaults.update(kwargs)

    data.plot.bar(ax=ax, **plot_defaults)

    # Clean up X
    plt.setp(ax.get_xticklabels(), rotation=0)
    ax.set_xlabel("")

    # Clean up Y
    ax.set_ylim(0, 1)
    ax.set_ylabel("")

    # Clean up other parts
    if title is not None:
        ax.set_title(title)

    if not legend:
        ax.legend_ = None

    sns.despine(ax=ax, left=True, bottom=True)

    return ax
