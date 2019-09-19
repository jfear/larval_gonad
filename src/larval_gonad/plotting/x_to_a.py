import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import shelve_load
from larval_gonad.plotting.common import cluster_boxplot
from larval_gonad.plotting.stats import add_pvals


def plot_x_to_a(shelve: str, cluster_color: list, cluster_order: list, ax: plt.Axes = None):
    """Boxplot of X:A ratios by cluster.

    Example
    -------
    >>> shelve = "output/x-to-a-wf/db/commonly_expressed.bak"
    >>> from larval_gonad.config import read_config
    >>> config = read_config("config/common.yaml")
    >>> color_config = read_config("config/colors.yaml")
    >>> cluster_color = color_config['clusters']
    >>> cluster_order = config['cluster_order']

    """
    df, pvals = _get_data(shelve, "x_to_a_ratio")

    if ax is None:
        _, ax = plt.subplots()

    cluster_boxplot("ratio", df, cluster_color, cluster_order, ax)
    ax.set(xlabel="", ylabel="X / A", axisbelow=True)
    ax.axhline(1, ls="--", color="grey")
    add_pvals(pvals.x, pvals.y, pvals.pvalue, ax)


def plot_4_to_a(shelve: str, cluster_color: list, cluster_order: list, ax: plt.Axes = None):
    """Boxplot of X:A ratios by cluster.

    Example
    -------
    >>> shelve = "output/x-to-a-wf/db/commonly_expressed.bak"
    >>> from larval_gonad.config import read_config
    >>> config = read_config("config/common.yaml")
    >>> color_config = read_config("config/colors.yaml")
    >>> cluster_color = color_config['clusters']
    >>> cluster_order = config['cluster_order']

    """
    df, pvals = _get_data(shelve, "fourth_to_a_ratio")

    if ax is None:
        _, ax = plt.subplots()

    cluster_boxplot("ratio", df, cluster_color, cluster_order, ax)
    ax.set(xlabel="", ylabel="4 / A", axisbelow=True)
    ax.axhline(1, ls="--", color="grey")
    add_pvals(pvals.x, pvals.y, pvals.pvalue, ax)


def _get_data(shelve, ratio_type):
    """Load data from shelves"""
    db = shelve_load(shelve)
    df = db["data"].pipe(lambda x: x[x.ratio_type == ratio_type])
    pvals = db["pvalues"].pipe(lambda x: x[x.ratio_type == ratio_type])
    return df, pvals
