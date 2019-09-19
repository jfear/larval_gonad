import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from larval_gonad.plotting.common import cluster_boxplot


def plot_log_expression_by_cluster(
    feather: str, cluster_color: list, cluster_order: list, ax: plt.Axes = None, fbgns: list = None
):
    """Boxplot of total TPM expression
    
    Parameters
    ----------
    feather : str
        File name with TPM or UMP arranged in tidy format.
    cluster_color : list
    cluster_order : list
    ax : plt.Axes, optional
    fbgns : list, optional
        Subsets the data by these FBgns if provided (i.e., commonly expressed genes)

    Example
    -------
    >>> feather = "output/seurat3-cluster-wf/tpm_by_cluster.feather"
    >>> from larval_gonad.config import read_config
    >>> config = read_config("config/common.yaml")
    >>> color_config = read_config("config/colors.yaml")
    >>> cluster_color = color_config['clusters']
    >>> cluster_order = config['cluster_order']
    >>> from larval_gonad.io import pickle_load
    >>> fbgns = pickle_load("output/cellselection-wf/commonly_expressed_genes.pkl")

    """

    df = pd.read_feather(feather).set_index("FBgn")
    target_col = [x for x in df.columns if x != "cluster"][0]

    # Use log values for plotting
    log_values = f"Log({target_col})"
    df[log_values] = np.log1p(df[target_col])

    if fbgns:
        df = df[df.index.isin(fbgns)]

    if ax is None:
        _, ax = plt.subplots()

    cluster_boxplot(log_values, df, cluster_color, cluster_order, ax)

