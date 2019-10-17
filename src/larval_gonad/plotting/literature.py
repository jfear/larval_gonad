from typing import List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from larval_gonad.io import feather_to_cluster_matrix
from larval_gonad.plotting.common import get_fbgn2symbol
from larval_gonad.plotting.biomarkers import _cleanup_xaxis as _cleanup_xaxis_rep


def plot_lit_evidence_profile(
    gene_metadata: str,
    lit_evidence: str,
    tpm_by_cluster: str,
    germ_clusters: list,
    axes: list = None,
):
    """Plot heatmap of evidence and expression patterns from the literature.

    Most of the evidence patterns are protein based.

    Example
    -------
    >>> from larval_gonad.config import read_config
    >>> config = read_config("config/common.yaml")
    >>> gene_metadata = f"references/gene_annotation_dmel_{config['tag']}.feather"
    >>> lit_evidence = "data/external/miriam/lit_gene_dummy_vars.tsv"
    >>> tpm_by_cluster = "output/seurat3-cluster-wf/tpm_by_cluster.feather"
    >>> germ_clusters = config["germ"]
    >>> plot_lit_expression_profile(gene_metadata, lit_evidence)

    """
    fbgn2symbol = get_fbgn2symbol(gene_metadata)
    germ_evidence = _get_lit_evidence(lit_evidence)[["SP", "ES", "MS", "LS"]]
    this_study = list(map(lambda x: fbgn2symbol[x], _genes_w_ptraps(lit_evidence)))
    binned_expression = _get_binned_expression(tpm_by_cluster)[germ_clusters]
    df = germ_evidence.join(binned_expression).rename(fbgn2symbol)

    if axes is None:
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(4, 8))
    else:
        ax1, ax2 = axes

    defaults = dict(square=True, linewidths=0.01, linecolor="k", yticklabels=True, cbar=False)
    sns.heatmap(data=df.iloc[:, :4], cmap=["#d3d3d3", "#450457", "#f8e621"], ax=ax1, **defaults)
    sns.heatmap(data=df.iloc[:, 4:], cmap=["#450457", "#ff7800", "#f8e621"], ax=ax2, **defaults)
    _cleanup_xaxis(ax1), _cleanup_yaxis(ax1, this_study)
    _cleanup_xaxis(ax2), _cleanup_yaxis(ax2, this_study)
    _add_legend(ax2)

def plot_lit_evidence_zscore_profile(
    gene_metadata: str,
    lit_evidence: str,
    zscore_by_cluster_rep: str,
    germ_clusters: list,
    axes: list = None,
):
    """Plot heatmap of evidence and expression patterns from the literature.

    Most of the evidence patterns are protein based.

    Example
    -------
    >>> from larval_gonad.config import read_config
    >>> config = read_config("config/common.yaml")
    >>> gene_metadata = f"references/gene_annotation_dmel_{config['tag']}.feather"
    >>> lit_evidence = "data/external/miriam/lit_gene_dummy_vars.tsv"
    >>> zscore_by_cluster_rep = "output/seurat3-cluster-wf/zscore_by_cluster_rep.feather"
    >>> germ_clusters = config["germ"]
    >>> plot_lit_expression_profile(gene_metadata, lit_evidence)

    """
    fbgn2symbol = get_fbgn2symbol(gene_metadata)
    germ_evidence = _get_lit_evidence(lit_evidence)[["SP", "ES", "MS", "LS"]]
    this_study = list(map(lambda x: fbgn2symbol[x], _genes_w_ptraps(lit_evidence)))
    zscore_expression = _get_zscore_expression(zscore_by_cluster_rep).loc[:, (germ_clusters, slice(None))]
    df = germ_evidence.join(zscore_expression).rename(fbgn2symbol)

    if axes is None:
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(4, 8))
    else:
        ax1, ax2 = axes

    defaults = dict(square=True, linewidths=0.01, linecolor="k", yticklabels=True, xticklabels=True, cbar=False)
    sns.heatmap(data=df.iloc[:, 4:], cmap='viridis', vmin=-3, vmax=3, ax=ax1, **defaults)
    sns.heatmap(data=df.iloc[:, :4], cmap=["#d3d3d3", "#450457", "#f8e621"], ax=ax2, **defaults)
    _cleanup_xaxis_rep(ax1, germ_clusters), _cleanup_yaxis(ax1, this_study)
    _cleanup_xaxis(ax2), _cleanup_yaxis(ax2, this_study)
    _add_legend(ax2)


def plot_lit_evidence_soma_profile(
    gene_metadata: str,
    lit_evidence: str,
    tpm_by_cluster: str,
    soma_clusters: list,
    axes: list = None,
):
    """Plot heatmap of evidence and expression patterns from the literature.

    Most of the evidence patterns are protein based.

    Example
    -------
    >>> from larval_gonad.config import read_config
    >>> config = read_config("config/common.yaml")
    >>> gene_metadata = f"references/gene_annotation_dmel_{config['tag']}.feather"
    >>> lit_evidence = "data/external/miriam/lit_gene_dummy_vars.tsv"
    >>> tpm_by_cluster = "output/seurat3-cluster-wf/tpm_by_cluster.feather"
    >>> soma_clusters = config["soma"]
    >>> plot_lit_expression_soma_profile(gene_metadata, lit_evidence)

    """
    fbgn2symbol = get_fbgn2symbol(gene_metadata)
    soma_evidence = _get_lit_evidence(lit_evidence)[["C", "EC", "MC", "LC", "PC", "TE"]]
    this_study = list(map(lambda x: fbgn2symbol[x], _genes_w_ptraps(lit_evidence)))
    binned_expression = _get_binned_expression(tpm_by_cluster)[soma_clusters]
    df = soma_evidence.join(binned_expression).rename(fbgn2symbol)

    if axes is None:
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(4, 8))
    else:
        ax1, ax2 = axes

    defaults = dict(square=True, linewidths=0.01, linecolor="k", yticklabels=True, cbar=False)
    sns.heatmap(data=df.iloc[:, :7], cmap=["#d3d3d3", "#450457", "#f8e621"], ax=ax1, **defaults)
    sns.heatmap(data=df.iloc[:, 7:], cmap=["#450457", "#ff7800", "#f8e621"], ax=ax2, **defaults)
    _cleanup_xaxis(ax1), _cleanup_yaxis(ax1, this_study)
    _cleanup_xaxis(ax2), _cleanup_yaxis(ax2, this_study)
    _add_legend(ax2)


def plot_lit_evidence_zscore_soma_profile(
    gene_metadata: str,
    lit_evidence: str,
    zscore_by_cluster_rep: str,
    soma_clusters: list,
    axes: list = None,
):
    """Plot heatmap of evidence and expression patterns from the literature.

    Most of the evidence patterns are protein based.

    Example
    -------
    >>> from larval_gonad.config import read_config
    >>> config = read_config("config/common.yaml")
    >>> gene_metadata = f"references/gene_annotation_dmel_{config['tag']}.feather"
    >>> lit_evidence = "data/external/miriam/lit_gene_dummy_vars.tsv"
    >>> zscore_by_cluster_rep = "output/seurat3-cluster-wf/zscore_by_cluster_rep.feather"
    >>> soma_clusters = config["soma"]
    >>> plot_lit_expression_soma_profile(gene_metadata, lit_evidence)

    """
    fbgn2symbol = get_fbgn2symbol(gene_metadata)
    soma_evidence = _get_lit_evidence(lit_evidence)[["C", "EC", "MC", "LC", "PC", "TE"]]
    this_study = list(map(lambda x: fbgn2symbol[x], _genes_w_ptraps(lit_evidence)))
    zscore_expression = _get_zscore_expression(zscore_by_cluster_rep).loc[:, (soma_clusters, slice(None))]
    df = soma_evidence.join(zscore_expression).rename(fbgn2symbol)

    if axes is None:
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(4, 8))
    else:
        ax1, ax2 = axes

    defaults = dict(square=True, linewidths=0.01, linecolor="k", yticklabels=True, xticklabels=True, cbar=False)
    sns.heatmap(data=df.iloc[:, 6:], cmap='viridis', vmin=-3, vmax=3, ax=ax1, **defaults)
    sns.heatmap(data=df.iloc[:, :6], cmap=["#d3d3d3", "#450457", "#f8e621"], ax=ax2, **defaults)
    _cleanup_xaxis_rep(ax1, soma_clusters), _cleanup_yaxis(ax1, this_study)
    _cleanup_xaxis(ax2), _cleanup_yaxis(ax2, this_study)
    _add_legend(ax2)


def _get_lit_evidence(lit_evidence):
    return pd.read_csv(lit_evidence, sep="\t", index_col=0).drop("References", axis=1)


def _genes_w_ptraps(lit_evidence):
    return (
        pd.read_csv(lit_evidence, sep="\t", index_col=0)
        .query("References == 'This study'")
        .index.tolist()
    )


def _get_binned_expression(tpm_by_cluster):
    """Bin expression.
    {
        0 not expressed: TPM < 1,
        1 low expression: 1 ≤ TPM < 5,
        2 expressed: 5 ≤ TPM < Inf
    }
    """
    return feather_to_cluster_matrix(tpm_by_cluster).apply(
        lambda x: pd.cut(x, [0, 1, 5, np.inf], labels=[0, 1, 2], right=False, include_lowest=True),
        axis=1,
    )


def _get_zscore_expression(zscore_by_cluster_rep):
    return pd.read_feather(zscore_by_cluster_rep).set_index(["FBgn", "cluster", "rep"]).squeeze().unstack([-2, -1])

def _cleanup_xaxis(ax):
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    return ax


def _cleanup_yaxis(ax, this_study):
    ax.set_ylabel("")
    labels = []
    for l in ax.get_yticklabels():
        l.set(fontstyle="italic")
        if l.get_text() in this_study:
            l.set(fontweight="bold")
        labels.append(l)
    ax.set_yticklabels(labels)
    return ax


def _add_legend(ax):
    off = mpatches.Patch(color="#450457", label="absent")
    # low = mpatches.Patch(color="#ff7800", label="low expression")
    high = mpatches.Patch(color="#f8e621", label="present")
    none = mpatches.Patch(color="#d3d3d3", label="not analyzed")
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1], handles=[off, high, none])
    return ax
