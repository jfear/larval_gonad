"""Common elements used for the scRNA-Seq data.

This is a place to store functions and variables that are used over and over
again with the scRNA-seq data.

"""
import os
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from .config import memory, config
from .plotting import figure_element


@memory.cache
def read_counts(fname):
    df = pd.read_csv(fname, sep="\t")
    df.index.name = "FBgn"
    df.columns.name = "cell_id"
    return df


class Seurat(object):
    """Class that stores basic paths for files from Seurat analysis."""

    def __init__(self, path=None):
        """Create a list of Seurat paths.
        raw : str
            Path to raw.
        metadata : str
            Path to metadata.
        scaled : str
            Path to scaled data.
        dispersion : str
            Path to dispersion estimates.
        var_genes : str
            Path to variable gene list.
        normalized_read_counts : str
            Path to normalized read counts.
        principal_components_cell : str
            Path to PCA cell loadings.
        principal_components_gene : str
            Path to PCA gene loadings.
        tsne : str
            Path to tsne.
        biomarkers : str
            Path to biomarkers list.
        clusters : str
            Path to clusters identities.
        robj : str
            Path to seurat R object.

        """

        self.path = path

        if path is None:
            self.raw = None
            self.metadata = None
            self.scaled = None
            self.dispersion = None
            self.var_genes = None
            self.normalized_read_counts = None
            self.principal_components_cell = None
            self.principal_components_gene = None
            self.principal_components_stdev = None
            self.tsne = None
            self.biomarkers = None
            self.clusters = None
            self.robj = None
        else:
            self.raw = os.path.join(path, "raw.tsv")
            self.metadata = os.path.join(path, "metadata.tsv")
            self.scaled = os.path.join(path, "scaled.tsv")
            self.dispersion = os.path.join(path, "dispersion.tsv")
            self.var_genes = os.path.join(path, "var_genes.txt")
            self.normalized_read_counts = os.path.join(path, "normalized_read_counts.tsv")
            self.principal_components_cell = os.path.join(path, "principal_components_cell.tsv")
            self.principal_components_gene = os.path.join(path, "principal_components_gene.tsv")
            self.principal_components_stdev = os.path.join(path, "principal_components_stdev.tsv")
            self.tsne = os.path.join(path, "tsne.tsv")
            self.biomarkers = os.path.join(path, "biomarkers.tsv")
            self.clusters = os.path.join(path, "clusters.tsv")
            self.robj = os.path.join(path, "seurat.Robj")

    def get_raw(self):
        return read_counts(self.raw)

    def get_metadata(self):
        df = pd.read_csv(self.metadata, sep="\t")
        df.index.name = "cell_id"
        return df

    def get_scaled(self):
        df = pd.read_csv(self.scaled, sep="\t")
        df.index.name = "FBgn"
        df.columns.name = "cell_id"
        return df

    def get_dispersion(self):
        df = pd.read_csv(self.dispersion, sep="\t")
        df.index.name = "FBgn"
        return df

    def get_var_genes(self):
        with open(self.var_genes) as fh:
            return [x.strip() for x in fh.readlines()]

    def get_normalized_read_counts(self):
        return read_counts(self.normalized_read_counts)

    def get_principal_components_cell(self):
        df = pd.read_csv(self.principal_components_cell, sep="\t")
        df.index.name = "cell_id"
        df.columns.name = "PC"
        return df

    def get_principal_components_gene(self):
        df = pd.read_csv(self.principal_components_gene, sep="\t")
        df.index.name = "FBgn"
        df.columns.name = "PC"
        return df

    def get_principal_components_stdev(self):
        df = pd.read_csv(self.principal_components_stdev, sep="\t")
        df.index.name = "PC"
        df.columns.name = "stdev"
        return df

    def get_tsne(self):
        df = pd.read_csv(self.tsne, sep="\t")
        df.index.name = "FBgn"
        return df

    def get_biomarkers(self, resolution):
        fname = Path(self.path, f"biomarkers_{resolution}.tsv")
        df = pd.read_csv(fname, sep="\t", index_col="primary_FBgn")
        df.index.name = "FBgn"
        return df

    def get_clusters(self, resolution=None):
        df = pd.read_csv(self.clusters, sep="\t", index_col=0)

        if resolution is None:
            return df

        clusters = df[resolution]
        clusters.name = "cluster"
        return clusters


# Data Munging Functions
@memory.cache
def raw_data(seurat_dir, cluster=None, resolution=None):
    s = Seurat(seurat_dir)

    if cluster is None:
        return s.get_raw()

    if resolution is None:
        resolution = config["resolution"]

    clusters = s.get_clusters(resolution)
    cells = clusters.index[clusters == cluster].tolist()
    return s.get_raw()[cells]


@memory.cache
def norm_data(seurat_dir, cluster=None, resolution=None):
    s = Seurat(seurat_dir)

    if cluster is None:
        return s.get_normalized_read_counts()

    if resolution is None:
        resolution = config["resolution"]

    clusters = s.get_clusters(resolution)
    cells = clusters.index[clusters == cluster].tolist()
    return s.get_normalized_read_counts()[cells]


def seurat_or_data(func):
    def wrapper(*args, **kwargs):
        data = kwargs.get("data", None)
        seurat_dir = kwargs.get("seurat_dir", None)
        if (data is None) & (seurat_dir is None):
            raise ValueError("You must provide either data or a seurat_dir")

        return func(*args, **kwargs)

    return wrapper


# Plotting Functions
def TSNEPlot(
    x="tSNE_1",
    y="tSNE_2",
    data=None,
    hue=None,
    cmap=None,
    palette=None,
    ax=None,
    class_names=None,
    legend_kws=None,
    cbar=True,
    **kwargs,
):
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
    defaults = {"vmin": 0, "vmax": 5, "edgecolor": "k"}
    defaults.update(kwargs)

    df = data.copy()
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    if palette is None:
        palette = sns.color_palette()

    legend_defaults = {"loc": "center left", "bbox_to_anchor": (1, 0.5)}
    if legend_kws is not None:
        legend_defaults.update(legend_kws)

    if isinstance(hue, pd.Series):
        df["on"] = hue.astype(int).apply(lambda x: str(x))
        hue = "on"

    values = sorted(df[hue].unique())
    if (isinstance(values[0], (int, np.integer)) & (len(values) > 20)) | isinstance(
        values[0], (float, np.float64)
    ):
        if cmap is None:
            cmap = mpl.colors.ListedColormap(palette)
        zeros = df[df[hue] == 0]
        if len(zeros) > 0:
            zeros.plot.scatter(
                x, y, c=zeros[hue], cmap=cmap, ax=ax, colorbar=False, zorder=1, **defaults
            )

        expressed = df[df[hue] > 0]
        if len(expressed) > 0:
            expressed.plot.scatter(
                x, y, c=expressed[hue], cmap=cmap, ax=ax, colorbar=cbar, zorder=3, **defaults
            )
    else:
        if cmap is None:
            cmap = {k: v for k, v in zip(values, palette)}

        for l, dd in df.groupby(hue):
            if class_names is None:
                _class = str(l).title()
            else:
                _class = class_names[l]

            try:
                dd.plot.scatter(x, y, c=cmap[l], label=_class, ax=ax, **defaults)
            except KeyError as e:
                print("Try Setting Palette with the correct number of colors.")
                print(len(cmap), len(values))
                print(type(values[0]))
                raise e

        ax.legend(**legend_defaults)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_yticks([])


def plot_confusion_matrix(
    cm, classes, normalize=False, title="Confusion matrix", cmap=plt.cm.Blues
):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print("Confusion matrix, without normalization")

    plt.imshow(cm, interpolation="nearest", cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45)
    plt.yticks(tick_marks, classes)

    fmt = ".2f" if normalize else "d"
    thresh = cm.max() / 2.0
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(
            j,
            i,
            format(cm[i, j], fmt),
            horizontalalignment="center",
            color="white" if cm[i, j] > thresh else "black",
        )

    plt.tight_layout()
    plt.ylabel("True label")
    plt.xlabel("Predicted label")


@seurat_or_data
@figure_element
def fe_tsne(data=None, seurat_dir=None, ax=None, resolution=None, **kwargs):
    if resolution is None:
        resolution = config["resolution"]

    if data is None:
        s = Seurat(seurat_dir)
        tsne = s.get_tsne()
        clusters = s.get_clusters(resolution)
        data = tsne.join(clusters)

    TSNEPlot(
        data=data,
        hue="cluster",
        class_names=config["cluster_annot"],
        palette=config["colors"]["clusters"],
        ax=ax,
        **kwargs,
    )

    return ax


@seurat_or_data
@figure_element
def fe_heatmap_all(data=None, seurat_dir=None, ax=None, resolution=None, **kwargs):
    if resolution is None:
        resolution = config["resolution"]

    if data is None:
        s = Seurat(seurat_dir)
        norm = s.get_normalized_read_counts()
        clusters = s.get_clusters(resolution)
        data = tsne.join(clusters)

    TSNEPlot(
        data=data,
        hue="cluster",
        class_names=config["cluster_annot"],
        palette=config["colors"]["clusters"],
        ax=ax,
        **kwargs,
    )

    return ax
