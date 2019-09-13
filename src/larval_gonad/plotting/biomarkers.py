from itertools import chain

from more_itertools import flatten
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

MULTI_ORDER = [
    "G|EPS",
    "Germline",
    "Spermatocytes",
    "Cyst Cells",
    "Germline and Somatic",
    "Somatic",
]


def plot_all_biomarkers(
    biomarkers: str,
    clusters: str,
    zscores: str,
    lit_genes: dict,
    cluster_annot: dict,
    cluster_order: list,
    ax: plt.Axes = None,
    cax: plt.Axes = None,
):
    pass


def plot_unique_biomarkers(
    biomarkers: str,
    zscores: str,
    lit_genes: dict,
    ax: plt.Axes = None,
    cax: plt.Axes = None,
    fig: plt.Figure = None,
):
    """Plot heatmap of unqiue biomarkers

    Example
    -------
    >>> biomarkers = "output/seurat3-cluster-wf/biomarkers.feather"
    >>> zscores = "output/seurat3-cluster-wf/zscore_by_cluster_rep.feather"
    >>> from larval_gonad.config import read_config
    >>> lit_genes = read_config("config/literature_genes.yaml")
    >>> plot_unique_biomarkers(biomarkers, zscores, lit_genes)

    """
    df = _read_biomarkers(biomarkers)
    df_unique = _unique_biomarkers(df)

    df_zscore = _read_zscores(zscores)
    zscore_ordered = _order_zscores(df_zscore, df_unique.index.unique().tolist())

    if ax is None:
        ax, cax = _make_panel(fig)

    _make_heatmap(zscore_ordered, ax, cax)
    _cleanup_xaxis(ax, df_unique)
    _cleanup_yaxis(ax, df_unique)
    _add_annotations(ax, df_unique, lit_genes, zscore_ordered.shape[1])


def plot_multi_biomarkers(
    biomarkers: str,
    zscores: str,
    lit_genes: dict,
    ax: plt.Axes = None,
    cax: plt.Axes = None,
    fig: plt.Figure = None,
):
    """Plot heatmap of multi biomarkers

    Example
    -------
    >>> biomarkers = "output/seurat3-cluster-wf/biomarkers.feather"
    >>> zscores = "output/seurat3-cluster-wf/zscore_by_cluster_rep.feather"
    >>> from larval_gonad.config import read_config
    >>> lit_genes = read_config("config/literature_genes.yaml")
    >>> plot_multi_biomarkers(biomarkers, zscores, lit_genes)

    """
    df = _read_biomarkers(biomarkers)
    df_multi = _multi_biomarkers(df)
    df_grouped = _order_multi_groups(_collapse_cluster(df_multi))

    df_zscore = _read_zscores(zscores)
    zscore_ordered = _order_zscores(df_zscore, df_grouped.index.unique().tolist())

    if ax is None:
        ax, cax = _make_panel(fig)

    _make_heatmap(zscore_ordered, ax, cax)
    _cleanup_xaxis(ax, df_grouped)
    _cleanup_yaxis(ax, df_grouped)
    _add_annotations(ax, df_grouped, lit_genes, zscore_ordered.shape[1])


def _add_annotations(ax, biomarkers, lit_genes, cols=30):
    """Adds annotation on the left and right sides.

    Annotation includes the cluster name, the number of genes, and the gene
    names of literature genes found in that set.

    """
    lit_fbgns = list(flatten(lit_genes.values()))

    loc = 0
    xloc_odd, xloc_even = (-5, cols + cols * 0.01)
    for i, (clus, dd) in enumerate(biomarkers.groupby("cluster")):
        loc, mid = _take_step(loc, dd.shape[0])
        txt = f"{clus} ({dd.shape[0]:,})"
        txt += _check_for_lit_genes(dd, lit_fbgns)

        if i % 2 == 0:
            xloc = xloc_even
        else:
            xloc = xloc_odd

        ax.text(xloc, mid, txt, ha="left", va="center", fontsize=8)


def _check_for_lit_genes(df, lit_fbgns):
    res = []
    for fbgn, dd in df.sort_values("gene_symbol").iterrows():
        if fbgn in lit_fbgns:
            res.append(dd.gene_symbol)
    return "\n" + "\n".join(res)


def _cleanup_xaxis(ax, biomarkers):
    cluster_order = biomarkers.cluster.cat.categories

    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(
        list(chain.from_iterable([("", x, "") for x in cluster_order])),
        ha="center",
        va="bottom",
    )

    # Add lines separating cell types
    for i in range(1, len(cluster_order)):
        ax.axvline(i * 3, color="w", ls="--", lw=0.5)

    return ax


def _cleanup_yaxis(ax, biomarkers):
    ax.set_ylabel("")

    # Add lines separating biomarker groups
    loc = 0
    for _, dd in biomarkers.groupby("cluster"):
        loc += dd.shape[0]
        ax.axhline(loc, color="w", ls="--", lw=0.5)

    return ax


def _collapse_cluster(df):
    """Collapse multi biomakers into groups."""
    return (
        df.groupby("FBgn")
        .apply(lambda x: "|".join(x.cluster.sort_values()))
        .apply(_group_logic)
        .rename("cluster")
        .to_frame()
        .join(df.gene_symbol.drop_duplicates())
    )


def _group_logic(x):
    """Label multi biomarker groups.

    Some biomarkers correspond to multiple clusters. This sets up the logic to 
    assign these clusters to individual groups.
    """
    if x == "G|EPS":
        return "G|EPS"
    elif (
        (x == "G|EPS|MPS|LPS")
        | (x == "G|EPS|MPS")
        | (x == "G|EPS|LPS")
        | (x == "G|MPS|LPS")
        | (x == "G|MPS")
        | (x == "G|LPS")
    ):
        return "Germline"
    elif (x == "EPS|MPS|LPS") | (x == "EPS|MPS") | (x == "EPS|LPS") | (x == "MPS|LPS"):
        return "Spermatocytes"
    elif (
        (x == "C1|C2|C3|C4")
        | (x == "C1|C2|C3")
        | (x == "C1|C2|C4")
        | (x == "C1|C3|C4")
        | (x == "C2|C3|C4")
        | (x == "C1|C2")
        | (x == "C1|C3")
        | (x == "C1|C4")
        | (x == "C2|C3")
        | (x == "C2|C4")
        | (x == "C3|C4")
    ):
        return "Cyst Cells"
    elif (("G" in x) | ("EPS" in x) | ("MPS" in x) | ("LPS" in x)) & (
        ("C1" in x) | ("C2" in x) | ("C3" in x) | ("P" in x) | ("T" in x)
    ):
        return "Germline and Somatic"
    else:
        return "Somatic"


def _make_panel(fig=None):
    from matplotlib.gridspec import GridSpec

    if fig is None:
        fig = plt.figure(figsize=(2, 4))

    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    return ax, cax


def _make_heatmap(df, ax, cax=None, **kwargs):
    defaults = dict(
        data=df,
        xticklabels=True,
        yticklabels=False,
        vmin=-3,
        vmax=3,
        rasterized=True,
        cmap="viridis",
        ax=ax,
        cbar_ax=cax,
        cbar_kws=dict(label="Z-Score (TPM)", ticks=[-3, 3], orientation="horizontal"),
    )
    defaults.update(kwargs)
    if df.shape[0] < 100:
        defaults["yticklabels"] = True

    return sns.heatmap(**defaults)


def _multi_biomarkers(df):
    """Keep only duplicated biomarkers"""
    return df[df.index.duplicated(keep=False)].sort_values("cluster")


def _order_multi_groups(biomarkers):
    return biomarkers.assign(
        cluster=lambda x: pd.Categorical(x.cluster, categories=MULTI_ORDER, ordered=True)
    ).sort_values("cluster")


def _order_zscores(zscores, fbgns):
    return zscores.reindex(fbgns)


def _read_biomarkers(biomarkers):
    """Read biomarkers feather file."""
    return pd.read_feather(biomarkers).set_index("FBgn").sort_index()


def _read_zscores(zscores):
    return pd.pivot_table(
        pd.read_feather(zscores),
        index="FBgn",
        columns=["cluster", "rep"],
        values="zscore",
        aggfunc="first",
    )


def _take_step(loc, step_size):
    prev = loc
    loc += step_size
    mid = loc - ((loc - prev) / 2)
    return loc, mid


def _unique_biomarkers(df):
    """Remove duplicated biomarkers"""
    return df[~df.index.duplicated(keep=False)].sort_values("cluster")

