import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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
    clusters: str,
    zscores: str,
    lit_genes: dict,
    cluster_annot: dict,
    cluster_order: list,
    ax: plt.Axes = None,
    cax: plt.Axes = None,
):
    df = _read_biomarkers(biomarkers, cluster_annot)
    df_unique = _unique_biomarkers(df)

    if ax is None:
        ax, cax = _make_panel()

    _make_heatmap(zscore, ax, cax)

    _cleanup_xaxis(ax)
    _cleanup_yaxis(ax, biomarkers)
    _add_annotations(ax, biomarkers, lit_genes, zscore.shape[1])


def plot_multi_biomarkers(
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


def _make_panel():
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(2, 4))
    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    return ax, cax


def _make_heatmap(df, ax, cax=None, **kwargs):
    defaults = dict(
        df,
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


def _read_biomarkers(biomarkers, cluster_annot):
    """Read biomarkers feather file."""
    return (
        pd.read_feather(biomarkers)
        .sort_values("FBgn")
        .assign(cluster=lambda x: x.cluster.map(cluster_annot))
    )


def _unique_biomarkers(df):
    """Remove duplicated biomarkers"""
    return (
        df[~df.duplicated(subset="FBgn", keep=False)]
        .sort_values("cluster")
        .set_index("FBgn")
    )


def _multi_biomarkers(df):
    """Keep only duplicated biomarkers"""
    return (
        df[~df.duplicated(subset="FBgn", keep=False)]
        .sort_values("cluster")
        .set_index("FBgn")
    )


def _collapse_cluster(df):
    """Collapse multi biomakers into groups."""
    return (
        df.groupby("FBgn")
        .apply(lambda x: "|".join(x.cluster.sort_values()))
        .apply(_group_logic)
        .rename("group")
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


def _cleanup_xaxis(ax):
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(
        list(
            chain.from_iterable([("", x, "") for x in snakemake.params.cluster_order])
        ),
        ha="center",
        va="bottom",
    )

    # Add lines separating cell types
    for i in range(1, len(snakemake.params.cluster_order)):
        ax.axvline(i * 3, color="w", ls="--", lw=0.5)

    return ax


def _cleanup_yaxis(ax, df):
    ax.set_ylabel("")

    # Add lines separating biomarker groups
    loc = 0
    for _, dd in biomarkers.groupby("cluster"):
        loc += dd.shape[0]
        ax.axhline(loc, color="w", ls="--", lw=0.5)

    return ax


def _add_annotations(ax, biomarkers, lit_genes, cols=30):
    """Adds annotation on the left and right sides.

    Annotation includes the cluster name, the number of genes, and the gene
    names of literature genes found in that set.

    """

    loc = 0
    xloc_odd, xloc_even = (-5, cols + cols * 0.01)
    for i, (clus, dd) in enumerate(biomarkers.groupby("cluster")):
        loc, mid = take_step(loc, dd.shape[0])
        txt = f"{clus} ({dd.shape[0]:,})"
        txt += check_for_lit_genes(dd, lit_genes)

        if i % 2 == 0:
            xloc = xloc_even
        else:
            xloc = xloc_odd

        ax.text(xloc, mid, txt, ha="left", va="center", fontsize=8)


def _take_step(loc, step_size):
    prev = loc
    loc += step_size
    mid = loc - ((loc - prev) / 2)
    return loc, mid


def _check_for_lit_genes(df, lit_genes):
    res = []
    for fbgn, dd in df.sort_values("gene_symbol").iterrows():
        if fbgn in lit_genes:
            res.append(dd.gene_symbol)
    return "\n" + "\n".join(res)
