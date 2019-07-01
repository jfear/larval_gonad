"""Functions for use with X To A Analysis."""
import os
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, scoreatpercentile
import matplotlib.pyplot as plt
import seaborn as sns

from pysam import AlignmentFile

from .config import memory, config
from .plotting import figure_element
from .scRNAseq import norm_data, raw_data, seurat_or_data

CHROMS = ["X", "2L", "2R", "3L", "3R", "4", "Y"]
CHROMS_CHR = ["chrX", "chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrY"]

AUTOSOMES = ["2L", "2R", "3L", "3R", "4"]
AUTOSOMES_CHR = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4"]

MAJOR_ARMS = ["2L", "2R", "3L", "3R"]
MAJOR_ARMS_CHR = ["chr2L", "chr2R", "chr3L", "chr3R"]

fbgn2chrom = pd.read_csv("../output/fbgn2chrom.tsv", sep="\t", index_col=0)


# Data Munging Functions
def agg_all(kind, seurat_dir, cluster=None, resolution=None):
    try:
        assert kind in ["raw", "norm"]
    except AssertionError as e:
        print('kind must be "raw" or "norm"')
        raise e

    if kind == "raw":
        df = raw_data(seurat_dir, cluster=cluster, resolution=resolution).sum(axis=1)
    elif kind == "norm":
        df = norm_data(seurat_dir, cluster=cluster, resolution=resolution).sum(axis=1)

    df.name = "cnts"
    return df


@seurat_or_data
def x_to_a(cluster, data=None, seurat_dir=None):
    """Calculates the X:A ratio of medians for a given cluster."""
    if data is None:
        cnts = agg_all("norm", seurat_dir=seurat_dir, cluster=cluster)
        data = cnts.to_frame().join(fbgn2chrom)

    med_by_chrom = data.groupby("chrom").median().loc[CHROMS_CHR[:-1]]
    autosome_median = data[data.chrom.isin(MAJOR_ARMS_CHR)].cnts.median()
    return med_by_chrom / autosome_median


@seurat_or_data
def commonly_expressed(data=None, seurat_dir=None, read_cutoff=0):
    """Create list of genes expressed in 1/3 of cells."""

    if data is None:
        data = norm_data(seurat_dir)

    mask = (data > read_cutoff).sum(axis=1) > (data.shape[1] / 3)
    return data.index[mask].tolist()


def clean_pvalue(pval, use_text=True):
    if (pval < 0.0001) & use_text:
        pvalue = "P<0.0001"
    elif pval < 0.0001:
        pvalue = "***"
    elif (pval < 0.001) & use_text:
        pvalue = "P<0.001"
    elif pval < 0.001:
        pvalue = "**"
    elif (pval < 0.01) & use_text:
        pvalue = "P<0.01"
    elif pval < 0.01:
        pvalue = "**"
    else:
        pvalue = "N.S."

    return pvalue


def strip_chr(ax):
    labels = ax.get_xticklabels()
    new = []
    for l in labels:
        new.append(l.get_text().strip("chr"))
    ax.set_xticklabels(new)


@memory.cache
def idx_stats_by_cluster(bam: str, cluster: str) -> pd.DataFrame:
    """Function to count the number of reads by cluster by chromosome.

    Parameters
    ----------
    bam : str
        The path to a bam file.
    cluster : str
        The path to a cluster file.

    Returns
    -------
    pd.DataFrame:
        DataFrame where rows are chromosome, columns are clusters, and values
        are number of aligned reads.

    """
    cluster = pd.read_csv(cluster, sep="\t")
    lookup = cluster.ident.to_dict()

    dat = AlignmentFile(bam, "rb")
    results = defaultdict(lambda: defaultdict(int))

    for read in dat.fetch():
        if read.is_unmapped:
            continue

        chrom = read.reference_name
        if chrom not in CHROMS:
            continue

        try:
            cell = read.get_tag("CB").replace("-1", "")
            clus = lookup[cell]
            results[clus][chrom] += 1
        except KeyError:
            pass

    df = pd.DataFrame(results)
    df.index.name = "chrom"
    return df


def x_autosome_boxplot(x, y, pvalue_cutoff=0.001, x_chrom=None, autosomes=None, **kwargs):
    """Boxplot to compare X to Autsome"""
    if autosomes is None:
        autosomes = AUTOSOMES_CHR

    if x_chrom is None:
        x_chrom = "chrX"

    _dat = kwargs.pop("data")
    chrX = _dat[_dat[x] == x_chrom][y]
    chrA = _dat[_dat[x].isin(autosomes)][y]
    stat, pvalue = mannwhitneyu(chrX, chrA, alternative="two-sided")

    med_x = chrX.median()
    med_a = chrA.median()
    diff = round(med_a - med_x, 2)

    # Pepare for plotting
    chrX = chrX.to_frame()
    chrX["class"] = "X"

    chrA = chrA.to_frame()
    chrA["class"] = "A"

    # Plot
    ax = sns.boxplot("class", y, data=pd.concat([chrX, chrA]), order=["X", "A"], **kwargs)
    ax.axhline(med_a, c=config["colors"]["c2"], ls="--", label="Median Autsome Expression")
    ax.axhline(med_x, c=config["colors"]["c3"], ls="--", label="Median X Expression")
    ax.legend()
    ax.text(0.5, -0.2, f"A - X = {diff}", transform=ax.transAxes, ha="center")

    # Add p-value
    if pvalue <= pvalue_cutoff:
        pvalue = clean_pvalue(pvalue)
        x1, x2 = 0, 1
        iqr = sns.utils.iqr(_dat[y])
        y, h, col = iqr + iqr * 2, 0.6, "k"
        plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=col)
        plt.text((x1 + x2) * 0.5, y + h + 0.01, f"{pvalue}", ha="center", va="bottom", color=col)

    return ax


def estimate_dcc(chrom_col, count_col, df):
    """Estimate Dosage Compensation.

    Parameters
    ----------
    chrom_col: The column name with chromosomes info.
    count_col: The column name with the count info.
    df: Dataframe with chrom_col and count_col.

    Returns
    -------
    median of X, median of major autosomes, proporition compensation.

    """
    med_major = df.loc[df[chrom_col].isin(MAJOR_ARMS_CHR), count_col].median()
    med_x = df.loc[df[chrom_col] == "chrX", count_col].median()

    try:
        prop_dcc = np.round((med_x / med_major) * 2, 2)
    except ZeroDivisionError:
        prop_dcc = np.nan

    return med_x, med_major, prop_dcc


def multi_chrom_boxplot(x, y, use_text=True, multiplier=2, ax=None, **kwargs):
    """Designed to be used with Seaborn.FacetGrid"""
    _dat = kwargs["data"]
    _dat = _dat[_dat[x].isin(CHROMS_CHR)]
    meds = _dat.groupby(x).median()

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    sns.boxplot(x, y, order=CHROMS_CHR[:-1], ax=ax, **kwargs)
    ax.scatter(
        range(meds.shape[0]), meds.loc[CHROMS_CHR[:-1]].values, color="w", edgecolor="k", marker="d"
    )

    med_x, med_major, prop_dcc = estimate_dcc(x, y, _dat)
    ax.axhline(
        med_major,
        ls="--",
        color=sns.xkcd_rgb["coral"],
        lw=2,
        label="Median Autosomal Expression",
        zorder=3,
    )

    # Clean up the pvalue for plotting
    pvalues = {}
    iqr = 0
    chromX = _dat[_dat[x] == "chrX"]
    for g, df in _dat.groupby(x):
        _iqr = sns.utils.iqr(df[y])
        if _iqr > iqr:
            iqr = _iqr
        if g == "chrX":
            continue
        _, pval = mannwhitneyu(chromX[y], df[y], alternative="two-sided")
        if pval <= 0.001:
            pvalues[g] = pval

    if isinstance(multiplier, tuple):
        multiplier, increment = multiplier
    else:
        increment = 1

    xloc = CHROMS_CHR.index("chrX")
    for k, v in pvalues.items():
        oloc = CHROMS_CHR.index(k)
        pval = clean_pvalue(v, use_text)
        y, h, col = iqr + iqr * multiplier, 0.4, "k"
        ax.plot([xloc, xloc, oloc, oloc], [y, y + h, y + h, y], lw=1, c=col)
        ax.text(
            (xloc + oloc) * 0.5,
            y + h + 0.01,
            f"{pval}",
            ha="center",
            va="bottom",
            color=col,
            fontsize=22,
        )
        multiplier += increment

    try:
        curr_cluster = _dat["cluster"].values[0]
        num_cells = kwargs.pop("num_cells")[curr_cluster]
        num_genes = _dat.shape[0]

        ax.text(
            0.05,
            0.85,
            (f"X Compensation: {prop_dcc}\n" f"(n={num_cells} cells)\n" f"(g={num_genes:,} genes)"),
            transform=ax.transAxes,
        )
    except KeyError:
        pass


def multi_chrom_boxplot2(x, y, use_text=True, multiplier=2, ax=None, **kwargs):
    """Designed to be used with Seaborn.FacetGrid"""
    _dat = kwargs["data"]
    _dat = _dat[_dat[x].isin(CHROMS_CHR)]
    meds = _dat.groupby(x).median()

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    sns.boxplot(x, y, order=CHROMS_CHR[:-1], ax=ax, **kwargs)
    ax.scatter(
        range(meds.shape[0]), meds.loc[CHROMS_CHR[:-1]].values, color="w", edgecolor="k", marker="d"
    )

    med_x, med_major, prop_dcc = estimate_dcc(x, y, _dat)
    ax.axhline(
        med_major,
        ls="--",
        color=sns.xkcd_rgb["coral"],
        lw=2,
        label="Median Autosomal Expression",
        zorder=3,
    )

    # Clean up the pvalue for plotting
    pvaluesX = False
    pvalues4 = False
    chromX = _dat[_dat[x] == "chrX"]
    chrom4 = _dat[_dat[x] == "chr4"]
    for g, df in _dat.groupby(x):
        if g == "chrX":
            q3X = scoreatpercentile(np.asarray(df[y]), 75)
            iqrX = sns.utils.iqr(df[y])
            continue

        if g == "chr4":
            q34 = scoreatpercentile(np.asarray(df[y]), 75)
            iqr4 = sns.utils.iqr(df[y])
            continue

        _, pval = mannwhitneyu(chromX[y], df[y], alternative="two-sided")
        if pval <= 0.01:
            pvaluesX = True

        _, pval = mannwhitneyu(chrom4[y], df[y], alternative="two-sided")
        if pval <= 0.01:
            pvalues4 = True

    xloc = CHROMS_CHR.index("chrX")
    if pvaluesX:
        ax.text(xloc, q3X + 1.51 * iqrX, "*", ha="center", fontsize=12, fontweight="bold")

    _4loc = CHROMS_CHR.index("chr4")
    if pvalues4:
        ax.text(_4loc, q34 + 1.51 * iqr4, "*", ha="center", fontsize=12, fontweight="bold")


def plot_cluster_x2a(dat, fbgns, cluster_id, ax1, ax2, fbgn2chrom):
    reds = sns.color_palette("Reds")
    boxplot_colors = [
        reds[-1],  # X
        "#ffffff",  # 2L
        "#ffffff",  # 2R
        "#ffffff",  # 3L
        "#ffffff",  # 3R
        reds[-1],  # 4
    ]

    idx = dat.query(f"cluster == {cluster_id}").index
    dat.drop("cluster", axis=1)

    dataM = dat.loc[idx, fbgns].median().to_frame().join(fbgn2chrom).query('chrom != "chrY"')
    dataM.columns = ["Normalized Expression (Median)", "Chromosome Arm"]

    dataS = dat.loc[idx, fbgns].sum().to_frame().join(fbgn2chrom).query('chrom != "chrY"')

    dataS.columns = ["Normalized Expression (Sum)", "Chromosome Arm"]

    multi_chrom_boxplot2(
        "Chromosome Arm",
        "Normalized Expression (Sum)",
        data=dataS,
        notch=True,
        palette=boxplot_colors,
        use_text=False,
        multiplier=(1, 0.1),
        showfliers=False,
        ax=ax1,
    )

    multi_chrom_boxplot2(
        "Chromosome Arm",
        "Normalized Expression (Median)",
        data=dataM,
        notch=True,
        palette=boxplot_colors,
        use_text=False,
        multiplier=(1, 0),
        showfliers=False,
        ax=ax2,
    )

    strip_chr(ax1)
    strip_chr(ax2)

    ax1.set_xlabel("")
    ax2.set_xlabel("")

    ax1.set_ylabel("")
    ax2.set_ylabel("")


def _get_mapper():
    REF_DIR = os.environ["REFERENCES_DIR"]
    assembly = config["assembly"]
    tag = config["tag"]
    fname = Path(REF_DIR, assembly, tag, "fb_annotation", f"dmel_{tag}.fb_annotation")

    annot = pd.read_csv(fname, sep="\t")
    mapper = {}
    for i, row in annot.iterrows():
        fbgn = row.primary_FBgn
        mapper[fbgn] = fbgn

        if isinstance(row.annotation_ID, str):
            mapper[row.annotation_ID] = fbgn

        if isinstance(row.secondary_FBgn, str):
            for sec in row.secondary_FBgn.split(","):
                mapper[sec] = fbgn

        if isinstance(row.secondary_annotation_ID, str):
            for acn in row.secondary_annotation_ID.split(","):
                mapper[acn] = fbgn

    return mapper


MAPPER = _get_mapper()


def cleanup_FBgn(dat):
    """Convert FBgn or Accn to current FBgn."""
    res = []
    for k in dat:
        try:
            res.append(MAPPER[k])
        except KeyError:
            if isinstance(k, str):
                print(f"{k} not found in current annotation.")
    return np.array(res)


def get_gene_sets():
    PROJ_DIR = "../"

    genes = {}
    # Tau
    tau = pd.read_csv(
        Path(PROJ_DIR, "output/2018-02-05_tau_haiwang_male_tau.tsv"),
        sep="\t",
        index_col="FBgn",
        squeeze=True,
    ).dropna()
    tau_genes = cleanup_FBgn(tau[tau <= 0.4].index)
    genes["Haiwang_male_tau"] = tau_genes

    # TSPS
    tsps = pd.read_csv(
        Path(PROJ_DIR, "output/2018-02-05_tau_haiwang_male_tsps.tsv"),
        sep="\t",
        index_col="FBgn",
        squeeze=True,
    ).dropna()
    tsps_genes = cleanup_FBgn(tsps[tsps < 1.0].index)
    genes["Haiwang_male_tsps"] = tsps_genes

    # Naieve Bayes
    with open(Path(PROJ_DIR, "data/external/Ferrari_et_al_2006_housekeeping_FBgn.txt")) as fh:
        housekeeping = np.array(fh.read().splitlines())

    bayes_genes = cleanup_FBgn(housekeeping)
    genes["Ferrari_housekeeping"] = bayes_genes

    # Protein groups from DroID and DPiM
    protein = pd.read_csv(
        Path(PROJ_DIR, "data/external/DroID_DPiM_2018-03-29.txt"), sep="\t", low_memory=False
    )

    network = defaultdict(set)
    for i, (bait, inter) in protein[["FBGN_BAIT", "FBGN_INTERACTOR"]].iterrows():
        network[bait].add(inter)
        network[inter].add(bait)

    # Sanity check to make sure I only have genes with at least one interaction
    dpim = set()
    for k, v in network.items():
        if len(v) > 1:
            dpim.add(k)
            dpim.union(set(v))

    dpim_genes = cleanup_FBgn(dpim)
    genes["DPiM_protein_complex"] = dpim_genes

    # Parse FlyBase Gene Groups
    header = [
        "FB_group_id",
        "FB_group_symbol",
        "FB_group_name",
        "Parent_FB_group_id",
        "Parent_FB_group_symbol",
        "FBgn",
        "gene_symbol",
    ]

    # Genes in any gene group
    fb_groups = pd.read_csv(
        Path(PROJ_DIR, "data/external/gene_group_data_fb_2017_03.tsv"),
        sep="\t",
        comment="#",
        names=header,
    )
    genes["flybase_groups"] = {}
    for i, group in fb_groups.groupby("FB_group_name"):
        fb_genes = cleanup_FBgn(group.FBgn)
        genes["flybase_groups"][i] = fb_genes

    return genes
