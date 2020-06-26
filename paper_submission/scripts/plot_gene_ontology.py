from typing import List, Dict

import numpy as np
import pandas as pd
import joblib
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

import larval_gonad.plotting
from larval_gonad.gene_ontology import (
    run_flyslim,
    flybase_molecular_function,
    flybase_biological_process,
    flybase_cellular_component,
)

plt.style.use(["1c", "science_base"])

def main():
    res, res_mapper = go_analysis()

    df = pd.concat(
        [
            molecular_function(res, res_mapper),
            biological_process(res, res_mapper),
            cellular_component(res, res_mapper),
        ],
        sort=False,
        ignore_index=True,
    )
    df

    g = sns.FacetGrid(df, row="asset", sharex=False, aspect=1.5)
    g.map_dataframe(
        draw_heatmap,
        vmin=-np.log10(0.05),
        vmax=-np.log10(1e-11),
        square=True,
        cmap="Blues",
        linewidths=0.2,
        linecolor="gray",
        yticklabels=False,
        xticklabels=True,
        cbar=False,
    )
    tweak_axes(g)
    add_colorbar(g, vmin=-np.log10(0.05), vmax=-np.log10(1e-11), cmap="Blues")

    g.savefig(snakemake.output.fig)


def _save_go(res: List[GOEnrichmentStudyNS]):
    df = pd.DataFrame(
        [(r.GO, r.NS, r.goterm.name, r.p_fdr_bh) for r in res],
        columns=["GO", "Asset", "Title", "FDR Adj P-value"],
    ).sort_values(["Asset", "GO"])

    df.to_csv(snakemake.output.summary, sep="\t", index=False)


def _load_gene_set(file_name: str):
    if file_name.endswith("pkl"):
        return joblib.load(file_name)

    if file_name.endswith("tsv"):
        return pd.read_csv(file_name, sep="\t").FBgn.unique().tolist()


def go_analysis():
    background = _load_gene_set(snakemake.input.background)
    target_genes = _load_gene_set(snakemake.input.target_genes)
    _, res_obj = run_flyslim(target_genes, background, cutoff=0.05, return_obj=True)
    res = res_obj.run_study(target_genes)
    res_mapper = {r.GO: r.p_fdr_bh for r in res}
    _save_go(res)
    return res, res_mapper


def draw_heatmap(*args, **kwargs):
    df = kwargs.pop("data").drop("asset", axis=1).set_index("title").T
    sns.heatmap(df, **kwargs)


def _add_other(df: pd.DataFrame, res: List[GOEnrichmentStudyNS], asset: str):
    other_p = 0
    for r in res:
        if r.NS == asset and r.p_fdr_bh <= 0.05 and r.GO not in df.index:
            other_p = max(-np.log10(r.p_fdr_bh), other_p)
    df.loc["GO:999999", "title"] = "other"
    df.loc["GO:999999", "qval"] = other_p
    return df


def molecular_function(res: List[GOEnrichmentStudyNS], res_mapper: Dict[str, float]):
    order = list(flybase_molecular_function.values()) + ["other"]
    return (
        pd.Series(flybase_molecular_function, name="title")
        .rename_axis("GO")
        .to_frame()
        .assign(qval=lambda sr: -np.log10(sr.index.map(res_mapper)))
        .pipe(_add_other, res=res, asset="MF")
        .groupby("title")
        .qval.max()
        .to_frame()
        .assign(asset="Molecular Function")
        .reindex(order)
        .reset_index()
    )


def biological_process(res: List[GOEnrichmentStudyNS], res_mapper: Dict[str, float]):
    order = list(flybase_biological_process.values()) + ["other"]
    return (
        pd.Series(flybase_biological_process, name="title")
        .rename_axis("GO")
        .to_frame()
        .assign(qval=lambda sr: -np.log10(sr.index.map(res_mapper)))
        .pipe(_add_other, res=res, asset="BP")
        .groupby("title")
        .qval.max()
        .to_frame()
        .assign(asset="Biological Process")
        .reindex(order)
        .reset_index()
    )


def cellular_component(res: List[GOEnrichmentStudyNS], res_mapper: Dict[str, float]):
    order = list(flybase_cellular_component.values()) + ["other"]
    return (
        pd.Series(flybase_cellular_component, name="title")
        .rename_axis("GO")
        .to_frame()
        .assign(qval=lambda sr: -np.log10(sr.index.map(res_mapper)))
        .pipe(_add_other, res=res, asset="CC")
        .groupby("title")
        .qval.max()
        .to_frame()
        .assign(asset="Cellular Component")
        .reindex(order)
        .reset_index()
    )


def _tweak_ax(ax: plt.Axes, xlabel: str):
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="left")
    ax.set(xlabel=xlabel)


def tweak_axes(g: sns.FacetGrid):
    g.set_titles("")
    for label, ax in zip(g.row_names, g.axes.ravel()):
        _tweak_ax(ax, label)
    return g


def add_colorbar(g: sns.FacetGrid):
    g.fig.subplots_adjust(bottom=0.2)
    ax = g.fig.add_subplot()


def add_colorbar(g: sns.FacetGrid, label="-Log10(P-Value)", **kwargs):
    g.fig.subplots_adjust(bottom=0.08)
    cax = g.fig.add_axes([0.20, 0.10, 0.60, 0.02])
    points = plt.scatter([], [], c=[], **kwargs)
    g.fig.colorbar(points, cax=cax, label=label, orientation="horizontal")
    return g


if __name__ == "__main__":
    main()
