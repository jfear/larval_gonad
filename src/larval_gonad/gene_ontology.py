from pathlib import Path
from textwrap import dedent
from collections import defaultdict

import matplotlib.pyplot as plt

from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

HERE = Path(__file__).absolute().parent

oboDag = GODag(HERE / "../../data/external/go-basic.obo")
slimDag = GODag(HERE / "../../data/external/goslim_drosophila_fb_ribbons.obo")
association = HERE / "../../data/external/gene_association.fb"

flybase_molecular_function = {
    "GO:0003824": "enzyme",
    "GO:0098772": "regulator",
    "GO:0038023": "receptor",
    # "GO:0005102": "receptor binding",
    "GO:0005215": "transporter",
    "GO:0005198": "structural molecule",
    "GO:0008092": "cytoskeleton binding",
    "GO:0003677": "DNA binding",
    "GO:0003723": "RNA binding",
    "GO:0140110": "transcription factor",
    "GO:0036094": "small molecule binding",
    "GO:0046872": "metal ion binding",
    # "GO:0008289": "lipid binding",
    # "GO:0030246": "carbohydrate binding",
}

flybase_biological_process = {
    "GO:0008283": "cell cycle/proliferation",
    "GO:0007049": "cell cycle/proliferation",
    "GO:0071840": "cellular organization/biogenesis",
    "GO:0051234": "cellular transport/localization",
    "GO:0033036": "cellular transport/localization",
    "GO:0032502": "development",
    # "GO:0000003": "reproduction",
    "GO:0048232": "male gamete generation",
    "GO:0007292": "female gamete generation",
    "GO:0002376": "immune system",
    # "GO:0050877": "nervous system process",
    # "GO:0007610": "behavior",
    "GO:0050896": "response to stimulus",
    "GO:0023052": "signaling",
    # "GO:0010467": "gene expression",
    "GO:0019538": "protein metabolism",
    "GO:0006259": "DNA metabolism",
    "GO:0044281": "small molecule metabolism",
}

flybase_cellular_component = {
    "GO:0005576": "extracellular",
    "GO:0005829": "cytosol",
    "GO:0005856": "cytoskeleton",
    "GO:0005739": "mitochondrion",
    "GO:0005634": "nucleus",
    "GO:0005694": "chromosome",
    "GO:0016020": "membrane",
    # "GO:0071944": "cell periphery",
    "GO:0042995": "cell projection",
    "GO:0005773": "endomembrane system",
    "GO:0012505": "endomembrane system",
    "GO:0045202": "synapse",
    "GO:0032991": "macromolecular complex",
}

def associate_fbgn_to_go_term(type_filter=None):
    fly = defaultdict(lambda: defaultdict(set))
    with open(association) as fh:
        for row in fh.readlines():
            if row.startswith("!"):
                continue

            cols = row.split("\t")
            fbgn = cols[1]
            goterm = cols[4]
            goaspect = {"P": "BP", "F": "MF", "C": "CC"}.get(cols[8], "")
            gtype = cols[11]

            if type_filter is not None:
                if gtype != type_filter:
                    continue

            if goterm not in oboDag:
                continue

            fly[goaspect][fbgn].add(goterm)

    return fly


def get_fly_go_full(type_filter=None):
    fly = associate_fbgn_to_go_term(type_filter=type_filter)
    go2fly = defaultdict(lambda: defaultdict(set))
    for aspect, genes in fly.items():
        for fbgn, goterms in genes.items():
            for goterm in goterms:
                go2fly["aspect"][fbgn].add(goterm)
    return go2fly


def get_fly_go_slim(type_filter=None):
    fly = associate_fbgn_to_go_term(type_filter=type_filter)
    flyslim = defaultdict(lambda: defaultdict(set))
    for aspect, genes in fly.items():
        for fbgn, goterms in genes.items():
            all_direct_anc = set()
            all_covered_anc = set()
            all_all_anc = set()
            for goterm in goterms:
                direct_anc, all_anc = mapslim(goterm, oboDag, slimDag)
                all_all_anc |= all_anc
                # collect all covered ancestors, so the direct ancestors
                # can be calculated afterwards
                all_covered_anc |= all_anc - direct_anc
            all_direct_anc = all_all_anc - all_covered_anc

            if len(all_direct_anc) > 0:
                flyslim[aspect][fbgn] |= all_direct_anc
    return flyslim


def print_example_command():
    print(
        dedent(
            """
        goObj = GOEnrichmentStudy(
            background,
            flyslim,
            slimDag,
            propagate_counts=False,
            alpha=0.001,
            methods=['fdr_bh']
        )
        """
        )
    )


def run_fly(genes, background, cutoff=0.001, return_obj=False):
    aspect2fbgn2go = get_fly_go_full()
    goObj = GOEnrichmentStudyNS(
        background,
        aspect2fbgn2go,
        oboDag,
        propagate_counts=False,
        alpha=cutoff,
        methods=["fdr_bh"],
    )

    if return_obj:
        return [res for res in goObj.run_study(genes) if res.p_fdr_bh <= cutoff], goObj
    else:
        return [res for res in goObj.run_study(genes) if res.p_fdr_bh <= cutoff]


def run_flyslim(genes, background, cutoff=0.001, return_obj=False):
    aspect2fbgn2go = get_fly_go_slim()
    goObj = GOEnrichmentStudyNS(
        background,
        aspect2fbgn2go,
        slimDag,
        propagate_counts=False,
        alpha=cutoff,
        methods=["fdr_bh"],
    )

    if return_obj:
        return [res for res in goObj.run_study(genes) if res.p_fdr_bh <= cutoff], goObj
    else:
        return [res for res in goObj.run_study(genes) if res.p_fdr_bh <= cutoff]


def go_wordcloud(results):
    from wordcloud import WordCloud

    freqs = {}
    for r in results:
        if r.name in ["biological_process", "cellular_component", "molecular_function"]:
            continue

        if r.study_count > 20:
            continue

        freqs[r.name] = r.study_count
        print(r.GO, r.name, r.study_count)

    wc = WordCloud(max_font_size=40).generate_from_frequencies(freqs)
    plt.imshow(wc, interpolation="bilinear")
    plt.axis("off")

    return plt.gca()
