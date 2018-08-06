from textwrap import dedent
from collections import defaultdict

import matplotlib.pyplot as plt

from goatools import GOEnrichmentStudy
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj
from wordcloud import WordCloud

oboDag = GODag('../data/external/go-basic.obo')
slimDag = GODag("../data/external/goslim_generic.obo")
association = '../data/external/gene_association.fb'


def get_fly_assoc(type_filter=None):
    fly = defaultdict(set)
    with open(association) as fh:
        for row in fh.readlines():
            if row.startswith('!'):
                continue

            cols = row.split('\t')
            fbgn = cols[1]
            goterm = cols[4]
            gtype = cols[11]

            if type_filter is not None:
                if gtype != type_filter:
                    continue

            if goterm not in oboDag:
                continue

            fly[fbgn].add(goterm)

    return fly


def get_go2fly(type_filter=None):
    fly = get_fly_assoc(type_filter=type_filter)
    go2fly = defaultdict(set)
    for k, v in fly.items():
        for t in v:
            if oboDag[t].namespace != 'biological_process':
                continue

            go2fly[t].add(k)

    return go2fly


def get_flyslim(type_filter=None):
    fly = get_fly_assoc(type_filter=type_filter)
    flyslim = defaultdict(set)
    for fbgn, goterms in fly.items():
        all_direct_anc = set()
        all_covered_anc = set()
        all_all_anc = set()
        for goterm in goterms:
            direct_anc, all_anc = mapslim(goterm, oboDag, slimDag)
            all_all_anc |= all_anc
            # collect all covered ancestors, so the direct ancestors
            # can be calculated afterwards
            all_covered_anc |= (all_anc - direct_anc)
        all_direct_anc = all_all_anc - all_covered_anc
        flyslim[fbgn] |= (all_direct_anc)

    return flyslim


def print_example_command():
    print(dedent("""
        goObj = GOEnrichmentStudy(
            background,
            flyslim,
            slimDag,
            propagate_counts=False,
            alpha=0.001,
            methods=['fdr_bh']
        )
        """))


def run_fly(genes, background, cutoff=0.001, return_obj=False):
    fly = get_fly_assoc()
    goObj = GOEnrichmentStudy(
        background,
        fly,
        oboDag,
        propagate_counts=False,
        alpha=cutoff,
        methods=['fdr_bh']
    )

    if return_obj:
        return [res for res in goObj.run_study(genes) if res.p_fdr_bh <= cutoff], goObj
    else:
        return [res for res in goObj.run_study(genes) if res.p_fdr_bh <= cutoff]


def run_flyslim(genes, background, cutoff=0.001, return_obj=False):
    flyslim = get_flyslim()
    goObj = GOEnrichmentStudy(
        background,
        flyslim,
        slimDag,
        propagate_counts=False,
        alpha=cutoff,
        methods=['fdr_bh']
    )

    if return_obj:
        return [res for res in goObj.run_study(genes) if res.p_fdr_bh <= cutoff], goObj
    else:
        return [res for res in goObj.run_study(genes) if res.p_fdr_bh <= cutoff]


def go_wordcloud(results):
    freqs = {}
    for r in results:
        if r.name in ['biological_process', 'cellular_component', 'molecular_function']:
            continue

        if r.study_count > 20:
            continue

        freqs[r.name] = r.study_count
        print(r.GO, r.name, r.study_count)

    wc = WordCloud(max_font_size=40).generate_from_frequencies(freqs)
    plt.imshow(wc, interpolation='bilinear')
    plt.axis("off")

    return plt.gca()
