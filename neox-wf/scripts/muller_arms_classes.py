"""Gene movement between D. mel and D. pse

Use orhtologous genes between D. mel and D. pse to figure out gene movement. Maps each gene
to their corresponding Muller element and then determines if the genes have moved or become lost
during evolution focusing on the Neo-X chromosome (D) in D. pse.

Here we classify genes as:
* Conserved
* MovedOn
* MovedOff
* Lost

"""
import numpy as np
import pandas as pd

INPUT_FILE = snakemake.input[0]
OUTPUT_A = snakemake.output["muller_a"]
OUTPUT_D = snakemake.output["muller_d"]
OUTPUT_E = snakemake.output["muller_e"]

# Debug settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'neox-wf/scripts'))
#     print(os.getcwd())
# except:
#     pass
# INPUT_FILE = "../../output/neox-wf/muller_arm_assignment.feather"


def main():
    muller = pd.read_feather(INPUT_FILE).set_index("FBgn")
    muller_classes = pd.Series(index=muller.index)
    other_species = ["dmel", "dvir", "dmoj", "dwil"]

    # Muller arm is the same for all species
    _logic = muller.apply(lambda x: len(set(x)) == 1, axis=1)
    muller_classes[_logic] = "conserved"

    # Muller arm is recently conserved in D. mel, D. pse, and D. wil
    _logic = (
        # conserved in these three species
        muller.apply(lambda x: len(set(x[["dmel", "dpse", "dwil"]])) == 1, axis=1)
        # not conserved in all species
        & (muller_classes != "conserved")
    )
    muller_classes[_logic] = "recent_conserved"

    # Gene Death: gene is present in all species but missing in D. pse
    _logic = muller[other_species].notna().any(axis=1) & muller.dpse.isna()
    muller_classes[_logic] = "gene_death"

    # Moved: The gene is on a different arm in D. pse vs other species
    _logic = (
        # same in other species
        muller.apply(lambda x: len(set(x[other_species])) == 1, axis=1)
        # different in D. dpse
        & (muller.dmel != muller.dpse)
        # Ignore if the gene is there but on an unknown scaffold
        & (muller.dpse != "unknown_scaffold")
    )
    muller_classes[_logic] = "moved"

    # Other types of gene classes
    muller_classes.fillna("other", inplace=True)

    # Zoom in on specific arms and classify
    muller_a = classify_specific_arm(muller, muller_classes, "A")
    muller_a.reset_index().to_feather(OUTPUT_A)

    muller_d = classify_specific_arm(muller, muller_classes, "D")
    muller_d.reset_index().to_feather(OUTPUT_D)

    muller_e = classify_specific_arm(muller, muller_classes, "E")
    muller_e.reset_index().to_feather(OUTPUT_E)


def classify_specific_arm(muller, muller_classes, arm="A"):
    df = muller_classes.copy()
    df.name = f"muller_{arm}"
    # If gene was on in other species but moved off in D. pse
    df.loc[(muller.dmel == arm) & (df == "moved")] = "moved_off"
    # If gene was on another arm in other species but moved on in D. pse
    df.loc[(muller.dpse == arm) & (df == "moved")] = "moved_on"
    df = df.reindex(muller.query(f'dpse == "{arm}" | dmel == "{arm}"').index)
    return df.to_frame()


if __name__ == "__main__":
    main()