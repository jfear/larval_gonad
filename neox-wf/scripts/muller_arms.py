"""Gene movement between D. mel and D. pse

Use orhtologous genes between D. mel and D. pse to figure out gene movement. Maps each gene 
to their corresponding Muller element and then determines if the genes have moved or become lost 
during evolution focusing on the Neo-X chromosome (D) in D. pse. 

Here we classify genes as:
* Conserved
* MovedOn
* MovedOff
* Lost

We then group these classed into:
* Selection for (SF)
* Selection against (SA)

"""
import numpy as np
import pandas as pd

INPUT_FILE = snakemake.input[0]
OUTPUT_FILE = snakemake.output[0]

# Debug settings
# INPUT_FILE = "../../output/neox-wf/ortholog_annotation.feather"

# Muller arm mappings from Maria
MAPPERS = {
    "dmel": {
        "X": "A",
        "2L": "B",
        "2R": "C",
        "3L": "D",
        "3R": "E",
        "4": "F",
        "Unmapped": np.nan,
        "211000022280328": np.nan,
    },
    "dpse": {"XL": "A", "4": "B", "3": "C", "XR": "D", "2": "E", "5": "F", "Unknown": np.nan},
}


def main():
    orthologs = (
        pd.read_feather(INPUT_FILE)
        .set_index("FBgn")
        .reindex(["dmel", "dpse", "dmoj", "dvir"], axis=1)
    )

    muller = orthologs.apply(extract_and_map, axis=0).dropna(how="all", axis=1)

    # Classify genes as conserved, moved off, moved on, or lost from dmel to dpse
    muller["classA"] = np.nan
    muller.loc[(muller.dmel == "A") & (muller.dpse == "A"), "classA"] = "ConservedA"
    muller.loc[
        ~((muller.dmel == "A") | muller.dmel.isna()) & (muller.dpse == "A"), "classA"
    ] = "MovedOnA"
    muller.loc[
        (muller.dmel == "A") & ~((muller.dpse == "A") | muller.dpse.isna()), "classA"
    ] = "MovedOffA"
    muller.loc[(muller.dmel == "A") & muller.dpse.isna(), "classA"] = "LostA"

    muller["classD"] = np.nan
    muller.loc[(muller.dmel == "D") & (muller.dpse == "D"), "classD"] = "ConservedD"
    muller.loc[
        ~((muller.dmel == "D") | muller.dmel.isna()) & (muller.dpse == "D"), "classD"
    ] = "MovedOnD"
    muller.loc[
        (muller.dmel == "D") & ~((muller.dpse == "D") | muller.dpse.isna()), "classD"
    ] = "MovedOffD"
    muller.loc[(muller.dmel == "D") & muller.dpse.isna(), "classD"] = "LostD"

    muller["classE"] = np.nan
    muller.loc[(muller.dmel == "E") & (muller.dpse == "E"), "classE"] = "ConservedE"
    muller.loc[
        ~((muller.dmel == "E") | muller.dmel.isna()) & (muller.dpse == "E"), "classE"
    ] = "MovedOnE"
    muller.loc[
        (muller.dmel == "E") & ~((muller.dpse == "E") | muller.dpse.isna()), "classE"
    ] = "MovedOffE"
    muller.loc[(muller.dmel == "E") & muller.dpse.isna(), "classE"] = "LostE"

    # Group classes as selection for A|D|E (SFD), selection against A|D|E (SAD)
    muller["groupA"] = np.nan
    muller.loc[
        muller.classA.str.contains("Conserved") | muller.classA.str.contains("MovedOn"), "groupA"
    ] = "SFA"
    muller.loc[muller.classA.str.contains("Off").fillna(False), "groupA"] = "SAA"
    muller.loc[muller.classA.str.contains("Lost").fillna(False), "groupA"] = "SAA?"

    muller["groupD"] = np.nan
    muller.loc[
        muller.classD.str.contains("Conserved") | muller.classD.str.contains("MovedOn"), "groupD"
    ] = "SFD"
    muller.loc[muller.classD.str.contains("Off").fillna(False), "groupD"] = "SAD"
    muller.loc[muller.classD.str.contains("Lost").fillna(False), "groupD"] = "SAD?"

    muller["groupE"] = np.nan
    muller.loc[
        muller.classE.str.contains("Conserved") | muller.classE.str.contains("MovedOn"), "groupE"
    ] = "SFE"
    muller.loc[muller.classE.str.contains("Off").fillna(False), "groupE"] = "SAE"
    muller.loc[muller.classE.str.contains("Lost").fillna(False), "groupE"] = "SAE?"

    # Create output datasets
    (
        orthologs.merge(muller, right_index=True, left_index=True, suffixes=["", "_muller"])
        .reset_index()
        .to_feather(OUTPUT_FILE)
    )


def extract_and_map(x):
    """Extract out chromosome name and map muller arm.
    
    Some of the chromosome names had suffixes, this function pulls our the
    chromosome name and then maps it to the correct muller arm.

    """
    _name = x.name
    mapper = MAPPERS.get(_name, {"scaffold": np.nan})
    return x.str.extract("([0-9a-zA-Z]+)", expand=False).replace(mapper)


if __name__ == "__main__":
    main()
