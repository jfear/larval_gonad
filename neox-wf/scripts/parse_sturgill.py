import pandas as pd

XLS = snakemake.input[0]
MULLERD = snakemake.output.mullerD
MULLERE = snakemake.output.mullerE

# Debug Settings
# XLS = 'output/neox-wf/sturgill_2007.xls'

def main():
    pd.read_excel(
            XLS,
            sheet_name="TableS3-muller D",
            skiprows=3,
            nrows=243,
            usecols=["D.mel", "Gene fate"],
            na_values="."
    ).dropna().reset_index(drop=True).to_feather(MULLERD)


    pd.read_excel(
            XLS,
            sheet_name="TableS4-muller E",
            skiprows=3,
            nrows=427,
            usecols=["D.mel", "Fate"],
            na_values="."
    ).dropna().reset_index(drop=True).to_feather(MULLERE)


if __name__ == "__main__":
    main()
