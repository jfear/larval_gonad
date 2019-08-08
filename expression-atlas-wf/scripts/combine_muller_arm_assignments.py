import pandas as pd

YO_FILE = snakemake.input["yo"]
FB_FILE = snakemake.input["fb"]
OUTPUT_FILE = snakemake.output[0]

# Debug settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), "expression-atlas-wf/scripts"))
#     print(os.getcwd())
# except:
#     pass
# YO_FILE = "../../output/expression-atlas-wf/YO_muller_arm_assignment.feather"
# FB_FILE = "../../output/expression-atlas-wf/FB_muller_arm_assignment.feather"


def main():
    yo = pd.read_feather(YO_FILE).set_index("FBgn")
    fb = pd.read_feather(FB_FILE).set_index("FBgn")
    combined = yo.combine_first(fb)

    combined.reset_index().to_feather(OUTPUT_FILE)


if __name__ == "__main__":
    main()
