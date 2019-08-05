"""Protein trap image scores.

Miriam was incharge of running the protein traps and summarizing several scores
results. She emailed me a excel sheet whose contents are just pated in ths
script. I added the FBgns by hand.

"""

from io import StringIO

import pandas as pd

oname = snakemake.output[0]

# This data is from a spreadsheet that Miriam gave me. I had to add on FBgns
dat = """FBgn,SP,EPS,MPS,LPS,C,ECY,MCY,LCY,PC,TE
FBgn0026573,1,2,3,3,1,0,0,0,0,1
FBgn0012037,1,2,2,2,1,2,2,2,1,1
FBgn0011206,1,1,2,3,0,0,0,0,0,0
FBgn0264494,0,0,0,0,1,1,1,1,1,1
FBgn0039754,0,1,1,1,0,1,1,1,1,1
FBgn0038180,0,0,0,0,0,1,1,2,0,0
FBgn0027598,1,1,1,1,1,1,1,1,2,2
FBgn0037015,0,0,0,0,0,1,1,1,0,0
FBgn0026533,2,2,2,2,1,1,1,1,2,2
FBgn0051361,1,1,1,0,0,0,0,0,0,0
FBgn0087008,2,2,2,2,1,1,1,1,2,2
FBgn0085420,2,2,1,1,0,0,0,0,1,2
FBgn0051158,1,1,1,1,1,2,2,2,1,1
FBgn0051158,2,1,1,1,2,2,2,2,1,1
FBgn0000636,1,1,1,1,1,1,1,1,1,3
FBgn0000636,0,0,0,0,0,0,0,0,0,2
FBgn0005633,1,1,1,1,2,2,2,2,1,1
FBgn0038197,2,2,2,2,1,1,1,1,1,2
FBgn0262743,2,2,2,2,1,1,1,1,1,1
FBgn0034282,2,2,2,1,1,1,1,1,2,2
FBgn0265487,1,2,2,2,1,0,0,0,1,3
FBgn0083963,1,1,1,1,1,2,2,1,0,1
FBgn0083963,0,0,0,0,0,1,1,1,0,0
FBgn0050418,1,1,1,1,0,0,0,0,0,2
FBgn0264975,2,1,1,1,2,2,2,2,2,2
FBgn0261885,1,1,1,1,2,2,2,2,2,2
FBgn0039044,2,2,1,0,0,0,0,0,0,1
FBgn0030520,3,3,3,3,1,2,2,2,3,3
FBgn0264953,0,0,0,0,0,0,0,0,0,3
FBgn0243486,1,1,1,1,2,2,2,2,1,1
FBgn0243486,1,1,1,1,2,2,3,3,1,1
FBgn0243486,1,1,1,1,2,2,3,3,1,1
FBgn0083981,1,1,1,1,1,1,1,1,1,1
FBgn0000416,1,1,1,1,2,3,3,3,2,3
FBgn0039232,1,1,1,1,1,2,2,2,1,1
FBgn0003475,1,1,1,1,0,0,0,0,0,0
FBgn0026370,2,2,2,2,0,0,0,0,1,1
FBgn0266521,3,2,2,2,2,2,2,2,1,3
FBgn0004575,1,1,1,1,1,1,1,2,1,1
FBgn0004575,1,1,1,1,1,1,1,1,1,1
FBgn0041182,0,1,1,1,0,2,2,2,2,0
FBgn0004885,1,1,1,1,1,0,0,0,1,2
FBgn0010473,1,1,1,1,0,0,0,0,1,1
FBgn0011725,1,1,1,1,0,0,0,0,0,0"""

df = pd.read_csv(StringIO(dat), sep=',')
df.to_parquet(oname)
