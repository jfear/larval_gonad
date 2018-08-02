"""Create dataframe with cluster information.

During the analysis I generate slightly different cluster sets depending on the
clustering settings and method. Here I pull out the clusters based on the
provided resolution. I then add on the cluster annotations for easier reading.

"""

import pandas as pd

from larval_gonad.config import config

fname = snakemake.input[0]
res = snakemake.params.resolution
oname = snakemake.output[0]

# import cluster calls for specific resolution
sr = pd.read_csv(fname, sep='\t')[res]
sr.name = 'cluster_num'
df = sr.to_frame()

# Add on cluster annotations for easier understanding
df['cluster_name'] = df['cluster_num'].map(config['cluster_annot'])
df.to_parquet(oname)
