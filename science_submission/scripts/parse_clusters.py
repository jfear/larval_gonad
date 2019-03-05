"""Create dataframe with cluster information.

During the analysis I generate slightly different cluster sets depending on the
clustering settings and method. Here I pull out the clusters based on the
provided resolution. I then add on the cluster annotations for easier reading.

"""
import pandas as pd

fname = snakemake.input[0]
resolution = snakemake.params.resolution
annotation = snakemake.params.annotation
oname = snakemake.output[0]

df = (
    pd.read_csv(fname, sep='\t')[resolution]
    .rename('cluster')
    .map(lambda x: annotation[x])
    .to_frame()
)

df.index.name = 'cell_id'
df.to_parquet(oname)

