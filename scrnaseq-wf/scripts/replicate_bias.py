import pandas as pd

from larval_gonad.stats import run_chisq


file_name = snakemake.input.num_cells
output_name = snakemake.output[0]


def main():
    df = pd.read_parquet(file_name)
    ct = df.unstack().T
    ct.index = ct.index.droplevel(0)
    res = run_chisq(ct.T.fillna(0)).loc[(slice(None), ['adj std residual', 'flag_sig']), :]

    summary = {'flag_sig_bias': {}}
    for clus, dd in res.groupby('cluster'):
        dd.index = dd.index.droplevel(0)
        for rep, ddd in dd.T.iterrows():
            summary['flag_sig_bias'][(clus, rep)] = 0
            if ddd.flag_sig:
                if ddd['adj std residual'] < 0:
                    summary['flag_sig_bias'][(clus, rep)] = -1
                elif ddd['adj std residual'] > 0:
                    summary['flag_sig_bias'][(clus, rep)] = 1

    flag_bias = (
        pd.DataFrame(summary)
        .rename_axis(['cluster', 'rep'])
    )
    (
        df.join(flag_bias, on=['cluster', 'rep'])
        .to_parquet(output_name)
    )


if __name__ == '__main__':
    main()
