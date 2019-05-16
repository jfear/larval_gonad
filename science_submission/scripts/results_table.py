"""Results Tables Aggregation

This script takes all results tables and puts them into a single XLSX workbook.
"""
import pickle
import yaml

import pandas as pd


results_config = yaml.load(open(snakemake.input.config_name))
fbgn2symbol = pickle.load(open(snakemake.input.fbgn2symbol, 'rb'))
fbgn2chrom = pd.read_parquet(snakemake.input.fbgn2chrom)

cluster_annot = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
legend_names = snakemake.params.legend_names

xlsx_name = snakemake.output[0]


def main():
    # Set-up Workbook
    writer = pd.ExcelWriter(xlsx_name)
    workbook = writer.book

    cell_format = workbook.add_format({'valign': 'top', 'align': 'left', 'font_size': 12})
    cell_format.set_text_wrap()
    header = workbook.add_format({'valign': 'vcenter', 'bold': True, 'font_size': 18})
    bold = workbook.add_format({'bold': True, 'font_size': 12})

    readme = workbook.add_worksheet('README')

    readme.set_column(0, 0, width=50, cell_format=cell_format)
    readme.set_column(1, 1, width=120, cell_format=cell_format)

    # Write out Cell Type IDs
    readme.write(0, 0, 'Cell Type', header)
    readme.write(0, 1, '\n'.join(legend_names))

    # Write out other table descriptions from results_table_config.yaml
    row = 1
    for block_name, block in results_config.items():
        # Write README information
        readme.write(row, 0, block_name, header)
        readme.write_rich_string(row, 1, *format_block(block['description'], cell_format, bold))
        readme.set_row(row, len(block['description'].split('\n')) * 12)
        row += 1

        # Write data as a worksheet
        if block.get('file_name', False):
            df = globals()[block['function']](block['file_name'])
            df.to_excel(writer, sheet_name=block_name)

            # Format indexed columns so they look nicer
            _sheet = workbook.get_worksheet_by_name(block_name)
            for i, idx in enumerate(df.index.names):
                max_len = df.index.get_level_values(idx).str.len().max()
                _sheet.set_column(i, i, width=max_len + 1)

    writer.close()


def format_block(block_text, default, bold):
    txts = block_text.split('**')

    if len(txts) == 1:
        return txts

    fmt = []
    for i, txt in enumerate(txts):
        if i % 2 == 1:
            fmt.append(bold)
            fmt.append(txt)
        else:
            fmt.append(default)
            fmt.append(txt)
    return fmt


def counts_table(file_name):
    w_rep = (
        pd.read_parquet(file_name)
        .assign(gene_symbol=lambda df: df.index.map(fbgn2symbol))
        .assign(chrom=lambda df: df.index.map(fbgn2chrom.iloc[:, 0].to_dict()))
        .set_index(['gene_symbol', 'chrom'], append=True)
        .assign(cluster=lambda df: pd.Categorical(df.cluster.map(cluster_annot), ordered=True, categories=cluster_order))
        .assign(rep=lambda df: pd.Categorical(df.rep, ordered=True, categories=['rep1', 'rep2', 'rep3']))
        .sort_values(by=['cluster', 'rep'])
        .set_index(['rep', 'cluster'], append=True)
        .unstack().unstack()
    )

    w_rep.columns = w_rep.columns.droplevel(0)
    return w_rep


def biomarkers(file_name):
    return (
        pd.read_parquet(file_name)
        .assign(gene_symbol=lambda df: df.index.map(fbgn2symbol))
        .assign(chrom=lambda df: df.index.map(fbgn2chrom.iloc[:, 0].to_dict()))
        .assign(cluster=lambda df: pd.Categorical(df.cluster.map(cluster_annot), ordered=True, categories=cluster_order))
        .sort_values(by='cluster')
        .reset_index()
        .melt(id_vars=['FBgn', 'gene_symbol', 'chrom', 'cluster'], var_name='seurat_output', value_name='value')
        .set_index(['FBgn', 'gene_symbol', 'chrom', 'seurat_output', 'cluster'])
        .unstack().unstack()
        .sort_index(level=0)
    )


def tsne(file_name):
    return pd.read_parquet(file_name)


def deg_res(file_name):
    return (
        pd.read_csv(file_name, sep='\t', index_col=0)
            .rename_axis('FBgn')
            .join(fbgn2chrom)
            .set_index(['gene_symbol', 'chrom'], append=True)
            .sort_index(level=0)
    )


def autosome_ratios(file_name):
    return (
        pd.read_parquet(file_name[0])
        .rename_axis('cell_id')
        .reindex(['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrY'], axis=1)
        .join(
            pd.read_parquet(file_name[1])
            .rename_axis('cell_id')
        )
        .assign(cluster=lambda df: pd.Categorical(df.cluster.map(cluster_annot), ordered=True, categories=cluster_order))
        .assign(rep=lambda df: pd.Categorical(df.rep, ordered=True, categories=['rep1', 'rep2', 'rep3']))
        .set_index(['cluster', 'rep'], append=True)
    )


def autosome_permutation(file_name):
    return (
        pd.read_parquet(file_name)
        .reset_index()
        .assign(cluster=lambda df: pd.Categorical(df.cluster.map(cluster_annot), ordered=True, categories=cluster_order))
        .set_index('cluster')
        .sort_index()
    )


if __name__ == '__main__':
    main()
