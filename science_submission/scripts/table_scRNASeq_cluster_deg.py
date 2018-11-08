from textwrap import dedent

import pandas as pd

from larval_gonad.config import config

from common import fbgn2symbol

def _add_sheet(writer : pd.ExcelWriter,
               sheet_name : str,
               fname : str,
               cell_format,
               alpha=0.01,
               comment=None,
               headers=None):

    sheet = writer.book.add_worksheet(sheet_name)
    writer.sheets[sheet_name] = sheet
    sheet.set_column(0, 1, 20)

    df = pd.read_csv(fname, sep='\t').query(f'p_val_adj <= {alpha}')
    df.sort_values(by='avg_logFC', ascending=False, inplace=True)

    if 'cluster' in df.columns:
        df.sort_values(by='cluster', inplace=True)
        df.cluster.replace(config['cluster_annot'], inplace=True)
        idx = df.columns.tolist().index('cluster')
        sheet.set_column(idx, idx, 20)

    if headers:
        df.rename({'pct.1': headers[0], 'pct.2': headers[1]}, inplace=True, axis=1)
    df.to_excel(writer, sheet_name=sheet_name, index=False, startrow=1, freeze_panes=(2, 2))

    if comment:
        sheet.set_row(0, 100, cell_format)
        sheet.merge_range('A1:G1', dedent(comment))


def _add_data(writer : pd.ExcelWriter,
               sheet_name : str,
               fname : str,
               cell_format,
               comment=None):

    sheet = writer.book.add_worksheet(sheet_name)
    writer.sheets[sheet_name] = sheet
    sheet.set_column(0, 1, 20)

    df = pd.read_parquet(fname)[config['sel_cluster_order_w_rep']]
    df['gene'] = df.index.map(lambda x: fbgn2symbol[x])
    df.set_index('gene', append=True, inplace=True)
    df.to_excel(writer, sheet_name=sheet_name, index=True, startrow=1, freeze_panes=(2, 2))

    if comment:
        sheet.set_row(0, 100, cell_format)
        sheet.merge_range('A1:G1', dedent(comment))


def main():
    writer = pd.ExcelWriter(snakemake.output[0])
    cell_format = writer.book.add_format({'valign': 'top'})
    cell_format.set_text_wrap()

    comment = "Identification of Bio markers. " \
              "Here we compare expression of each cluster to all other cells. " \
              "This creates a list of genes that are up-regulated in each cluster. " \
              "This table is grouped by clusters and sorted by avg_logFC."

    _add_sheet(writer, 'One vs Rest (biomarkers)', snakemake.input.biomarkers, cell_format, comment=comment)

    comment = "Differential expression between germ cell and somatic cell clusters. " \
              "I combine all of the germ cell clusters (6, 3, 2, 0) vs all of the somatic cell "\
              "clusters (5, 1, 4, 7, 8).\n" \
              "Positive avg_logFC are germ biased genes.\n" \
              "Negative avg_logFC are soma biased genes.\n"

    _add_sheet(writer, 'Germ Cells vs Somatic Cells', snakemake.input.germ_soma, cell_format, comment=comment,
               headers=('pct.germ', 'pct.soma'))

    comment = "Differential expression of spermatogonia vs 1º spermatocytes. " \
              "Spermatogonia cluster and compared it to all spermatocyte clusters combined together.\n\n" \
              "Positve avg_logFC are spermatogonia biased genes.\n" \
              "Negative avg_logFC are 1º spermatocyte biased genes."

    _add_sheet(writer, 'Gonia vs Cytes', snakemake.input.gonia_cytes, cell_format, comment=comment,
               headers=('pct.gonia', 'pct.cytes'))

    comment = "Differential expression of Early 1º spermatocytes vs Mid and Late 1º spermatocytes.\n" \
              "Positve avg_logFC are early 1º spermatocyte biased genes.\n" \
              "Negative avg_logFC are mid and late 1º spermatocyte biased genes."

    _add_sheet(writer, 'Early cytes vs Mid and Late', snakemake.input.gonia_early, cell_format, comment=comment,
               headers=('pct.early', 'pct.midLate'))

    comment = "Differential expression of Mid 1º spermatocytes vs Early and Late 1º spermatocytes.\n" \
              "Positve avg_logFC are mid 1º spermatocyte biased genes.\n" \
              "Negative avg_logFC are early and late 1º spermatocyte biased genes."

    _add_sheet(writer, 'Mid cytes vs Early and Late', snakemake.input.gonia_mid, cell_format, comment=comment,
               headers=('pct.mid', 'pct.earlyLate'))

    comment = "Differential expression of Late 1º spermatocytes vs Early and Mid 1º spermatocytes.\n" \
              "Positve avg_logFC are late 1º spermatocyte biased genes.\n" \
              "Negative avg_logFC are early and mid 1º spermatocyte biased genes."

    _add_sheet(writer, 'Late cytes vs Early and Mid', snakemake.input.gonia_late, cell_format, comment=comment,
               headers=('pct.late', 'pct.earlyMid'))

    comment = "These are raw counts aggregated by Cluster:Rep using Sum."

    _add_data(writer, 'Raw Counts (Sum)', snakemake.input.raw, cell_format, comment=comment)

    comment = "These are tpm normalized counts by Cluster:Rep."

    _add_data(writer, 'TPM', snakemake.input.tpm, cell_format, comment=comment)

    comment = "These are z-scores of tpm normalized counts by Cluster:Rep."

    _add_data(writer, 'Z-scores', snakemake.input.zscore, cell_format, comment=comment)

    writer.close()


if __name__ == '__main__':
    main()