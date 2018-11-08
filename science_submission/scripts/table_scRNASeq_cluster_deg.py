import pandas as pd


def _add_sheet(writer : pd.ExcelWriter, sheet_name : str, fname : str, alpha=0.01, comment=None, headers=None):
    sheet = writer.book.add_worksheet(sheet_name)
    writer.sheets[sheet_name] = sheet
    sheet.set_column(0, 1, 20)

    df = pd.read_csv(Path(seurat_dir, fname), sep='\t').query(f'p_val_adj <= {alpha}')
    df.sort_values(by='avg_logFC', ascending=False, inplace=True)

    if 'cluster' in df.columns:
        df.sort_values(by='cluster', inplace=True)
        df.cluster.replace(CLUSTER_ANNOT, inplace=True)
        idx = df.columns.tolist().index('cluster')
        sheet.set_column(idx, idx, 20)

    if headers:
        df.rename({'pct.1': headers[0], 'pct.2': headers[1]}, inplace=True, axis=1)
    df.to_excel(writer, sheet_name=sheet_name, index=False, startrow=1, freeze_panes=(2, 2))

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

    _add_sheet(writer, 'One vs Rest (biomarkers)', f'biomarkers_{resolution}.tsv', comment=comment)

    comment = "Differential expression between germ cell and somatic cell clusters. " \
              "For this analysis I combine all of the germ cell clusters (6, 3, 2, 0) vs all of the somatic cell "\
              "clusters (5, 1, 4, 7, 8).\n" \
              "Positive avg_logFC are germ biased genes.\n" \
              "Negative avg_logFC are soma biased genes.\n"

    _add_sheet(writer, 'Germ Cells vs Somatic Cells', '2018-05-16_scrnaseq_germ_vs_soma_biomarkers.tsv',
               comment=comment, headers=('pct.germ', 'pct.soma'))

    comment = "Here I have done a differential expression of spermatogonia vs 1º spermatocytes. " \
              "For this analysis I took the spermatogonia cluster and compared it to all spermatocyte " \
              "clusters combined together.\n\n" \
              "Positve avg_logFC are spermatogonia biased genes.\n" \
              "Negative avg_logFC are 1º spermatocyte biased genes."

    _add_sheet(
        writer,
        'Gonia vs Cytes',
        '2018-05-16_scrnaseq_spermatogonia_vs_spermatocytes_biomarkers.tsv',
        comment=comment,
        headers=('pct.gonia', 'pct.cytes')
    )

    comment = "Here I have done a differential expression of Early 1º spermatocytes vs Mid and Late 1º spermatocytes.\n" \
              "Positve avg_logFC are early 1º spermatocyte biased genes.\n" \
              "Negative avg_logFC are mid and late 1º spermatocyte biased genes."

    _add_sheet(
        writer,
        'Early cytes vs Mid and Late',
        '2018-05-16_scrnaseq_early_spermatocytes_vs_spermatocytes_biomarkers.tsv',
        comment=comment,
        headers=('pct.early', 'pct.midLate')
    )

    comment = "Here I have done a differential expression of Mid 1º spermatocytes vs Early and Late 1º spermatocytes.\n" \
              "Positve avg_logFC are mid 1º spermatocyte biased genes.\n" \
              "Negative avg_logFC are early and late 1º spermatocyte biased genes."

    _add_sheet(
        writer,
        'Mid cytes vs Early and Late',
        '2018-05-16_scrnaseq_mid_spermatocytes_vs_spermatocytes_biomarkers.tsv',
        comment=comment,
        headers=('pct.mid', 'pct.earlyLate')
    )

    comment = "Here I have done a differential expression of Late 1º spermatocytes vs Early and Mid 1º spermatocytes.\n" \
              "Positve avg_logFC are late 1º spermatocyte biased genes.\n" \
              "Negative avg_logFC are early and mid 1º spermatocyte biased genes."

    _add_sheet(
        writer,
        'Late cytes vs Early and Mid',
        '2018-05-16_scrnaseq_late_spermatocytes_vs_spermatocytes_biomarkers.tsv',
        comment=comment,
        headers=('pct.late', 'pct.earlyMid')
    )

    comment = "Cluster eleven is a little bit of a mystery to us. " \
              "It behaves kind of like a 1º spermatocyte, but has very low expression. " \
              "Here I run a differential expression between cluster eleven and the 1º spermatocyte clusters.\n" \
              "Positve avg_logFC are cluster 11 biased genes.\n" \
              "Negative avg_logFC are 1º spermatocyte biased genes."

    _add_sheet(
        writer,
        'Cluster 11 vs cytes',
        '2018-05-16_scrnaseq_eleven_vs_spermatocytes_biomarkers.tsv',
        comment=comment,
        headers=('pct.eleven', 'pct.cytes')
    )

    writer.close()

if __name__ == '__main__':
    main()