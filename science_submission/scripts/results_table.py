"""Results Tables Aggregation

This script takes all results tables and puts them into a single XLSX workbook.
"""
import yaml
from itertools import zip_longest

import pandas as pd


config_name = snakemake.input['config']
xlsx_name = snakemake.output[0]


def main():
    config = get_config(config_name)

    # Set-up Workbook
    writer = pd.ExcelWriter(xlsx_name)
    workbook = writer.book

    cell_format = workbook.add_format({'valign': 'top', 'align': 'left'})
    cell_format.set_text_wrap()
    bold = workbook.add_format({'bold': True})

    readme = workbook.add_worksheet('README')

    readme.set_column(0, 0, width=30, cell_format=cell_format)
    readme.set_colum(1, 1, width=120, cell_format=cell_format)

    row = 0
    for block_name, block in config.items():
        # Write README information
        readme.write(row, 0, block_name)
        txt = zip_longest([block['description'].split('**')], [bold], fillvalue=bold)
        readme.write(row, 1, block['description'])
        row += 1

        # Write data as a worksheet
        if block.get('file_name', False):
            df = block['function'](block['file_name'])
            df.to_excel(writer, sheet_name=block_name)

    writer.close()


def get_config(file_name):
    with open(file_name) as fh:
        return yaml.load(fh.read())


if __name__ == '__main__':
    main()
