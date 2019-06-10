import re

import pandas as pd

barcodes = snakemake.input[0]
cellids = snakemake.output[0]
sample = snakemake.params.sample

def main():
    rep = 'rep' + re.findall('.*(\d+)', sample)[0]
    with open(barcodes) as fin, open(cellids, 'w') as fout:
        for row in fin.read().strip().split('\n'):
            cell = rep + '_' + re.sub('-1', '', row) + '\n'
            fout.write(cell)

if __name__ == '__main__':
    main()