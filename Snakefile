import yaml
from pathlib import Path

import pandas as pd

from larval_gonad.config import read_config

# Build References
common_config = read_config('config/common.yaml')
assembly = common_config['assembly']
tag = common_config['tag']


rule targets:
    input:
        expand('references/{assembly}-all-{tag}.modFix.gtf', assembly=assembly, tag=tag),
        expand('references/fbgn2symbol-{assembly}-all-{tag}.feather', assembly=assembly, tag=tag),
        expand('references/fbgn2chrom-{assembly}-all-{tag}.feather', assembly=assembly, tag=tag),


rule fix_gtf:
    input: 'lcdb-references/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
    output: 'references/{assembly}-all-{tag}.modFix.gtf'
    shell: """
        awk 'BEGIN{{FS="\t"}}{{
            if (!($9 ~ /FBgn0002781/ && ($7 == "+" || $7 == "."))){{print $0}}
        }}' {input} > {output[0]}
    """


rule fbgn2symbol:
    input: rules.fix_gtf.output[0]
    output: 'references/fbgn2symbol-{assembly}-all-{tag}.tsv'
    script: 'scripts/fbgn2symbol.py'


rule fbgn2symbol_feather:
    input: rules.fbgn2symbol.output[0]
    output: 'references/fbgn2symbol-{assembly}-all-{tag}.feather'
    run:
        df = pd.read_csv(input[0], sep='\t').sort_values('FBgn').reset_index(drop=True)
        df.to_feather(output[0])


rule fbgn2chrom:
    input: rules.fix_gtf.output[0]
    output: 'references/fbgn2chrom-{assembly}-all-{tag}.tsv'
    script: 'scripts/fbgn2chrom.py'


rule fbgn2chrom_feather:
    input: rules.fbgn2chrom.output[0]
    output: 'references/fbgn2chrom-{assembly}-all-{tag}.feather'
    run:
        df = pd.read_csv(input[0], sep='\t').sort_values('FBgn').reset_index(drop=True)
        df.to_feather(output[0])