import pandas as pd

from larval_gonad.config import read_config

# Build References
common_config = read_config('config/common.yaml')
assembly = common_config['assembly']
tag = common_config['tag']


rule targets:
    input:
        expand('references/{assembly}-all-{tag}.modFix.gtf', assembly=assembly, tag=tag),
        expand('references/gene_annotation_{assembly}_{tag}.feather', assembly=assembly, tag=tag)


rule fix_gtf:
    input: 'lcdb-references/{assembly}/{tag}/gtf/{assembly}_{tag}.gtf'
    output: 'references/{assembly}-all-{tag}.modFix.gtf'
    shell: """
        awk 'BEGIN{{FS="\t"}}{{
            if (!($9 ~ /FBgn0002781/ && ($7 == "+" || $7 == "."))){{print $0}}
        }}' {input} > {output[0]}
    """


rule gene_annotation:
    input: rules.fix_gtf.output[0]
    output: 'references/gene_annotation_{assembly}_{tag}.tsv'
    script: 'scripts/gene_annotation.py'


rule gene_annotation_feather:
    input: rules.gene_annotation.output[0]
    output: 'references/gene_annotation_{assembly}_{tag}.feather'
    run:
        df = pd.read_csv(input[0], sep='\t', keep_default_na=False, na_values=['99999999999']).sort_values('FBgn').reset_index(drop=True)
        df.to_feather(output[0])
