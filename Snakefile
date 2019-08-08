import pandas as pd

from larval_gonad.config import read_config

# Build References
common_config = read_config('config/common.yaml')
assembly = common_config['assembly']
tag = common_config['tag']
release_date = common_config['release_date']


rule targets:
    input:
        f'references/{assembly}-all-{tag}.modFix.gtf',
        f'references/gene_annotation_{assembly}_{tag}.feather',
        f'references/primary2secondary_{assembly}_{tag}.pkl',
        'expression-atlas-wf/config/sampletable.tsv'


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


rule primary_secondary_fbgn_mapper:
    params: f"ftp://ftp.flybase.net/releases/FB{release_date}/precomputed_files/genes/fbgn_annotation_ID_fb_{release_date}.tsv.gz"
    output: f'references/primary2secondary_{assembly}_{tag}.pkl'
    script: 'scripts/primary_secondary_fbgn_mapper.py'


rule generate_haiwang_sampletable:
    output: 'expression-atlas-wf/config/sampletable.tsv'
    params: 'GSE99574'
    script: 'scripts/parse_geo_to_sampletable.py'