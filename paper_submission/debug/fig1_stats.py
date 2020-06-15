import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input="../../output/paper_submission/fig1_data_avg_tpm_per_chrom.feather"
)
