import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

builtins.snakemake = MockSnake(
    input="../../output/paper_submission/fig1_data_avg_tpm_per_chrom.feather",
    params=dict(colors=read_config("../../config/colors.yaml")),
)
