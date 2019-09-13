"""Shared bits of code."""

from pathlib import Path
from pickle import load

import matplotlib.pyplot as plt
import pandas as pd

from larval_gonad.config import CONFIG_DIR
import larval_gonad.plotting


# Setup general plotting info
plt.style.use(['common', 'paper'])

# make mappers
try:
    with Path('../output/science_submission/fbgn2symbol.pkl').open('rb') as fh:
        fbgn2symbol = load(fh)
except FileNotFoundError:
    print('Please make sure to run the rules fbgn2symbol first.')
    raise

try:
    with Path('../output/science_submission/symbol2fbgn.pkl').open('rb') as fh:
        symbol2fbgn = load(fh)
except FileNotFoundError:
    print('Please make sure to run the rules symbol2fbgn first.')
    raise

try:
    fbgn2chrom = pd.read_parquet('../output/science_submission/fbgn2chrom.parquet')
except FileNotFoundError:
    print('Please make sure to run the rules fbgn2chrom first.')
    raise
