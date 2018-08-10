"""Shared bits of code."""

from pathlib import Path
from pickle import load

import matplotlib as mpl
from svgutils.compose import Text
import pandas as pd

from larval_gonad.config import CONFIG_DIR
from larval_gonad.plotting import add_styles


def label(text):
    return Text(text, 5, 10, size=12, weight="bold", font="sans-serif")


# Setup general plotting info
add_styles(Path(CONFIG_DIR, 'stylelib'))
mpl.style.use(['common', 'paper'])

# make mappers
try:
    with Path('data/fbgn2symbol.pkl').open('rb') as fh:
        fbgn2symbol = load(fh)
except FileNotFoundError:
    print('Please make sure to run the rules fbgn2symbol first.')
    raise

try:
    with Path('data/symbol2fbgn.pkl').open('rb') as fh:
        symbol2fbgn = load(fh)
except FileNotFoundError:
    print('Please make sure to run the rules symbol2fbgn first.')
    raise

try:
    fbgn2chrom = pd.read_parquet('data/fbgn2chrom.parquet')
except FileNotFoundError:
    print('Please make sure to run the rules fbgn2chrom first.')
    raise
