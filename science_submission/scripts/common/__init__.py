"""Shared bits of code."""

from pathlib import Path

import matplotlib as mpl
import seaborn as sns
from svgutils.compose import Text

from larval_gonad.config import config, CONFIG_DIR
from larval_gonad.plotting import add_styles


def label(text):
    return Text(text, 5, 15, size=12, weight="bold", font="sans-serif")


# Setup general plotting info
add_styles(Path(CONFIG_DIR, 'stylelib'))
mpl.style.use(['common', 'paper'])
#sns.set(style='ticks', context='paper', palette="tab10")
#sns.set_palette('colorblind')
