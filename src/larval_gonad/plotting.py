from functools import wraps

from numpy import arange
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from .config import config

# colormaps
cluster_cmap = dict(zip(config["cluster_order"], config["colors"]["clusters"]))
chrom_cmap = dict(zip(config["chrom_order"], config["colors"]["chrom"]))

# I have a separate color scheme for boxplots, that does not contain Y.
chrom_boxplot_cmap = dict(zip(config["chrom_order"][:-1], config["colors"]["chrom_boxplot"]))


def make_ax(*args, **kwargs):
    fig, ax = plt.subplots(*args, **kwargs)
    return ax


def centerify(text, width=-1):
    """Center multiline text."""
    lines = text.split(" ")
    width = max(map(len, lines)) if width == -1 else width
    return "\n".join(line.center(width) for line in lines)
