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


def figure_element(func):
    """Decorate figure elements to create an Axes if none is given.

    I define a figure element as a function that imports, munges, and plots. A
    figure panel can be built by combining multiple figure elements using
    matplotlib.Axes. However, some times I will only want to plot the element,
    so this wrapper will generate a new Axes and pass it to the function for
    me.

    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        ax = kwargs.get("ax", None)
        if ax is None:
            fig, ax = plt.subplots()
            kwargs.update({"ax": ax})
        return func(*args, **kwargs)

    return wrapper


def centerify(text, width=-1):
    """Center multiline text."""
    lines = text.split(" ")
    width = max(map(len, lines)) if width == -1 else width
    return "\n".join(line.center(width) for line in lines)
