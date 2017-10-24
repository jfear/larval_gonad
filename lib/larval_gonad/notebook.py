"""Set of helper functions for the notebook."""
import os
from datetime import datetime
from yaml import load

import numpy as np
import seaborn as sns
from IPython import get_ipython

from .plotting import add_styles


def setup_notebook(watermark=True):
    """Helper function to consistently setup notebooks.

    Parameters
    ----------
    watermark : bool
        If truen then output watermark information.

    Returns
    ----------
    dict
        Useful config options.

    """
    # add notebook magics
    mgc = get_ipython().magic

    ## Activate the autoreload extension for easy reloading of external packages
    mgc('reload_ext autoreload')
    mgc('autoreload 2')

    ## Trun on the water mark
    if watermark:
        mgc('reload_ext watermark')
        mgc('watermark -u -d -g')

    ## Plot inline
    mgc('matplotlib inline')

    # Figure out current project folder
    prj = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
    cfg = os.path.join(prj, 'config')

    # Set-up useful configs
    with open(os.path.join(cfg, 'common.yml')) as fh:
        config = load(fh)

    ## Add useful directories
    config['PROJECT_DIR'] = prj

    config['FIGURES'] = os.path.join(prj, 'output/figures')
    os.makedirs(config['FIGURES'], exist_ok=True)

    config['TABLES'] = os.path.join(prj, 'output/tables')
    os.makedirs(config['TABLES'], exist_ok=True)

    config['CACHE'] = os.path.join(prj, 'output/cache')
    os.makedirs(config['CACHE'], exist_ok=True)

    ## Current date
    nn = datetime.now()
    config['date'] = nn.strftime("%Y-%m-%d")

    ## figure formats
    config['formats'] = ['png', 'pdf']

    ## common figure styles
    config['formats'] = ['notebook', 'paper', 'talk', 'poster']

    # Set up plotting
    add_styles(os.path.join(prj, 'config/stylelib'))
    sns.set_context('notebook')

    # Turn off scientific notation
    np.set_printoptions(precision=5, suppress=True)

    return config
