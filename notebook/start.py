# Imports
import os
import sys

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
sys.path.insert(0, '../lib')
from larval_gonad.notebook import setup_notebook

# Setup notebook
nbconfig = setup_notebook()

# Turn on cache
from joblib import Memory
memory = Memory(cachedir=nbconfig['CACHE'], verbose=0)
