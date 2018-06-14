import os
from pathlib import Path
from yaml import load

from joblib import Memory

# Useful directories
PROJECT_DIR = Path(__file__).absolute().parents[1].as_posix()
CONFIG_DIR = Path(PROJECT_DIR, 'config').as_posix()
CACHE_DIR = Path(PROJECT_DIR, 'output/cache').as_posix()

# Make missing dirs
Path(CACHE_DIR).mkdir(exist_ok=True, parents=True)
Path(CONFIG_DIR).mkdir(exist_ok=True, parents=True)

# Load config file
with Path(CONFIG_DIR, 'common.yaml').open() as fh:
    config = load(fh)

REFERENCES_DIR = config.get('references_dir',
                            os.environ.get('REFERENCES_DIR', None))

# Trun on caching
memory = Memory(cachedir=CACHE_DIR, verbose=0)
