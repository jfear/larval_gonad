"""Collection of io related items."""
from pathlib import Path
from joblib import Memory
from yaml import load

# Create Cache
cache_dir = Path(Path(__file__).parent, '../../output/cache').resolve()
cache_dir.mkdir(exist_ok=True)
memory = Memory(cachedir=cache_dir, verbose=0)

# Get config
cname = Path(Path(__file__).parent, '../../config/common.yml').resolve()
with cname.open() as fh:
    config = load(fh)
