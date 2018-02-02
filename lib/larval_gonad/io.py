"""Collection of io related items."""
from pathlib import Path
from joblib import Memory

cache_dir = Path(Path(__file__).parent, '../../output/cache').resolve()
cache_dir.mkdir(exist_ok=True)

memory = Memory(cachedir=cache_dir, verbose=0)
