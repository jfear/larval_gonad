import numpy as np

# Make a beep in the notebook
# >>> from IPython.display import Audio
# >>> Audio(**beep)
beep = {
    "data": np.sin(2 * np.pi * 400 * np.arange(10000 * 2) / 20000),
    "rate": 20000,
    "autoplay": True,
}
