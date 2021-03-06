#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'notebook'))
	print(os.getcwd())
except:
	pass

#%%
import numpy as np
import scipy.io
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')

import scrublet as scr


#%%
cell_ids = '../output/cellranger-wf/testis1/outs/raw_gene_bc_matrices/dmelr6-26/cell_ids.tsv'
mtx = '../output/cellranger-wf/testis1/outs/raw_gene_bc_matrices/dmelr6-26/matrix.mtx' 
cell_calls = '../output/cellselection-wf/testis1_combined_cell_calls.feather'


#%%
# Get list of all cell_ids
with open(cell_ids) as fh:
	cells = np.array(fh.read().strip().split("\n"))		# (737280,)
	print(cells.shape)


#%%
# Get a filtered set of cell calls
calls = pd.read_feather(cell_calls)		# (14169, 5)
filtered_cells = calls[calls.sum(axis=1) >= 3].cell_id.values		# (2717,)
print(filtered_cells.shape)

#%%
# Get a list of index locations of filtered cells
cell_index = np.in1d(cells, filtered_cells)		# (737280,)
print(cell_index.shape)

#%%
# Import matrix and filter cells
counts_matrix = scipy.io.mmread(mtx).tocsc()[:, cell_index].T		# (2717, 17508)
print(counts_matrix.shape)

#%%
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

#%%
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=20)
flag_doublet = scrub.call_doublets(threshold=0.25, verbose=True) # (2717,)

#%% 
scrub.plot_histogram()

#%%
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)

#%%
dubs = cells[cell_index][flag_doublet] 		# 42
print(dubs.shape)

#%%
scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, perplexity=30))
scrub.plot_embedding('tSNE', order_points=True)

#%%
