## Imports
import scanpy as sc
import numpy as np
import pandas as pd


DATA_PATH = './data/500.h5ad'
SAVE_DIR = './data/'


# single cell data
data = sc.read_h5ad(DATA_PATH)

id_train = data.uns['celmod']['avg.predicted.prop']['train'].index.astype(int)
id_train = id_train.sort_values()
id_shared = data.uns['celmod']['shared.donors'].astype(int)
id_shared.sort()

# Filtering to only the most confident cell states (58)
cell_corr_data = data.uns['celmod']['test.corrs']
cell_state_filter = cell_corr_data.index[(cell_corr_data['adj.pval'] < 0.005) & (cell_corr_data['corr'] > 0)]

# This is the bulk samples that also appears on the single cell data and filtered cell types
shared_bulk_data = data.uns['celmod']['avg.predicted.prop']['train'][cell_state_filter]

X = data.layers['sqrt.prev'] # single cell data
features_names = shared_bulk_data.columns
# first response vector, this is the probability of each sample to be part of the W/S branch
y1 = data.uns['trajectories']['palantir']['branch.probs']
# second response vector, this is the psuedotime of each sample
y2 = data.uns['trajectories']['palantir']['pseudotime']

mask = ~np.isnan(y2)
X_mask = X[mask]
y1_mask = y1[mask]
y2_mask = y2[mask]
print("Sainty check, num of null values: ",sum(np.isnan(y2_mask)))
print(X.shape, X_mask.shape)
# adding ID's as index for y2
y2_mask = pd.DataFrame(y2_mask)
y2_mask.index = y1_mask.index
y2_mask.rename(columns = {0:'Pseudotime'}, inplace=True)
IDs = y2_mask.index # This is the none null patients.
# filtering out of the valid patients the shared patients.
mask = []
for id in shared_bulk_data.index:
    if id in IDs:
        mask.append(id)

shared_bulk_data_mask = shared_bulk_data.loc[mask]
print(shared_bulk_data.shape, shared_bulk_data_mask.shape) # sanity check


# filtering the results vectors
new_y1_mask = y1.loc[mask]
new_y2_mask = pd.DataFrame(y2_mask).loc[mask]
print("Sainty check, num of null values: ", sum(np.isnan(new_y2_mask).astype(int)['Pseudotime']))
print(new_y1_mask.shape, new_y2_mask.shape)
new_y1_mask['psuedotime'] = new_y2_mask
y_mask = new_y1_mask

# save the data
np.save(SAVE_DIR + 'X.npy', X_mask)
y_mask.to_csv(SAVE_DIR + 'y.csv')
# save the features names
features_names = pd.DataFrame(features_names)
features_names.to_csv(SAVE_DIR + 'features_names.csv', index=False)
# save the shared bulk data
shared_bulk_data_mask.to_csv(SAVE_DIR + 'shared_bulk_data_mask.csv')

