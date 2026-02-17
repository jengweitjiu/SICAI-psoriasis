"""
SICAI Coupling Metrics: CS, CSp, AS, and composite SICAI.
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.spatial import cKDTree

def _get_spatial_neighbours(coords, k=6):
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k+1)
    return indices[:, 1:]

def compute_cs(adata, fibroblast_key='cell2location', immune_keys=None, k_neighbours=6, coord_key='spatial'):
    coords = adata.obsm[coord_key]
    fib = pd.DataFrame(adata.obsm[fibroblast_key], columns=adata.uns.get(f'{fibroblast_key}_names', [f'F{i}' for i in range(adata.obsm[fibroblast_key].shape[1])]))
    imm = adata.obs[immune_keys].copy()
    nn = _get_spatial_neighbours(coords, k=k_neighbours)
    imm_neigh = pd.DataFrame(np.array([imm.iloc[nn[i]].mean(axis=0) for i in range(len(adata))]), columns=immune_keys)
    cs = pd.DataFrame(index=fib.columns, columns=immune_keys, dtype=float)
    for f in fib.columns:
        for im in immune_keys:
            rho, _ = spearmanr(fib[f].values, imm_neigh[im].values)
            cs.loc[f, im] = rho
    return cs

def compute_csp(cs_matrix):
    abs_cs = cs_matrix.abs()
    raw = abs_cs.div(abs_cs.sum(axis=1), axis=0)
    return raw.subtract(raw.mean(axis=1), axis=0).div(raw.std(axis=1), axis=0)

def compute_as(adata, fibroblast_key='cell2location', immune_keys=None, window_size=500, n_windows=20, k_neighbours=6, coord_key='spatial', seed=42):
    rng = np.random.RandomState(seed)
    coords = adata.obsm[coord_key]
    window_cs = []
    for _ in range(n_windows):
        cx = rng.uniform(coords[:,0].min()+window_size/2, coords[:,0].max()-window_size/2)
        cy = rng.uniform(coords[:,1].min()+window_size/2, coords[:,1].max()-window_size/2)
        mask = ((coords[:,0]-cx)**2 + (coords[:,1]-cy)**2) < (window_size/2)**2
        if mask.sum() < 20: continue
        window_cs.append(compute_cs(adata[mask].copy(), fibroblast_key, immune_keys, k_neighbours, coord_key))
    if len(window_cs) < 3:
        fib_names = adata.uns.get(f'{fibroblast_key}_names', [f'F{i}' for i in range(adata.obsm[fibroblast_key].shape[1])])
        return pd.DataFrame(0.5, index=fib_names, columns=immune_keys)
    stacked = np.stack([c.values for c in window_cs])
    means = np.abs(stacked.mean(axis=0))
    stds = stacked.std(axis=0)
    cv = np.divide(stds, means, where=means>0.01, out=np.ones_like(means))
    return pd.DataFrame(1-np.clip(cv,0,1), index=window_cs[0].index, columns=window_cs[0].columns)

def compute_sicai(cs, csp, as_mat):
    return cs * csp * as_mat
