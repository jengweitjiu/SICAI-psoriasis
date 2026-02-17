"""In silico L-R perturbation analysis."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm

def lr_perturbation(adata, fib_col, imm_col, ligand_genes, receptor_genes, lr_pairs, k_neighbours=6, coord_key='spatial'):
    from .coupling import _get_spatial_neighbours
    coords = adata.obsm[coord_key]
    nn = _get_spatial_neighbours(coords, k=k_neighbours)
    fib_vals = adata.obsm['cell2location'][:, list(adata.uns['cell2location_names']).index(fib_col)]
    imm_vals = adata.obs[imm_col].values
    imm_neigh = np.array([imm_vals[nn[i]].mean() for i in range(len(adata))])
    cs_full, _ = spearmanr(fib_vals, imm_neigh)
    contribs = {}
    for pair in lr_pairs:
        lig, rec = ligand_genes.get(pair), receptor_genes.get(pair)
        if not lig or not rec: contribs[pair] = 0.0; continue
        try:
            lx = adata[:,lig].X.toarray().ravel() if hasattr(adata[:,lig].X,'toarray') else adata[:,lig].X.ravel()
            rx = adata[:,rec].X.toarray().ravel() if hasattr(adata[:,rec].X,'toarray') else adata[:,rec].X.ravel()
            X = sm.add_constant(np.column_stack([lx, rx]))
            fr, ir = sm.OLS(fib_vals, X).fit().resid, sm.OLS(imm_neigh, X).fit().resid
            cp, _ = spearmanr(fr, ir)
            contribs[pair] = max(0, round((abs(cs_full)-abs(cp))/abs(cs_full)*100, 1)) if abs(cs_full)>0.01 else 0.0
        except: contribs[pair] = 0.0
    return pd.Series(contribs, name='contribution_pct')
