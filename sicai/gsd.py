"""Geometric Stability Decomposition."""
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity

def fit_landscape(cs_vectors, condition_labels, bandwidth=0.4):
    pca = PCA()
    pc_scores = pca.fit_transform(cs_vectors.values)
    cumvar = np.cumsum(pca.explained_variance_ratio_)
    n_pcs = max(2, int(np.argmax(cumvar >= 0.80) + 1))
    pc2 = pc_scores[:,:2]
    kde = gaussian_kde(pc2.T, bw_method=bandwidth)
    m = 1.0
    xg = np.linspace(pc2[:,0].min()-m, pc2[:,0].max()+m, 100)
    yg = np.linspace(pc2[:,1].min()-m, pc2[:,1].max()+m, 100)
    xx, yy = np.meshgrid(xg, yg)
    density = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    return dict(pca=pca, pc_scores=pc_scores, pc2=pc2, kde=kde, conditions=np.array(condition_labels), grid_x=xx, grid_y=yy, potential=-np.log(density+1e-10), n_pcs=n_pcs)

def compute_basin_metrics(landscape, condition_labels=None):
    if condition_labels is None: condition_labels = landscape['conditions']
    pc2, conds = landscape['pc2'], np.array(condition_labels)
    results = []
    for cond in np.unique(conds):
        mask = conds == cond; pts = pc2[mask]; n = mask.sum()
        if n < 2: continue
        d = landscape['kde'](pts.T)
        depth = -np.log(d.min()+1e-10) - (-np.log(d.max()+1e-10))
        curvature = 1.0 / (np.mean(np.var(pts, axis=0)) + 1e-6)
        cs = cosine_similarity(landscape['pc_scores'][mask])
        coherence = (cs.sum()-n)/(n*(n-1)) if n>1 else 1.0
        other = pc2[~mask]
        escape = np.mean([np.min(np.linalg.norm(other-p, axis=1)) for p in pts]) if len(other)>0 else np.nan
        results.append(dict(condition=cond, n_samples=n, depth=round(depth,4), curvature=round(curvature,4), coherence=round(coherence,4), escape_energy=round(escape,4)))
    return pd.DataFrame(results)
