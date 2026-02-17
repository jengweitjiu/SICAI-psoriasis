"""
SICAI: Stromal-Immune Coupled Attractor Index

A computational framework for quantifying stromal-immune coupling 
architecture from spatial transcriptomics data.
"""

__version__ = "1.0.0"
__author__ = "[Authors]"

from .coupling import compute_cs, compute_csp, compute_as, compute_sicai
from .gsd import fit_landscape, compute_basin_metrics
from .perturbation import lr_perturbation
