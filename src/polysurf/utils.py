from __future__ import annotations
import numpy as np
from typing import Tuple
from ase import Atoms

def principal_axis_pca(positions: np.ndarray) -> np.ndarray:
    """Return the principal direction (unit vector) of a cloud of points using PCA."""
    if positions.shape[0] < 2:
        raise ValueError("At least two positions are required for PCA.")
    centered = positions - positions.mean(axis=0)
    cov = np.dot(centered.T, centered) / (centered.shape[0] - 1)
    eigvals, eigvecs = np.linalg.eigh(cov)
    principal = eigvecs[:, np.argmax(eigvals)]
    return principal / np.linalg.norm(principal)

def detect_chain_axis(atoms: Atoms) -> Tuple[int, dict]:
    """Detect the axis index (0 for a, 1 for b) along which the polymer chain grows most."""
    cell = atoms.get_cell()
    a_vec, b_vec = cell[0], cell[1]
    principal = principal_axis_pca(atoms.get_positions())
    a_hat = a_vec / np.linalg.norm(a_vec)
    b_hat = b_vec / np.linalg.norm(b_vec)
    proj_a = abs(float(np.dot(principal, a_hat)))
    proj_b = abs(float(np.dot(principal, b_hat)))
    extent_a = float(np.ptp(atoms.get_positions().dot(a_hat)))
    extent_b = float(np.ptp(atoms.get_positions().dot(b_hat)))
    axis = 0 if proj_a >= proj_b else 1
    info = {
        "proj_a": proj_a,
        "proj_b": proj_b,
        "extent_a": extent_a,
        "extent_b": extent_b,
        "principal_vector": principal.tolist(),
    }
    return axis, info