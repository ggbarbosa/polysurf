from __future__ import annotations
import numpy as np
from typing import Tuple, Dict
from ase import Atoms
from ase.build import make_supercell
from ase.io import read
from .utils import detect_chain_axis

class AlignmentError(RuntimeError):
    """Raised when no acceptable alignment is found."""
    pass

def _vectorized_stretch_search(
    slab_cell: np.ndarray,
    poly_cell: np.ndarray,
    poly_axis_idx: int,
    max_n_slab: int = 50,
    max_n_poly: int = 5,
    stretch_tol: Tuple[float, float] = (0.8, 1.25),
) -> Dict:
    """
    Vectorized search considering replication of both slab and polymer.
    Returns best candidate minimizing |stretch - 1|.
    """
    a_slab, b_slab = float(slab_cell[0, 0]), float(slab_cell[1, 1])
    a_poly, b_poly = float(poly_cell[0, 0]), float(poly_cell[1, 1])
    p_arr = np.array([a_poly, b_poly])
    m_arr = np.array([a_slab, b_slab])

    n_slab = np.arange(1, max_n_slab + 1)
    n_poly = np.arange(1, max_n_poly + 1)

    # stretch shape: (2, n_slab, n_poly)
    stretch = (n_slab[None, :, None] * m_arr[:, None, None]) / (
        n_poly[None, None, :] * p_arr[:, None, None]
    )

    valid_mask = (stretch > stretch_tol[0]) & (stretch < stretch_tol[1])
    candidates = np.argwhere(valid_mask)
    if candidates.size == 0:
        raise AlignmentError("No valid stretch/replication found within given tolerances.")

    # choose candidate with minimal absolute deviation from 1.0
    best_idx = min(candidates, key=lambda idx: abs(stretch[tuple(idx)] - 1.0))
    axis_idx, i_slab, i_poly = int(best_idx[0]), int(best_idx[1]), int(best_idx[2])

    return {
        "axis": "a" if axis_idx == 0 else "b",
        "stretch": float(stretch[axis_idx, i_slab, i_poly]),
        "n_slab": int(n_slab[i_slab]),
        "n_poly": int(n_poly[i_poly]),
    }

def build_slab_polymer_system(
    slab_path: str,
    polymer_path: str,
    v_step: float = 3.0,
    z_void: float = 18.0,
    lateral_min: float = 10.0,
    max_n_slab: int = 50,
    max_n_poly: int = 5,
    stretch_tol: Tuple[float, float] = (0.8, 1.25),
) -> Tuple[Atoms, Atoms, Atoms, Dict]:
    """Build and return (slab_super, polymer_replicated, combined, info).

    slab_path: path to slab structure file (ASE-readable)
    polymer_path: path to polymer structure file (ASE-readable)
    v_step: vertical offset (Å) between slab top and polymer bottom
    z_void: vacuum thickness (Å)
    lateral_min: minimum lateral spacing in x/y (Å)
    max_n_slab: maximum replication factor to try for slab along chain axis
    max_n_poly: maximum replication factor to try for polymer along chain axis
    stretch_tol: acceptable stretch interval (min, max)
    """
    slab = read(slab_path)
    polymer = read(polymer_path)

    # detect polymer axis (0->a,1->b)
    poly_axis_idx, diag = detect_chain_axis(polymer)

    match = _vectorized_stretch_search(
        slab.cell.array, polymer.cell.array, poly_axis_idx, max_n_slab, max_n_poly, stretch_tol
    )

    idx = 0 if match["axis"] == "a" else 1
    orth_idx = 1 - idx

    # replicate slab
    rep_matrix = np.identity(3, dtype=int)
    rep_matrix[idx, idx] = match["n_slab"]
    # ensure lateral_min (Å) spacing along orthogonal direction
    base_len = float(slab.cell.array[orth_idx, orth_idx])
    rep_orth = max(1, int(np.ceil(lateral_min / base_len)))
    # make rep_orth even to preserve symmetry if needed
    if rep_orth % 2 != 0:
        rep_orth += 1
    rep_matrix[orth_idx, orth_idx] = rep_orth

    slab_super = make_supercell(slab, rep_matrix)

    # replicate polymer
    rep_poly = np.identity(3, dtype=int)
    rep_poly[idx, idx] = match["n_poly"]
    polymer_replicated = make_supercell(polymer, rep_poly)

    # stretch polymer cell along chain axis
    poly_cell = polymer_replicated.cell.array.copy()
    poly_cell[idx, idx] *= match["stretch"]
    polymer_replicated.set_cell(poly_cell, scale_atoms=True)

    # center-of-mass alignment in x,y
    com_slab = slab_super.get_center_of_mass()
    com_poly = polymer_replicated.get_center_of_mass()
    polymer_replicated.translate([com_slab[0] - com_poly[0], com_slab[1] - com_poly[1], 0.0])

    # vertical positioning with v_step (Å)
    z_slab_max = float(np.max(slab_super.positions[:, 2]))
    z_poly_min = float(np.min(polymer_replicated.positions[:, 2]))
    polymer_replicated.translate([0.0, 0.0, z_slab_max - z_poly_min + v_step])

    # combined cell with vacuum z_void (Å)
    z_poly_max = float(np.max(polymer_replicated.positions[:, 2]))
    cell = slab_super.cell.array.copy()
    cell[2, 2] = z_poly_max - float(np.min(slab_super.positions[:, 2])) + z_void

    # shift so min z = 0
    z_min_total = min(np.min(slab_super.positions[:, 2]), np.min(polymer_replicated.positions[:, 2]))
    slab_super.translate([0.0, 0.0, -z_min_total])
    polymer_replicated.translate([0.0, 0.0, -z_min_total])

    slab_super.set_cell(cell, scale_atoms=False)
    polymer_replicated.set_cell(cell, scale_atoms=False)

    combined = slab_super + polymer_replicated

    info = {
        "detection": diag,
        "match": match,
        "rep_matrix": rep_matrix.tolist(),
        "lateral_min": float(lateral_min),
        "v_step": float(v_step),
        "z_void": float(z_void),
    }

    return slab_super, polymer_replicated, combined, info