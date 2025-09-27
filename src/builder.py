# polysurf/builder.py
import numpy as np
from typing import Dict, Tuple
from ase import Atoms
from ase.io import read
from ase.build import make_supercell
from scipy.spatial.distance import cdist


def min_polymer_image_distance(polymer: Atoms, cell: np.ndarray, orth_idx: int) -> float:
    """
    Compute the minimum distance between polymer atoms and their periodic
    images along the orthogonal axis.

    Parameters
    ----------
    polymer : Atoms
        ASE Atoms object of the polymer (already placed in slab cell).
    cell : np.ndarray
        3x3 array representing the cell vectors.
    orth_idx : int
        Index of the orthogonal axis (0 for x, 1 for y).

    Returns
    -------
    float
        Minimum distance (in Angstroms) between polymer atoms and their
        periodic image along the orthogonal axis.
    """
    coords = polymer.positions
    shift_vec = np.zeros(3)
    shift_vec[orth_idx] = cell[orth_idx, orth_idx]

    coords_image = coords + shift_vec
    return float(np.min(cdist(coords, coords_image)))


def builder(
    slab_path: str,
    polymer_path: str,
    chain_axis: str,
    v_step: float,
    z_void: float,
    lateral_min: float,
    max_rep: int = 100,
    stretch_ratio_range: Tuple[float, float] = (0.75, 1.25),
) -> Tuple[Atoms, Atoms, Atoms, Dict]:
    """
    Build a combined slab + polymer system with automatic stretch,
    vertical placement, and lateral separation.

    Parameters
    ----------
    slab_path : str
        Path to the slab .xyz file.
    polymer_path : str
        Path to the polymer .xyz file.
    chain_axis : str
        Growth axis of the polymer chain ("x" or "y").
    v_step : float
        Vertical separation (Angstroms) between polymer and slab surface.
    z_void : float
        Extra vacuum (Angstroms) added along z.
    lateral_min : float
        Minimum lateral separation (Angstroms) between polymer and its periodic image.
    max_rep : int, optional
        Maximum allowed replication along the orthogonal axis (default: 100).
    stretch_ratio_range : Tuple[float, float], optional
        Allowed stretch ratio range (default: (0.75, 1.25)).

    Returns
    -------
    slab_super : Atoms
        Final slab after replication.
    polymer_final : Atoms
        Final stretched and centered polymer.
    combined : Atoms
        Combined slab + polymer system.
    info : dict
        Metadata about the final structure.
    """
    # Load input structures
    slab = read(slab_path)
    polymer = read(polymer_path)

    # Axis mapping
    axis_map = {"x": 0, "y": 1}
    chain_idx = axis_map[chain_axis]
    orth_idx = 1 if chain_idx == 0 else 0

    # Find best replication along chain axis and corresponding stretch
    rep_chain = 1
    best_stretch = None
    poly_len = polymer.cell[chain_idx, chain_idx]
    while rep_chain <= max_rep:
        rm = np.eye(3, dtype=int)
        rm[chain_idx, chain_idx] = rep_chain
        slab_candidate = make_supercell(slab, rm)

        slab_len = slab_candidate.cell[chain_idx, chain_idx]
        stretch_ratio = slab_len / poly_len

        if stretch_ratio_range[0] <= stretch_ratio <= stretch_ratio_range[1]:
            best_stretch = (stretch_ratio - 1.0) * 100.0
            slab_chain = slab_candidate
            break
        rep_chain += 1

    if best_stretch is None:
        raise RuntimeError("Unable to match slab and polymer within stretch ratio range.")

    # Apply stretch to polymer
    polymer_stretched = polymer.copy()
    polymer_stretched.positions[:, chain_idx] *= (1.0 + best_stretch / 100.0)

    # Replicate orthogonal axis until minimum polymer separation is satisfied
    rep_orth = 1
    slab_super, polymer_centered = None, None
    while rep_orth <= max_rep:
        rm = np.eye(3, dtype=int)
        rm[chain_idx, chain_idx] = rep_chain
        rm[orth_idx, orth_idx] = rep_orth
        slab_candidate = make_supercell(slab, rm)

        # Place polymer in this candidate cell
        poly_temp = polymer_stretched.copy()
        slab_com = slab_candidate.get_center_of_mass()
        poly_com = poly_temp.get_center_of_mass()
        poly_temp.translate([slab_com[0] - poly_com[0], slab_com[1] - poly_com[1], 0.0])

        # Check separation
        lateral_gap = min_polymer_image_distance(poly_temp, slab_candidate.cell.array, orth_idx)

        if lateral_gap >= lateral_min:
            slab_super = slab_candidate
            polymer_centered = poly_temp
            break
        rep_orth += 1

    if slab_super is None:
        raise RuntimeError("Unable to satisfy lateral separation condition.")

    # Vertical placement
    slab_zmax = np.max(slab_super.positions[:, 2])
    polymer_zmin = np.min(polymer_centered.positions[:, 2])
    polymer_centered.translate([0.0, 0.0, slab_zmax - polymer_zmin + v_step])

    # Combine slab + polymer
    combined = slab_super.copy()
    combined.extend(polymer_centered)

    # Shift minimum z to 0
    min_z = np.min(combined.positions[:, 2])
    combined.positions[:, 2] -= min_z

    # Extend vacuum along z
    cell = combined.cell.copy()
    z_span = np.ptp(combined.positions[:, 2])
    cell[2, 2] = z_span + z_void
    combined.set_cell(cell, scale_atoms=False)

    # Final supercell to ensure consistency
    slab_super.set_cell(cell, scale_atoms=False)
    polymer_centered.set_cell(cell, scale_atoms=False)
    slab_super.positions[:, 2] -= min_z
    polymer_centered.positions[:, 2] -= min_z

    # Build metadata
    cell = combined.cell
    info = {
        "slab": slab_path,
        "polymer": polymer_path,
        "a_cell_A": round(cell[0, 0], 3),
        "b_cell_A": round(cell[1, 1], 3),
        "c_cell_A": round(cell[2, 2], 3),
        "m": rep_chain,
        "n": rep_orth,
        "nat": len(combined),
        "chain_axis": chain_axis,
        "stretch_pct": round(best_stretch, 2),
        "v_step_A": round(v_step, 3),
        "z_void_A": round(z_void, 3),
        "lateral_min_A": round(lateral_min, 3),
    }

    return slab_super, polymer_centered, combined, info
