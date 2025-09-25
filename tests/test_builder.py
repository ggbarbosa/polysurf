import numpy as np
from ase import Atoms
from polysurf.utils import detect_chain_axis
from polysurf.builder import _vectorized_stretch_search

def test_detect_axis():
    positions = np.array([[0.0, i * 1.5, 0.0] for i in range(5)])
    atoms = Atoms(symbols=['C'] * 5, positions=positions)
    atoms.set_cell([[10.0, 0, 0], [0, 10.0, 0], [0, 0, 20.0]])
    axis, info = detect_chain_axis(atoms)
    assert axis == 1
    assert info['extent_b'] > info['extent_a']

def test_vectorized_search():
    slab_cell = np.diag([5.0, 6.0, 20.0])
    poly_cell = np.diag([3.0, 2.0, 10.0])
    result = _vectorized_stretch_search(slab_cell, poly_cell, poly_axis_idx=1, max_n_slab=10, max_n_poly=4)
    assert 'n_slab' in result and 'n_poly' in result