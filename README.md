# PolySurf

**PolySurf** is a open-source, small, focused Python package that builds slab + polymer systems using ASE and NumPy.

PolySurf's design principles:
- **Single responsibility:** only generates structures (ASE `Atoms`); no job submission or DFT inputs.
- **Framework-agnostic:** can be used standalone or integrated into workflows (AiiDA examples provided separately).
- **Robust matching:** vectorized search for best slab/polymer replication and stretch.

## Quick example (native usage)

```python
from polysurf.builder import build_slab_polymer_system
from ase.io import write
import json

# Input parameters
slab_file = 'slab.xyz'
polymer_file = 'polymer.xyz'
chain_axis = "x"   # 'x' or 'y'
v_step = 3.0       # Å
z_void = 15.0      # Å
lateral_min = 10.0 # Å
max_rep = 100 # default
stretch_ratio_range = (0.75, 1.25) # default

# Build system
slab, polymer, combined, info = build_slab_polymer_system(
    slab_file,
    polymer_file,
    chain_axis,
    v_step,
    z_void,
    lateral_min
)

# Save outputs
write("slab_super.xyz", slab)
write("polymer_stretched.xyz", polymer)
write("combined.xyz", combined)

# Save metadata
with open("metadata.json", "w") as f:
    json.dump(info, f, indent=2)

print("Test completed. Files generated: slab_super.xyz, polymer_stretched.xyz, combined.xyz, metadata.json")
```

## Integration with AiiDA (optional)

If you want to integrate PolySurf into AiiDA workflows, simply convert the ASE outputs to `StructureData` and the metadata dict to `Dict`. This example is intentionally outside the core package (PolySurf does not depend on AiiDA):

```python
# example_aiida_integration.py (external to the package)
from polysurf.builder import build_slab_polymer_system
from aiida.orm import StructureData, Dict

slab_super, polymer_stretched, combined, info = build_slab_polymer_system(
    slab_path="slab.xyz",
    polymer_path="polymer.xyz"
)

slab_node = StructureData(ase=slab_super).store()
polymer_node = StructureData(ase=polymer_stretched).store()
combined_node = StructureData(ase=combined).store()
info_node = Dict(dict=info).store()
```

## Running tests locally

```bash
pip install -e .
pip install pytest
pytest -v
```

## License

MIT
