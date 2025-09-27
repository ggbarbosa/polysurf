#test/test_run.py
import json
from ase.io import write
from polysurf import build

# Input parameters
slab_file = "graphene_rectangular.xyz"
polymer_file = "pe_2monomers_linear.xyz"
chain_axis = "y"
v_step = 2.5
z_void = 15.0
lateral_min = 10.0

# Build system
slab, polymer, combined, info = build(
    slab_file,
    polymer_file,
    chain_axis,
    v_step,
    z_void,
    lateral_min
)

# Save outputs
write("slab_final.xyz", slab)
write("polymer_final.xyz", polymer)
write("combined.xyz", combined)

# Save metadata
with open("metadata.json", "w") as f:

    json.dump(info, f, indent=2)