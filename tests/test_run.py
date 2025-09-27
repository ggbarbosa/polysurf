import json
from ase.io import write
from polysurf import build_slab_polymer_system

# Input parameters
slab_file = "ti3c2o2_rec_ucell.xyz"
polymer_file = "pe_2_monomers.xyz"
chain_axis = "y"
v_step = 2.5
z_void = 15.0
lateral_min = 10.0

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
write("slab_final.xyz", slab)
write("polymer_final.xyz", polymer)
write("combined.xyz", combined)

# Save metadata
with open("metadata.json", "w") as f:
    json.dump(info, f, indent=2)
