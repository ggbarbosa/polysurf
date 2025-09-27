import json
from ase.io import write
from polysurf import build_slab_polymer_system

# Input parameters
slab_file = "ti3c2o2_rec_ucell.xyz" # slab's 2D rectangular or square unit cell
polymer_file = "pe_2_monomers.xyz" # polymer's straight one-dimensional unit cell
chain_axis = "y"   # "x" or "y" polymer growth direction
v_step = 2.5       # (Å) polymer-slab vertical separation
z_void = 15.0      # (Å) vacuum in z
lateral_min = 10.0 # (Å) minimal lateral separation between polymer images

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

print("Test completed. Files generated: slab_final.xyz, polymer_final.xyz, combined.xyz, metadata.json")