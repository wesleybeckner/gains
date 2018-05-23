import gains as genetic
from gains.salt_generator import generate_solvent

target = [1000, 1000]
model_ID = ["cpt", "density"]
generate_solvent(target, model_ID, heavy_atom_limit=20, sim_bounds=[0.8, 1],
                 hits=5, write_file=True)


