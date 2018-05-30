import gains as genetic
from gains.salt_generator import generate_solvent
import salty

model_ID = ["cpt", "density"]
T = [298.1, 298.16] #select narrow state variable ranges
P = [101, 102] #we will set MD simulation to 101 kPa and 298 K
data = salty.aggregate_data(model_ID,T=T,P=P)
merged = salty.merge_duplicates(data)
to_hull = merged.iloc[:,2:4]

target = [1000, 1000]
generate_solvent(target, model_ID, heavy_atom_limit=20, sim_bounds=[0.8, 1],
                 hits=5, write_file=True, hull=to_hull)


