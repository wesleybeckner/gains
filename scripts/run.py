from gains.salt_generator import generate_solvent
import salty
from random import randint

model_ID = ["cpt", "density"]
T = [298.1, 298.16]
P = [101, 102]
exp_data = ["cpt", "density"]
data = salty.aggregate_data(exp_data, T=T, P=P)
merged = salty.merge_duplicates(data)
to_hull = merged.iloc[:, 2:4]
target = [0, 0]
simplex_id = 1
token_id = randint(1000, 9999)

generate_solvent(target, model_ID, heavy_atom_limit=20, sim_bounds=[0.55, 1],
                 hits=10, write_file=False, hull=to_hull, simplex=simplex_id,
                 exp_data=data, verbose=0,
                 gen_token=token_id, hull_bounds=[0,.1], inner_search=False,
                 parent_cap=1, mutation_cap=1000)


