import gains as genetic
from gains.salt_generator import generate_solvent
import salty

model_ID = ["cpt_1", "density_1"]
T = [298.1, 298.16] # select narrow state variable ranges
P = [101, 102] # we will set MD simulation to 101 kPa and 298 K
exp_data = ["cpt", "density"]
data = salty.aggregate_data(exp_data,T=T,P=P)
merged = salty.merge_duplicates(data)
to_hull = merged.iloc[:,2:4]
target = [1000, 1000]

### need to figure out convex strat
### something like convex_ID = some number
### that is attached during setup.py

simplex_id = 0 # from 0 to len(simplices)

####hull will now have to include the simplex ID

generate_solvent(target, model_ID, heavy_atom_limit=20, sim_bounds=[0.8, 1],
                 hits=5, write_file=True, hull=to_hull, simplex=simplex_id,
                 path='/home/wesleybeckner/Dropbox/Python/py3/gains/scripts'
                      '/adl_models', exp_data=data)


