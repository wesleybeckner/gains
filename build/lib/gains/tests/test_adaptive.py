import gains.adaptive as adl
import pandas as pd

viz = pd.read_csv("../../scripts/visualize/assets/vizapp.csv")
logs = "keras_8.0.0_salt_log.csv"
df = viz.loc[viz["Round"] == 5].reset_index(drop=True)
salt_log = pd.read_csv("../../scripts/salt_logs/{}".format(logs))
dens_model = adl.build_model_from_md(salt_log, "density")
cpt_model = adl.build_model_from_md(salt_log, "cpt")
