from rdkit.Chem import AllChem as Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import\
    MolecularDescriptorCalculator as calculator
from rdkit.Chem.rdmolfiles import MolToPDBFile
from numpy import array
import datetime
import salty
from math import exp
import random
import gains as genetic
import numpy as np
import re
import pandas as pd
import subprocess
import os
import shutil
targets = [[ 377, 1162],
       [ 811, 1238],
       [ 947, 1157],
#       [ 436, 1202],
#       [ 873,  995],
       [ 627, 1390],
#       [ 889, 1408],
       [ 769, 1232],
       [ 975, 1024],
       [ 499, 1272]]
model_ID = "cpt_density_m1"
pdb_files = 25

def get_fitness(anion, genes, target, model_ID):
    cation = Chem.MolFromSmiles(genes)
    model = genetic.load_data("{}.sav".format(model_ID), pickleFile=True)
    deslist = genetic.load_data("{}_descriptors.csv".format(model_ID))
    feature_vector = []

    for item in deslist:

        if "anion" in item:
            with genetic.suppress_stdout_stderr():
                feature_vector.append(calculator([item.partition('-')
                                      [0]]).CalcDescriptors(anion)[0])
        elif "cation" in item:
            with genetic.suppress_stdout_stderr():
                feature_vector.append(calculator([item.partition('-')
                                      [0]]).CalcDescriptors(cation)[0])
        elif "Temperature, K" in item:
            feature_vector.append(298.15)
        elif "Pressure, kPa" in item:
            feature_vector.append(101.325)
        else:
            print("unknown descriptor in list: %s" % item)
    features_normalized = (feature_vector - deslist.iloc[0].values) /\
        deslist.iloc[1].values
    prediction = np.round(np.exp(model.predict(np.array(features_normalized).
                          reshape(1, -1))[0]), decimals=2)
    error = abs((prediction - target) / target)
    error = np.average(error)

    return 1 - error, prediction

def xerror(x, y):
  """
  returns the error
  of x in regard to y 
  """
  return round (abs((x-y)/y) * 100, 2)

def add_cpt_dens(targets):
  """
  loop through target list, use as handle
  read salt_log
  read current/salt_log.csv
  loop through created salts: current/systems/C{}_A{}/production_md/analysis/*dens and *cpt
  read cpt and aggregate density (w/ std)
  add to new salt_log file and old one?
  """
  for i in range(len(targets)):
    current="{}_{}".format(targets[i][0],targets[i][1])
    os.chdir(current)
    log = pd.read_csv("salt_log.csv")
    for j in range(log.shape[0]):
      try:
        salt_ID = log["Salt ID"][j]
        path = ("systems/{}/production_md/analysis/{}".format(salt_ID,salt_ID))
        with open("{}.cpt".format(path), "rt") as fin:
          for line in fin:
            md_cpt = float(re.findall("\d+\.\d+", line)[0])
        md_dens = np.average(pd.read_csv("{}.dens".format(path), header=None)[1])
        log.iloc[j,7] = "{:.2f}, {:.2f}".format(md_cpt, md_dens)
      except:
        pass
    log.to_csv("salt_log.csv", index=False)
    os.chdir("../")

def add_error(targets):
  """
  read salt logs
  calculate differences between pred and actual
  drop nan rows and append to master log
  save master log file
  """
  for i in range(len(targets)):
    current="{}_{}".format(targets[i][0],targets[i][1])
    os.chdir(current)
    log = pd.read_csv("salt_log.csv")
    for j in range(log.shape[0]):
      try:
        salt_ID = log["Salt ID"][j]
        predictions = log["Model Prediction"][j]
        calculations = log["MD Calculation"][j]
        predictions = re.findall("\d+\.\d+", predictions)
        calculations = re.findall("\d+\.\d+", calculations)
        errors = []
        for k in range(len(predictions)):
            errors.append(xerror(float(predictions[k]), float(calculations[k])))
        log.iloc[j,8] = "{}".format(errors)
      except:
        pass
    log.to_csv("salt_log.csv", index=False)
    os.chdir("../")

def drop_nan_and_agg(targets):
  """
  loop through target directories, grab salt_log,
  drop nans, append to growing dataframe.
  save dataframe in current directory.
  """
  out = pd.DataFrame()
  for i in range(len(targets)):
    current="{}_{}".format(targets[i][0],targets[i][1])
    os.chdir(current)
    log = pd.read_csv("salt_log.csv")
    log = log.dropna()
    out = pd.concat([log, out], axis=0)
    os.chdir("../")
  out.reset_index(drop=True, inplace=True)
  out.to_csv("salt_log.csv", index=False)

def add_pdb_to_salt_log(targets, model_ID, write_file=False):
  """
  loop through target directories, find completed
  salt simulations in systems/{}/production_md/analysis/{}
  and make sure they've made it to the salt_log
  """
  for i in range(len(targets)):
    current="{}_{}".format(targets[i][0],targets[i][1])
    os.chdir(current)
    target = targets[i]
    summary = genetic.load_data("{}_summary.csv".format(model_ID))
    parent_candidates = None
    if parent_candidates is None:
        parent_candidates = eval(summary.iloc[1][1])
    anion_candidates = eval(summary.iloc[2][1])
    cols = ["Salt ID", "Salt Smiles", "Cation Heavy Atoms",
            "Tanimoto Similarity Score", "Molecular Relative", "Anion",
            "Model Prediction", "MD Calculation", "Error"]
    salts = pd.DataFrame(columns=cols)

    for i in range(1, pdb_files + 1):
        try:
          if i < 10:
              CAT_ID = "C0%s" % i
              AN_ID = "A0%s" % i
          else:
              CAT_ID = "C%s" % i
              AN_ID = "A%s" % i
          salt_ID = CAT_ID + "_" + AN_ID
          path = ("systems/{}/production_md/analysis/{}".format(salt_ID,salt_ID))
          anion = Chem.MolFromPDBFile("structures/{}.pdb".format(AN_ID))
          anion_smiles = Chem.MolToSmiles(anion)
          cation = Chem.MolFromPDBFile("structures/{}.pdb".format(CAT_ID))
          best = genetic.Chromosome(Chem.MolToSmiles(cation), 0)
          tan_sim_score, sim_index =\
              genetic.molecular_similarity(best, parent_candidates)
          cation_heavy_atoms = best.Mol.GetNumAtoms()
          salt_smiles = best.Genes + "." + Chem.MolToSmiles(anion)

          scr, pre = get_fitness(anion, best.Genes, target, model_ID)
          molecular_relative = salty.check_name(parent_candidates
                                                [sim_index])
          anion_name = salty.check_name(anion_smiles)
	  #if anion_name == 0:
          anion_chromosome = genetic.Chromosome(anion_smiles, 0)
          find_name_score, find_name_index =\
              genetic.molecular_similarity(anion_chromosome, anion_candidates)
          print(find_name_score)
          anion_name = salty.check_name(anion_candidates
                                              [find_name_index])
          new_entry = pd.DataFrame([[salt_ID, salt_smiles,
                                     cation_heavy_atoms,
                                     tan_sim_score,
                                     molecular_relative,
                                     anion_name, pre]],
                                   columns=cols[:-2])
          new = pd.DataFrame(pd.concat([salts, new_entry]),
                               columns=cols)
          salts = new
        except:
          pass
    if write_file:
      pd.DataFrame.to_csv(new, path_or_buf="salt_log.csv", index=False)   
    else:
      print(new)
    os.chdir("../")
add_pdb_to_salt_log(targets,model_ID=model_ID,write_file=True)
add_cpt_dens(targets)
add_error(targets)
drop_nan_and_agg(targets)
