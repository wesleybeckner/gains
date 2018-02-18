from __future__ import absolute_import, division, print_function
import gains as genetic
from rdkit.Chem import AllChem as Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import\
    MolecularDescriptorCalculator as calculator
from rdkit.Chem.rdmolfiles import MolToPDBFile
import numpy as np
from numpy import array
import pandas as pd
import datetime
import salty
from math import exp
import random


def guess_password(target, anion, parent_candidates, model_ID):
    startTime = datetime.datetime.now()

    def fnGetFitness(genes):
        return get_fitness(anion, genes, target, model_ID)

    def fnDisplay(candidate, mutation):
        display(candidate, mutation, startTime)

    def fnShowIon(genes, target, mutation_attempts, sim_score,
                  molecular_relative):
        show_ion(genes, target, mutation_attempts, sim_score,
                 molecular_relative, model_ID, anion)

    optimalFitness = 0.99
    geneSet = genetic.generate_geneset()
    best = genetic.get_best(fnGetFitness,
                            optimalFitness, geneSet, fnDisplay,
                            fnShowIon, target, parent_candidates)
    return best


def display(candidate, mutation, startTime):
    timeDiff = datetime.datetime.now() - startTime
    print("{}\t{}\t{}".format(
        candidate.Genes, candidate.Fitness, mutation, timeDiff))


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
    prediction = np.round(exp(model.predict(np.array(features_normalized).
                          reshape(1, -1))[0]), decimals=2)
    error = abs((prediction - target) / target)

    return 1 - error, prediction


def show_ion(genes, target, mutation_attempts, sim_score, molecular_relative,
             model_ID, anion):
    mol = Chem.MolFromSmiles(genes)
    fitness, mol_property = get_fitness(anion, genes, target, model_ID)
    print("{}\t{}".format("number of atoms: ", mol.GetNumAtoms()))
    print("{}\t{}".format("mutation attempts: ", mutation_attempts))
    print("with density: \t\t{0:1.2f} (kg/m)".format(mol_property))
    print("similarity score:  {0:10.3f}".format(sim_score))
    print("{}\t{}\n".format("molecular relative: ",
          salty.check_name(molecular_relative)))


def generate_solvent(target, model_ID, heavy_atom_limit=30,
                     sim_bounds=[0.6, 1.0], hits=1, write_file=False):

    parent_candidates = eval(genetic.load_data("{}_summary.csv".
                             format(model_ID)).loc[1][1])
    anion_candidates = eval(genetic.load_data("{}_summary.csv".
                            format(model_ID)).loc[2][1])
    cols = ["Salt ID", "Salt Smiles", "Cation Heavy Atoms",
            "Tanimoto Similarity Score", "Molecular Relative", "Anion",
            "Model Density", "MD Density", "Error"]
    salts = pd.DataFrame(columns=cols)
    for i in range(1, hits + 1):
        while True:
            anion_smiles = random.sample(list(anion_candidates), 1)[0]
            anion = Chem.MolFromSmiles(anion_smiles)
            best = guess_password(target, anion, parent_candidates, model_ID)
            tan_sim_score, sim_index =\
                genetic.molecular_similarity(best, parent_candidates)
            cation_heavy_atoms = best.Mol.GetNumAtoms()
            salt_smiles = best.Genes + "." + Chem.MolToSmiles(anion)
            if cation_heavy_atoms < heavy_atom_limit and\
                    tan_sim_score >= sim_bounds[0] and\
                    tan_sim_score < sim_bounds[1] and\
                    salt_smiles not in salts["Salt Smiles"]:
                scr, pre = get_fitness(anion, best.Genes, target, model_ID)
                if i < 10:
                    CAT_ID = "C0%s" % i
                    AN_ID = "A0%s" % i
                else:
                    CAT_ID = "C%s" % i
                    AN_ID = "A%s" % i
                salt_ID = CAT_ID + "_" + AN_ID
                molecular_relative = salty.check_name(parent_candidates
                                                      [sim_index])
                anion_name = salty.check_name(anion_smiles)
                new_entry = pd.DataFrame([[salt_ID, salt_smiles,
                                           cation_heavy_atoms,
                                           tan_sim_score,
                                           molecular_relative,
                                           anion_name, pre]],
                                         columns=cols[:-2])
                try:
                    cation = Chem.AddHs(best.Mol)
                    Chem.EmbedMolecule(cation, Chem.ETKDG())
                    Chem.UFFOptimizeMolecule(cation)
                    anion = Chem.AddHs(anion)
                    Chem.EmbedMolecule(anion, Chem.ETKDG())
                    Chem.UFFOptimizeMolecule(anion)
                    new = pd.DataFrame(pd.concat([salts, new_entry]),
                                       columns=cols)
                except BaseException:
                    continue
                if write_file:
                    MolToPDBFile(cation,
                                 "{}.pdb".format(CAT_ID))
                    MolToPDBFile(anion,
                                 "{}.pdb".format(AN_ID))
                break
            else:
                continue
        if write_file:
            pd.DataFrame.to_csv(new, path_or_buf="salt_log.csv", index=False)
        salts = new
    if not write_file:
        return new
