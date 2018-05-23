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


def generate_solvent(target, model_ID, heavy_atom_limit=50,
                     sim_bounds=[0.4, 1.0], hits=1, write_file=False,
                     seed=None):
    """
    the primary public function of the salt_generator module

    Parameters
    ----------
    target : array, float, or int
        the desired property value to be achieved by the engine, if
        an array, a multi-output model must be supplied to the engine
    model_ID : str
        the name of the model to be used by the engine. Gains has
        several built-in models to choose from
    heavy_atom_limit : int, optional
        the upper value for allowable heavy atoms in the returned
        candidate
    sim_bounds : array, optional
        the tanimoto similarity score between the returned candidate
        and its closest molecular relative in parent_candidates
    hits : int, optional
        the number of desired solutions
    write_file : boolean, optional
        defaults to False. if True will return the solutions and a
        csv log file

    Returns
    -------
    new : object
        default behavior is to return a pandas DataFrame. This is
        a log file of the solution(s). if write_file = True the
        function will also return pdb files of the cations/anions
    """
    parent_candidates = []
    anion_candidates = []
    for i, name in enumerate(model_ID):
        summary = genetic.load_data("{}_summ.csv".format(name))
        parents = eval(summary.iloc[1][1])
        anions = eval(summary.iloc[2][1])
        if i > 0:
            parent_candidates = np.concatenate((parents, parent_candidates))
            anion_candidates = np.concatenate((anions, anion_candidates))
        else:
            parent_candidates = parents
            anion_candidates = anions
    cols = ["Salt ID", "Salt Smiles", "Cation Heavy Atoms",
            "Tanimoto Similarity Score", "Molecular Relative", "Anion",
            "Model Prediction", "MD Calculation", "Error"]
    salts = pd.DataFrame(columns=cols)
    for i in range(1, hits + 1):
        while True:
            if seed:
                random.seed(seed)
            anion_smiles = random.sample(list(anion_candidates), 1)[0]
            anion = Chem.MolFromSmiles(anion_smiles)
            best = _guess_password(target, anion, parent_candidates, model_ID,
                                   seed=seed)
            tan_sim_score, sim_index =\
                genetic.molecular_similarity(best, parent_candidates)
            cation_heavy_atoms = best.Mol.GetNumAtoms()
            salt_smiles = best.Genes + "." + Chem.MolToSmiles(anion)
            if cation_heavy_atoms < heavy_atom_limit and\
                    tan_sim_score >= sim_bounds[0] and\
                    tan_sim_score < sim_bounds[1] and\
                    salt_smiles not in salts["Salt Smiles"]:
                scr, pre = _get_fitness(anion, best.Genes, target, model_ID)
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


def _guess_password(target, anion, parent_candidates, model_ID, seed=None):
    """
    for interacting with the main engine. Contains helper functions
    to pass to the engine what it expects
    """
    startTime = datetime.datetime.now()

    def fnGetFitness(genes):
        return _get_fitness(anion, genes, target, model_ID)

    def fndisplay(candidate, mutation):
        _display(candidate, mutation, startTime)

    def fnShowIon(genes, target, mutation_attempts, sim_score,
                  molecular_relative):
        _show_ion(genes, target, mutation_attempts, sim_score,
                  molecular_relative, model_ID, anion)

    optimalFitness = 0.99
    geneSet = genetic.generate_geneset()
    if seed:
        best = genetic.get_best(fnGetFitness, optimalFitness, geneSet,
                                fndisplay, fnShowIon, target,
                                parent_candidates, seed=seed)
    return best


def _display(candidate, mutation, startTime):
    """
    for printing results to the screen. _display is called for every
    accepted mutation
    """
    timeDiff = datetime.datetime.now() - startTime
    print("{}\t{}\t{}".format(
        candidate.Genes, candidate.Fitness, mutation, timeDiff))


def _get_fitness(anion, genes, target, model_ID):
    """
    the fitness function passed to the engine. In this case fitness
    is determined by a model developed by the salty module.
    The fitness function can handle multi-output models
    """
    predictions = []
    for i, name in enumerate(model_ID):
        cation = Chem.MolFromSmiles(genes)
        model = genetic.load_data("{}_qspr.h5".format(name), h5File=True)
        deslist = genetic.load_data("{}_desc.csv".format(name))
        feature_vector = []

        for item in deslist:

            if "anion" in item:
                with genetic.suppress_rdkit_sanity():
                    feature_vector.append(calculator([item.partition('-')
                                          [0]]).CalcDescriptors(anion)[0])
            elif "cation" in item:
                with genetic.suppress_rdkit_sanity():
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
        prediction = np.round(np.exp(model.predict(np.array(
                              features_normalized).reshape(1, -1))[0]),
                              decimals=2)
        predictions.append(prediction[0])
    predictions = np.array(predictions)
    error = abs((predictions - target) / target)
    error = np.average(error)

    return 1 - error, predictions


def _show_ion(genes, target, mutation_attempts, sim_score, molecular_relative,
              model_ID, anion):
    """
    for printing results to the screen. _show_ion is called when a candidate
    has achieved the desired fitness core and is returned by the engine
    """
    mol = Chem.MolFromSmiles(genes)
    fitness, mol_property = _get_fitness(anion, genes, target, model_ID)
    print("{}\t{}".format("number of atoms: ", mol.GetNumAtoms()))
    print("{}\t{}".format("mutation attempts: ", mutation_attempts))
    print("with prediction: \t{}".format(mol_property))
    print("similarity score:  {0:10.3f}".format(sim_score))
    print("{}\t{}\n".format("molecular relative: ",
          salty.check_name(molecular_relative)))
