from __future__ import absolute_import, division, print_function
from os.path import join
from keras.models import load_model
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
                     sim_bounds=[0, 1.0], hits=1, write_file=False,
                     seed=None, hull=None, simplex=None, path=None,
                     exp_data=None, verbose=0, gen_token=False,
                     hull_bounds=[0, 1], inner_search=True, parent_cap=25,
                     mutation_cap=1000):
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
    seed : int, optional
        optional randint seed for unittest consistency
    hull : pandas DataFrame, optional
        nxm pandas DataFrame to use convex hull search strategy. hull
        columns should be the same properties used in the genetic algorithm
        fitness test
    simplex : array, optional
        array to access boundary datapoints in the convex hull. This is used
        during target resampling defined by the convex hull/simplex
    path : str, optional
        absolute path to the qspr model used as the fitness function
    exp_data: salty devmodel obj, optional
        used during hull target reassignment search strategy. Salty devmodel
        object of the original experimental data
    verbose : int, optional, default 0
        0 : most verbose. Best child, parent/target resampling,
            sanitization failure
        1 : parent/target resampling, solution metadata, sanitization failure
        2 : solution metdata, sanitization failure
        3 : target resampling, csv-formatted solution metadata
        4 : csv-formatted solution metadata
    gen_token : int, str, optional
        a string or integer to append to file outputs. Useful in the case of
        parallel searches.
    hull_bounds : array, optional
        if hull and simplex are not none, hull_bounds describes the
        proximity convex_search should be to the simplex
    inner_search : bool, optional
        if hull and simplex are not none, inner_search specifies if
        convex_search should return values only within the convex hull

    Returns
    -------
    new : object
        default behavior is to return a pandas DataFrame. This is
        a log file of the solution(s). if write_file = True the
        function will also return pdb files of the cations/anions
    """
    parent_candidates = []
    anion_candidates = []
    models = []
    deslists = []
    for i, name in enumerate(model_ID):
        if path:
            model = np.array([load_model(join(path,
                                              '{}_qspr.h5'.format(name)))])
            with open(join(path, '{}_desc.csv'.format(name)),
                      'rb') as csv_file:
                deslist = list([pd.read_csv(csv_file, encoding='latin1')])
            with open(join(path, '{}_summ.csv'.format(name)),
                      'rb') as csv_file:
                summary = pd.read_csv(csv_file, encoding='latin1')
        else:
            model = np.array([genetic.load_data("{}_qspr.h5".format(name),
                                                h5File=True)])
            deslist = list([genetic.load_data("{}_desc.csv".format(name))])
            summary = genetic.load_data("{}_summ.csv".format(name))
        parents = eval(summary.iloc[1][1])
        anions = eval(summary.iloc[2][1])
        if i > 0:
            parent_candidates = np.concatenate((parents, parent_candidates))
            anion_candidates = np.concatenate((anions, anion_candidates))
            models = np.concatenate((model, models))
            deslists = list([deslist, deslists])
        else:
            parent_candidates = parents
            anion_candidates = anions
            models = model
            deslists = deslist
    cols = ["Salt ID", "Salt Smiles", "Cation Heavy Atoms",
            "Tanimoto Similarity Score", "Molecular Relative", "Anion",
            "Model Prediction", "MD Calculation", "Error"]
    salts = pd.DataFrame(columns=cols)
    if exp_data:
        anion_candidates = eval(exp_data.Data_summary.iloc[2][0])
    for i in range(1, hits + 1):
        while True:
            if seed:
                random.seed(seed)
            anion_smiles = random.sample(list(anion_candidates), 1)[0]
            anion = Chem.MolFromSmiles(anion_smiles)
            best = _guess_password(target, anion_smiles, parent_candidates,
                                   models, deslists, seed=seed, hull=hull,
                                   simplex=simplex, exp_data=exp_data,
                                   verbose=verbose, hull_bounds=hull_bounds,
                                   inner_search=inner_search,
                                   parent_cap=parent_cap,
                                   mutation_cap=mutation_cap)
            if exp_data:
                exp_parent_candidates = eval(exp_data.Data_summary.iloc[1][0])
                tan_sim_score, sim_index = \
                    genetic.molecular_similarity(best, exp_parent_candidates)
            else:
                tan_sim_score, sim_index = \
                    genetic.molecular_similarity(best, parent_candidates)
            cation_heavy_atoms = best.Mol.GetNumAtoms()
            salt_smiles = best.Genes + "." + Chem.MolToSmiles(anion)
            if cation_heavy_atoms < heavy_atom_limit and \
                    sim_bounds[0] <= tan_sim_score < sim_bounds[1] and\
                    salt_smiles not in salts["Salt Smiles"]:
                scr, pre = _get_fitness(anion, best.Genes, target, models,
                                        deslists)
                if i < 10:
                    CAT_ID = "C0%s" % i
                    AN_ID = "A0%s" % i
                else:
                    CAT_ID = "C%s" % i
                    AN_ID = "A%s" % i
                salt_ID = CAT_ID + "_" + AN_ID
                if exp_data:
                    molecular_relative = salty.check_name(exp_parent_candidates
                                                          [sim_index])
                else:
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
                    if verbose == any([0, 1, 2]):
                        print("molecule not sanitizable")
                    continue
                if write_file:
                    if verbose == any([3, 4]):
                        print(new)
                    if gen_token:
                        MolToPDBFile(cation,
                                     "{}_{}.pdb".format(gen_token, CAT_ID))
                        MolToPDBFile(anion,
                                     "{}_{}.pdb".format(gen_token, AN_ID))
                    else:
                        MolToPDBFile(cation,
                                     "{}.pdb".format(CAT_ID))
                        MolToPDBFile(anion,
                                     "{}.pdb".format(AN_ID))
                break
            else:
                continue
        if write_file:
            if gen_token:
                pd.DataFrame.to_csv(new, path_or_buf="{}_salt_log.csv".
                                    format(gen_token), index=False)
            else:
                pd.DataFrame.to_csv(new, path_or_buf="salt_log.csv",
                                    index=False)
        salts = new
    if not write_file:
        return new


def _guess_password(target, anion_smiles, parent_candidates, models, deslists,
                    seed=None, hull=None, simplex=None, exp_data=None,
                    verbose=0, hull_bounds=[0, 1], inner_search=True,
                    parent_cap=25, mutation_cap=1000):
    """
    for interacting with the main engine. Contains helper functions
    to pass to the engine what it expects
    """
    startTime = datetime.datetime.now()
    anion = Chem.MolFromSmiles(anion_smiles)

    def fnGetFitness(genes):
        return _get_fitness(anion, genes, target, models, deslists)

    def fndisplay(candidate, mutation):
        _display(candidate, mutation, startTime)

    def fnShowIon(genes, target, mutation_attempts, sim_score,
                  molecular_relative):
        _show_ion(genes, target, mutation_attempts, sim_score,
                  molecular_relative, models, deslists, anion_smiles,
                  exp_data)

    optimalFitness = 0.95
    geneSet = genetic.generate_geneset()
    best = genetic.get_best(fnGetFitness, optimalFitness, geneSet,
                            fndisplay, fnShowIon, target,
                            parent_candidates, seed=seed,
                            hull=hull, simplex=simplex,
                            verbose=verbose, hull_bounds=hull_bounds,
                            inner_search=inner_search,
                            parent_cap=parent_cap, mutation_cap=mutation_cap)
    return best


def _display(candidate, mutation, startTime):
    """
    for printing results to the screen. _display is called for every
    accepted mutation
    """
    timeDiff = datetime.datetime.now() - startTime
    print("{}\t{}\t{}".format(
        candidate.Genes, candidate.Fitness, mutation, timeDiff))


def _get_fitness(anion, genes, target, models, deslists):
    """
    the fitness function passed to the engine. In this case fitness
    is determined by a model developed by the salty module.
    The fitness function can handle multi-output models
    """
    predictions = []
    for i, name in enumerate(models):
        cation = Chem.MolFromSmiles(genes)
        model = name
        deslist = deslists[i]
        if isinstance(deslist, list):
            deslist = deslist[0]
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
    predictions = np.flip(np.array(predictions), axis=0)
    error = abs((predictions - target) / target)
    error = np.average(error)

    return 1 - error, predictions


def _show_ion(genes, target, mutation_attempts, sim_score, molecular_relative,
              models, deslists, anion_smiles, exp_data=None):
    """
    for printing results to the screen. _show_ion is called when a candidate
    has achieved the desired fitness core and is returned by the engine
    """
    mol = Chem.MolFromSmiles(genes)
    anion = Chem.MolFromSmiles(anion_smiles)
    fitness, mol_property = _get_fitness(anion, genes, target,
                                         models, deslists)
    anion_name = salty.check_name(anion_smiles)
    if exp_data:
        chrom = genetic.Chromosome(genes, fitness)
        exp_parent_candidates = eval(exp_data.Data_summary.iloc[1][0])
        tan_sim_score, sim_index = \
            genetic.molecular_similarity(chrom, exp_parent_candidates)
        molecular_relative = exp_parent_candidates[sim_index]
    print("{}\t{}".format("Salt Smiles: ", genes))
    print("{}\t{}".format("Cation Heavy Atoms: ", mol.GetNumAtoms()))
    print("Tanimoto Similarity Score: \t{0:10.3f}".format(sim_score))
    print("{}\t{}".format("Molecular Relative: ",
                          salty.check_name(molecular_relative)))
    print("{}\t{}".format("Anion: ", anion_name))
    print("{}\t{}".format("Model Prediction: ", mol_property))
    print("{}\t{}".format("Mutation Attempts: ", mutation_attempts))
