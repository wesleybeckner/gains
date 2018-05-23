from __future__ import absolute_import, division, print_function
import gains as genetic
from rdkit.Chem import AllChem as Chem
from rdkit.ML.Descriptors.MoleculeDescriptors import\
    MolecularDescriptorCalculator as calculator
import numpy as np
import unittest
import datetime
from math import exp
import random
import salty


class GuessIonTests(unittest.TestCase):
    geneSet = genetic.generate_geneset()
    df = salty.load_data("cationInfo.csv")
    parent_candidates = df['smiles'].unique()
    df = salty.load_data("anionInfo.csv")
    df = df['smiles'].unique()
    random.seed(1234)
    ohPickMe = random.sample(range(df.shape[0]), 1)
    anion = Chem.MolFromSmiles(df[ohPickMe[0]])

    def test_1_density(self):
        target = 1250
        self.guess_password(target)

    def test_benchmark(self):
        genetic.Benchmark.run(self.test_1_density)

    def guess_password(self, target):
        startTime = datetime.datetime.now()

        def fnGetFitness(genes):
            return get_fitness(self.anion, genes, target)

        def fnDisplay(candidate, mutation):
            display(candidate, mutation, startTime)

        def fnShowIon(genes, target, mutation_attempts, sim_score,
                      molecular_relative):
            show_ion(genes, target, mutation_attempts, sim_score,
                     molecular_relative, self.anion)

        optimalFitness = 0.9
        best = genetic.get_best(fnGetFitness, optimalFitness,
                                self.geneSet, fnDisplay,
                                fnShowIon, target, self.parent_candidates,
                                seed=123)
        cation = best.Mol
        anion = self.anion
        # Uncomment PDB lines to wrote PDB file
        # cation = Chem.AddHs(best.Mol)
        # Chem.EmbedMolecule(cation, Chem.ETKDG())
        # Chem.UFFOptimizeMolecule(cation)
        # Chem.rdmolfiles.MolToPDBFile(cation, "cation_test.pdb")
        # anion = Chem.AddHs(self.anion)
        # Chem.EmbedMolecule(anion, Chem.ETKDG())
        # Chem.UFFOptimizeMolecule(anion)
        # Chem.rdmolfiles.MolToPDBFile(anion, "anion_test.pdb")
        return cation, anion


def display(candidate, mutation, startTime):
    timeDiff = datetime.datetime.now() - startTime
    print("{}\t{}\t{}\t{}".format(candidate.Genes, candidate.Fitness,
                                  mutation, timeDiff))


class prod_model():
    def __init__(self, coef_data, model):
        self.Coef_data = coef_data
        self.Model = model


def get_fitness(anion, genes, target):
    model_ID = "density"
    cation = Chem.MolFromSmiles(genes)
    model = genetic.load_data("{}_qspr.h5".format(model_ID), h5File=True)
    deslist = genetic.load_data("{}_desc.csv".format(model_ID))
    feature_vector = []
    with genetic.suppress_rdkit_sanity():
        for item in deslist:
            if "anion" in item:
                feature_vector.append(calculator([item.partition('-')
                                      [0]]).CalcDescriptors(anion)[0])
            elif "cation" in item:
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
    prediction = exp(model.predict(np.array(features_normalized).
                     reshape(1, -1))[0])
    error = abs((prediction - target) / target)
    return 1 - error, prediction


def show_ion(genes, target, mutation_attempts, sim_score, molecular_relative,
             anion):
    mol = Chem.MolFromSmiles(genes)
    fitness, mol_property = get_fitness(anion, genes, target)
    print("{}\t{}".format("number of atoms: ", mol.GetNumAtoms()))
    print("{}\t{}".format("mutation attempts: ", mutation_attempts))
    print("with density: \t\t{0:1.2f} (kg/m)".format(mol_property))
    print("similarity score:  {0:10.3f}".format(sim_score))
    print("{}\t{}\n".format("molecular relative: ", molecular_relative))


if __name__ == '__main__':
    unittest.main()
