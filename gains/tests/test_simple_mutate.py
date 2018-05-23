from __future__ import absolute_import, division, print_function
import gains.engine as genetic
import unittest
import numpy as np
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem as Chem


class GuessIonTests(unittest.TestCase):

    def test_1_model(self):
        def fnGetFitness(genes):
            return get_fitness(genes, target)
        target = "CCCC"
        parent_candidates = np.array(["CCCO"])
        geneSet = genetic.generate_geneset()
        optimalFitness, prediction = get_fitness(target, target)
        genetic.get_best(fnGetFitness, optimalFitness, geneSet,
                         display, result_display, target,
                         parent_candidates, seed=123)

    def test_benchmark(self):
        genetic.Benchmark.run(self.test_1_model)


def get_fitness(genes, target):
    ms = [Chem.MolFromSmiles(target), Chem.MolFromSmiles(genes)]
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    return DataStructs.FingerprintSimilarity(fps[0], fps[1]), None


def display(candidate, mutation):
    print("{:>20}{:>15}{:>15}".format(mutation, "{:3.4f}".
                                      format(candidate.Fitness),
                                      candidate.Genes))


def result_display(genes, target, mutation_attempts, sim_score,
                   molecular_relative):
    mol = Chem.MolFromSmiles(genes)
    print("{:>20}{:>15}".format("number of atoms:", mol.GetNumAtoms()))
    print("{:>20}{:>15}".format("mutation attempts:", mutation_attempts))


if __name__ == '__main__':
    unittest.main()
