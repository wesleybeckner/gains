from __future__ import absolute_import, division, print_function
import gains as genetic
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem as Chem
import random
import unittest
import datetime
import salty


class GuessIonTests(unittest.TestCase):
    geneSet = genetic.generate_geneset()
    df = salty.load_data("cationInfo.csv")
    df = df.loc[df["name"].str.contains("imid", case=False)]
    df = df.loc[~df["name"].str.contains("phenyl", case=False)]
    df = df.loc[~df["name"].str.contains("benzyl", case=False)]
    df = df.loc[~df["name"].str.contains("azido", case=False)]
    df = df.loc[~df["name"].str.contains("cyan", case=False)]
    df = df.loc[~df["name"].str.contains("benz", case=False)]
    df = df.loc[~df["name"].str.contains("cyclo", case=False)]
    df = df.loc[~df["name"].str.contains("sulf", case=False)]
    df = df.loc[~df["name"].str.contains("azepinium", case=False)]
    parent_candidates = df['smiles'].unique()

    def test_2_similarity_map(self):
        best = genetic.Chromosome('CCCO', 0)
        genetic.molecular_similarity(best, self.parent_candidates, all=True)

    def test_1_similarity_map(self):
        df = self.parent_candidates
        random.seed(1234)
        ohPickMe = random.sample(range(df.shape[0]), 1)
        target = df[ohPickMe[0]]
        self.guess_password(target)

    def test_benchmark(self):
        genetic.Benchmark.run(self.test_1_similarity_map)
        genetic.Benchmark.run(self.test_2_similarity_map)

    def guess_password(self, target):
        startTime = datetime.datetime.now()

        def fnGetFitness(genes):
            return get_fitness(genes, target)

        def fnDisplay(candidate, mutation):
            display(candidate, mutation, startTime)

        def fnShowIon(genes, target, mutation_attempts, sim_score,
                      molecular_relative):
            show_ion(genes, target, mutation_attempts, sim_score,
                     molecular_relative)

        optimalFitness, prediction = get_fitness(target, target)
        best = genetic.get_best(fnGetFitness, optimalFitness, self.geneSet,
                                fnDisplay, fnShowIon, target,
                                self.parent_candidates, seed=123)
        return best


def display(candidate, mutation, startTime):
    timeDiff = datetime.datetime.now() - startTime
    print("{}\t{}\t{}\t{}".format(candidate.Genes, candidate.Fitness,
                                  mutation, timeDiff))


def get_fitness(genes, target):
    ms = [Chem.MolFromSmiles(target), Chem.MolFromSmiles(genes)]
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    return DataStructs.FingerprintSimilarity(fps[0], fps[1]), None


def show_ion(genes, target, mutation_attempts, sim_score, molecular_relative):
    mol = Chem.MolFromSmiles(genes)
    print("{}\t{}".format("number of atoms: ", mol.GetNumAtoms()))
    print("{}\t{}".format("mutation attempts: ", mutation_attempts))
    print("within 1%% of target density: %s (kg/m) " % target)
    print("{}\t{}".format("similarity score: ", sim_score))
    print("{}\t{}".format("with molecular relative: ", molecular_relative))


if __name__ == '__main__':
    unittest.main()
