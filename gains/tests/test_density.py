from __future__ import absolute_import, division, print_function
import os
import sys
import gains as genetic
import rdkit
from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog
from rdkit import RDConfig
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import ShowMol
from rdkit.ML.Descriptors.MoleculeDescriptors import\
    MolecularDescriptorCalculator as calculator
import pandas as pd
import numpy as np
import random
import unittest
import datetime
from math import exp

class GuessIonTests(unittest.TestCase):
    geneSet = genetic.generate_geneset()
    
    def test_1_density(self):
        df = genetic.load_data("saltInfo.csv")
        df = df.loc[df["cation_name"].str.contains("imid", case=False)]
        df = df['cation_SMILES'].unique()
        ohPickMe = random.sample(range(df.shape[0]),1)
        target = 800
        self.guess_password(target)

    def test_benchmark(self):
        genetic.Benchmark.run(self.test_1_density)

    def guess_password(self, target):
        startTime = datetime.datetime.now()

        def fnGetFitness(genes):
            return get_fitness(genes, target)

        def fnDisplay(candidate, mutation):
            display(candidate, mutation, startTime)

        def fnShowIon(genes, target, mutation_attempts):
            show_ion(genes, target, mutation_attempts)

        optimalFitness = 0.99
        best = genetic.get_best(fnGetFitness,\
			optimalFitness, self.geneSet, fnDisplay,\
                        fnShowIon, target)
    
def display(candidate, mutation, startTime):
    timeDiff = datetime.datetime.now() - startTime
    print("{}\t{}\t{}\t{}".format(
	candidate.Genes, candidate.Fitness, mutation, timeDiff))
    
def get_fitness(genes, target):
    cation = Chem.MolFromSmiles(genes)
    anion = Chem.MolFromSmiles("[B-](F)(F)(F)F")
    model = genetic.load_data("density_nn_model.sav", pickleFile=True)
    deslist = genetic.load_data("density_nn_model_descriptors.csv")
    feature_vector=[]
    with genetic.suppress_stdout_stderr():
        for item in deslist:
            if "anion" in item:
                feature_vector.append(calculator([item.partition('-')\
                    [0]]).CalcDescriptors(anion)[0])
            elif "cation" in item:
                feature_vector.append(calculator([item.partition('-')\
                [0]]).CalcDescriptors(cation)[0])          
            elif "Temperature_K" in item:
                feature_vector.append(298.15)
            elif "Pressure_kPa" in item:
                feature_vector.append(101.325)
            else:
                print("unknown descriptor in list: %s" % item)
    features_normalized = (feature_vector-deslist.iloc[0].values)/deslist.iloc[1].values
    prediction = exp(model.predict(np.array(features_normalized).reshape(1,-1))[0])
    error = abs((prediction-target)/target)
    return 1-error

def show_ion(genes, target,  mutation_attempts):
    mol = Chem.MolFromSmiles(genes)
    print("{}\t{}".format("number of atoms: ", mol.GetNumAtoms()))
    print("{}\t{}".format("mutation attempts: ", mutation_attempts))
    print("within 1%% of target density: %s (kg/m) " % target)



if __name__ == '__main__':
    unittest.main()
