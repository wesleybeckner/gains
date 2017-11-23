from __future__ import absolute_import, division, print_function
from os.path import dirname, join
import numpy as np
import pandas as pd
import scipy.optimize as opt
from scipy.special import erf
#from .due import due, Doi
import subprocess
import rdkit
from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import ShowMol
import pickle
import os
import statistics
import time
import random
import sys
from contextlib import contextmanager

__all__ = ["suppress_stdout_stderr", "Benchmark", "GeneSet", "Chromosome",\
        "generate_geneset", "_generate_parent", "_mutate", "get_best",\
        "load_data"]

# Use duecredit (duecredit.org) to provide a citation to relevant work to
# be cited. This does nothing, unless the user has duecredit installed,
# And calls this with duecredit (as in `python -m duecredit script.py`):
#due.cite(Doi("10.1167/13.9.30"),
#         description="Template project for small scientific Python projects",
#         tags=["reference-implementation"],
#         path='gains')

"""
This GA uses RDKit to make atomic mutations to a starting imidazole.
"""

def load_data(data_file_name, pickleFile=False, simpleList=False):
    """Loads data from module_path/data/data_file_name.
    Parameters
    ----------
    data_file_name : String. Name of csv file to be loaded from
    module_path/data/data_file_name. For example 'salt_info.csv'.
    Returns
    -------
    data : Pandas DataFrame
        A data frame. For example with each row representing one 
        salt and each column representing the features of a given 
        salt.
    """
    module_path = dirname(__file__)
    if pickleFile:
        with open(join(module_path, 'data', data_file_name), 'rb') as pickle_file:
            data = pickle.load(pickle_file, encoding='latin1')
    elif simpleList:
        with open(join(module_path, 'data', data_file_name)) as csv_file:
            data = csv_file.read().splitlines()
    else:
        with open(join(module_path, 'data', data_file_name)) as csv_file:
            data = pd.read_csv(csv_file)
    return data

class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in 
    Python, i.e. will suppress all print, even if the print originates in a 
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).      

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)

class Benchmark:
    @staticmethod
    def run(function):
        timings = []
        stdout = sys.stdout
        for i in range(5):
            sys.stdout = None
            startTime = time.time()
            function()
            seconds = time.time() - startTime
            sys.stdout = stdout
            timings.append(seconds)
            mean = statistics.mean(timings)
            print("{} {:3.2f} {:3.2f}".format(
                1 + i, mean,
                statistics.stdev(timings, mean) if i > 1 else 0))

class GeneSet():
    def __init__(self, atoms, rdkitFrags, customFrags):
        self.Atoms = atoms
        self.RdkitFrags = rdkitFrags
        self.CustomFrags = customFrags

class Chromosome(Chem.rdchem.Mol):
    def __init__(self, genes, fitness):
        Chem.rdchem.Mol.__init__(self)
        self.Genes = genes
        self.Fitness = fitness
        self.Mol = Chem.MolFromSmiles(genes)
        self.RWMol = Chem.MolFromSmiles(genes)
        self.RWMol = Chem.RWMol(Chem.MolFromSmiles(genes))

def generate_geneset():
    atoms = [6,7]
    fName = os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
    rdkitFrags = FragmentCatalog.FragCatParams(1,5,fName)
    customFrags = FragmentCatalog.FragCatalog(rdkitFrags)
    fcgen = FragmentCatalog.FragCatGenerator()
    m = Chem.MolFromSmiles('CCCC')
    fcgen.AddFragsFromMol(m,customFrags)
    return GeneSet(atoms, rdkitFrags, customFrags)

def _generate_parent(geneSet, get_fitness):
    df = load_data("saltInfo.csv")
    df = df.loc[df["cation_name"].str.contains("imid", case=False)]
    df = df['cation_SMILES'].unique()
    ohPickMe = random.sample(range(df.shape[0]),1)
    genes = df[ohPickMe[0]]
    fitness = get_fitness(genes)
    print(fitness)
    return Chromosome(genes, fitness)

def _mutate(parent, geneSet, get_fitness, target):
    def replace_atom(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.Atoms
        if childGenes.RWMol.GetAtomWithIdx(oldGene).IsInRing() == True:
            genes = Chem.MolToSmiles(parent.Mol)
            return Chromosome(genes, 0)
        newGene = random.sample(geneSet, 1)[0]
        childGenes.RWMol.GetAtomWithIdx(oldGene).SetAtomicNum(newGene) 
        return childGenes  
    def add_atom(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.Atoms
        newGeneNumber = childGenes.RWMol.GetNumAtoms()  
        newGene = random.sample(geneSet, 1)[0]
        childGenes.RWMol.AddAtom(Chem.Atom(newGene))
        childGenes.RWMol.AddBond(newGeneNumber,oldGene,Chem.BondType.SINGLE) 
        return childGenes
    def remove_atom(childGenes, GeneSet, oldGene):
        if childGenes.RWMol.GetAtomWithIdx(oldGene).GetExplicitValence() != 1:
            genes = Chem.MolToSmiles(parent.Mol)
            return Chromosome(genes, 0)
        childGenes.RWMol.RemoveAtom(oldGene)
        return childGenes
    def add_custom_fragment(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.CustomFrags
        newGene = Chromosome(geneSet.GetEntryDescription(\
            random.sample(range(geneSet.GetNumEntries()), 1)[0]),0)
        oldGene = oldGene + newGene.Mol.GetNumAtoms()
        combined = Chem.EditableMol(Chem.CombineMols(newGene.Mol,childGenes.Mol))
        combined.AddBond(0,oldGene,order=Chem.rdchem.BondType.SINGLE)
        childGenes = combined.GetMol()   
        try:
            childGenes = Chromosome(Chem.MolToSmiles(childGenes),0)  
            return childGenes
        except:
            return 0
    def add_rdkit_fragment(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.RdkitFrags
        try:
            newGene = Chromosome(Chem.MolToSmiles(geneSet.GetFuncGroup(\
                random.sample(range(geneSet.GetNumFuncGroups()), 1)[0])),0)
        except:
            return 0
        oldGene = oldGene + newGene.Mol.GetNumAtoms() 
        combined = Chem.EditableMol(Chem.CombineMols(newGene.Mol,childGenes.Mol))
        combined.AddBond(1,oldGene,order=Chem.rdchem.BondType.SINGLE)
        combined.RemoveAtom(0)
        try:
            childGenes = Chromosome(Chem.MolToSmiles(combined.GetMol()),0)
            return childGenes
        except:
            return 0
    def remove_custom_fragment(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.CustomFrags
        newGene = Chromosome(geneSet.GetEntryDescription(\
            random.sample(range(geneSet.GetNumEntries()), 1)[0]),0)
        try:
            truncate = Chem.DeleteSubstructs(childGenes.Mol,newGene.Mol)
            childGenes = truncate
            childGenes = Chromosome(Chem.MolToSmiles(childGenes),0)  
            return childGenes
        except:
            return 0
    def remove_rdkit_fragment(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.RdkitFrags
        try:
            newGene = Chromosome(Chem.MolToSmiles(geneSet.GetFuncGroup(\
                random.sample(range(geneSet.GetNumFuncGroups()), 1)[0])),0)
        except:
            return 0
        try:
            truncate = Chem.DeleteSubstructs(childGenes.Mol,newGene.Mol)
            childGenes = truncate
            childGenes = Chromosome(Chem.MolToSmiles(childGenes),0)  
            return childGenes
        except:
            return 0
    childGenes = Chromosome(parent.Genes,0)
    oldGene = random.sample(range(childGenes.RWMol.GetNumAtoms()), 1)[0]
    mutate_operations = [add_atom, remove_atom, remove_custom_fragment,\
	replace_atom, add_rdkit_fragment, add_custom_fragment, remove_rdkit_fragment]
    i = random.choice(range(len(mutate_operations)))
    mutation = mutate_operations[i].__name__
    childGenes = mutate_operations[i](childGenes, geneSet, oldGene)
    try:
        childGenes.RWMol.UpdatePropertyCache(strict=True)
        Chem.SanitizeMol(childGenes.RWMol)
        genes = Chem.MolToSmiles(childGenes.RWMol)
        fitness = get_fitness(genes)
        return Chromosome(genes, fitness), mutation
    except:
        return Chromosome(parent.Genes, 0), mutation
        

def get_best(get_fitness, optimalFitness, geneSet, display,\
        show_ion, target):
    mutation_attempts = 0
    attempts_since_last_adoption = 0
    random.seed()
    bestParent = _generate_parent(geneSet, get_fitness)
    display(bestParent, "starting structure")
    if bestParent.Fitness >= optimalFitness:
        return bestParent
    while True:
        with suppress_stdout_stderr():
            child, mutation = _mutate(bestParent, geneSet, get_fitness, target)
        mutation_attempts += 1
        attempts_since_last_adoption += 1
        if attempts_since_last_adoption > 500:
            bestParent = _generate_parent(geneSet, get_fitness)
            attempts_since_last_adoption = 0
            print("starting from new parent")
        if bestParent.Fitness >= child.Fitness:
            continue
        display(child, mutation)
        attempts_since_last_adoption = 0
        if child.Fitness >= optimalFitness:
            show_ion(child.Genes, target, mutation_attempts)
            return child
        bestParent = child
