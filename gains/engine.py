from __future__ import absolute_import, division, print_function
from os.path import dirname, join
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull
from sklearn.preprocessing import MinMaxScaler
from keras.models import load_model
from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from math import sqrt
import os
import statistics
import time
import random
import sys

__all__ = ["get_best", "molecular_similarity", "suppress_rdkit_sanity",
           "generate_geneset", "load_data", "Chromosome",
           "GeneSet", "Benchmark", "convex_search"]


"""
This GA uses RDKit to search molecular structure
"""


def get_best(get_fitness, optimalFitness, geneSet, display,
             show_ion, target, parent_candidates, seed=None,
             hull=None, simplex=None, verbose=0, hull_bounds=[0, 1],
             inner_search=True, parent_cap=25, mutation_cap=1000):
    """
    the primary public function of the engine

    Parameters
    ----------
    get_fitness : function
        the fitness function. Usually based on a molecular property.
        An example can be found in the salt_generator module
    optimalFitness : float
        0-1 the user specifies how close the engine should get to
        the target (1 = exact)
    geneSet : object
        consists of atomtypes (by periodic number), rdkit molecular
        fragments and custom fragments (that are currently hard
        coded into the engine). These are the building blocks that
        the engine can use to mutate the molecular candidate via the
        _mutate() function
    display : function
        for printing results to the screen. Display is called for
        every accepted mutation
    show_ion : function
        for printing results to the screen. show_ion is called when
        a candidate has achieved the desired fitness score and is
        returned by the engine
    target : array, float, or int
        the desired property value to be achieved by the engine.
        If an array, a model containing multi-output targets must
        be supplied to the engine
    parent_candidates : array
        an array of smiles strings that the engine uses to choose
        a starting atomic configuration
    seed : int, optional
        optional randint seed for unittest consistency
    hull : pandas DataFrame, optional
        nxm pandas DataFrame to use convex hull search strategy. hull
        columns should be the same properties used in the genetic algorithm
        fitness test
    simplex : array, optional
        array to access boundary datapoints in the convex hull. This is used
        during target resampling defined by the convex hull/simplex
    verbose : int, optional, default 0
        0 : most verbose. Best child, parent/target resampling,
            sanitization failure
        1 : parent/target resampling, solution metadata, sanitization failure
        2 : solution metdata, sanitization failure
        3 : target resampling, csv-formatted solution metadata
        4 : csv-formatted solution metadata
    hull_bounds : array, optional
        if hull and simplex are not none, hull_bounds describes the
        proximity convex_search should be to the simplex
    inner_search : bool, optional
        if hull and simplex are not none, inner_search specifies if
        convex_search should return values only within the convex hull

    Returns
    ----------
    child : Chromosome object
        the accepted molecular configuration. See Chromosome class
        for details
    """
    mutation_attempts = 0
    attempts_since_last_adoption = 0
    parents_sampled = 0
    targets_sampled = 0
    if seed:
        random.seed(seed)
    bestParent = _generate_parent(parent_candidates, get_fitness)
    if verbose == 0:
        display(bestParent, "starting structure")
    if bestParent.Fitness >= optimalFitness:
        return bestParent
    while True:
        with suppress_rdkit_sanity():
            child, mutation = _mutate(bestParent, geneSet, get_fitness, target)
        mutation_attempts += 1
        attempts_since_last_adoption += 1

        if attempts_since_last_adoption > mutation_cap:
            if hull is not None and parents_sampled > parent_cap:
                target = convex_search(hull, simplex, hull_bounds,
                                       inner_search)
                targets_sampled += 1
                if verbose in [0, 1, 3]:
                    print("assigning new target: {}".format(target))
            child = _generate_parent(parent_candidates, get_fitness)
            attempts_since_last_adoption = 0
            parents_sampled += 1
            if verbose in [0, 1]:
                print("starting from new parent")
        elif bestParent.Fitness >= child.Fitness:
            continue
        if verbose == 0:
            display(child, mutation)
        attempts_since_last_adoption = 0
        if child.Fitness >= optimalFitness:
            sim_score, sim_index = molecular_similarity(child,
                                                        parent_candidates)
            molecular_relative = parent_candidates[sim_index]
            if verbose in [0, 1, 2]:
                show_ion(child.Genes, target, mutation_attempts, sim_score,
                         molecular_relative)
            return child
        bestParent = child


def convex_search(data, simplex_id, constraints=[0, 1], inner=True):
    """
    need to make sure the minimum distance applies to all simplexes!!!
    then lets wrap this puppy bad boy up and compute some targets
    add to the genetic engine (replace/update convex_search)
    """

    def _grab_candidate(data, hull, simplex, constraints, inner=inner):
        """
        randomly selects candidate and checks w/in min/max of user defined
        boundaries from the hull
        """
        if str(simplex) != 'any':
            x1, x2 = data[simplex, 0]  # simplex is the
            y1, y2 = data[
                simplex, 1]  # index of the datapoint forming the hull
            segment_length = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            m = (y2 - y1) / (x2 - x1)  # slope
            b = y2 - (x2 * m)  # intercept
        while True:
            x3 = random.random()
            y3 = random.random()
            new_hull = ConvexHull(
                np.append(np.array([[x3, y3]]), data, axis=0))
            if str(simplex) != 'any':
                distance = float(abs(m * x3 - y3 + b)) / float(
                    sqrt(m ** 2 + 1))
                point_1_distance = sqrt((x3 - x1) ** 2 + (y3 - y1) ** 2)
                point_2_distance = sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)
                if segment_length > point_1_distance and segment_length >\
                        point_2_distance:

                    if inner and hull.area >= new_hull.area:

                        if distance < constraints[1]:
                            if distance > \
                                    constraints[0]:
                                return x3, y3
                    elif not inner and new_hull.area / hull.area - 1 < \
                            constraints[1]:
                        if distance < constraints[1]:
                            if distance > \
                                    constraints[0]:
                                return x3, y3
            elif inner and hull.area == new_hull.area:  # search inside hull
                return x3, y3

    def _check_minimum_distances(data, hull, x3, y3, constraints, inner=inner):
        """
        checks to see whether the candidate data point is more than
        the minimum distance from any of the boundaries
        """
        distances = []
        for simplex_all in hull.simplices:
            x1_a, x2_a = data[simplex_all, 0]
            y1_a, y2_a = data[simplex_all, 1]
            m_a = (y2_a - y1_a) / (x2_a - x1_a)  # slope
            b_a = y2_a - (x2_a * m_a)  # intercept
            distances.append(
                float(abs(m_a * x3 - y3 + b_a)) / float(sqrt(m_a ** 2 + 1)))
        if any(i < constraints[0] for i in
               distances):  # check minimum constraint
            return True
        else:
            return None

    instance = MinMaxScaler()
    data = instance.fit_transform(data)
    hull = ConvexHull(data)
    if str(simplex_id) != 'any':
        simplex = hull.simplices[simplex_id]
    else:
        simplex = 'any'
    while True:
        x3, y3 = _grab_candidate(data, hull, simplex, constraints, inner=inner)
        check_pass = _check_minimum_distances(data, hull, x3, y3, constraints)
        if not check_pass:
            break
    values = instance.inverse_transform([[x3, y3]])
    return values[0][1], values[0][0]


def molecular_similarity(best, parent_candidates, all=False):
    """
    returns a similarity score (0-1) of best with the
    closest molecular relative in parent_candidates

    Parameters
    ----------
    best : object
        Chromosome object, the current
        mutated candidate
    parent_candidates : array
        parent pool of molecules to compare with best.
        These are represented by SMILES
    all : boolean, optional, default = False
        default behavior is false and the tanimoto
        similarity score is returned. If True
        tanimoto, dice, cosine, sokal, kulczynski,
        and mcconnaughey similarities are returned

    Returns
    ----------
    similarity_score : float
    similarity_index : int
        if all=False the best tanimoto similarity score
        as well as the index of the closest molecular
        relative are returned
        if all=True an array of best scores and indeces
        of the closest molecular relative are returned
    """
    scores = []
    if all:
        indices = []
        metrics = [DataStructs.TanimotoSimilarity,
                   DataStructs.DiceSimilarity,
                   DataStructs.CosineSimilarity,
                   DataStructs.SokalSimilarity,
                   DataStructs.KulczynskiSimilarity,
                   DataStructs.McConnaugheySimilarity]

        for j in range(len(metrics)):

            scores_micro = []
            for i in range(len(parent_candidates)):
                ms = [best.Mol, Chem.MolFromSmiles(parent_candidates[i])]
                fps = [FingerprintMols.FingerprintMol(x) for x in ms]
                score = DataStructs.FingerprintSimilarity(fps[0], fps[1],
                                                          metric=metrics[j])
                scores_micro.append(score)
            scores.append(max(scores_micro))
            indices.append(scores_micro.index(max(scores_micro)))
        return scores, indices
    else:
        for i in range(len(parent_candidates)):
            ms = [best.Mol, Chem.MolFromSmiles(parent_candidates[i])]
            fps = [FingerprintMols.FingerprintMol(x) for x in ms]
            score = DataStructs.FingerprintSimilarity(fps[0], fps[1])
            scores.append(score)
        return max(scores), scores.index(max(scores))


def load_data(data_file_name, h5File=False):
    """
    Loads data from module_path/data/data_file_name.

    Parameters
    ----------
    data_file_name : string
        name of csv file to be loaded from module_path/data/
        data_file_name.
    h5File : boolean, optional, default = False
        if True opens hdf5 file

    Returns
    -------
    data : Pandas DataFrame
    """
    module_path = dirname(__file__)
    if h5File:
        data = load_model(join(module_path, 'data', data_file_name))
    else:
        with open(join(module_path, 'data', data_file_name), 'rb') as csv_file:
            data = pd.read_csv(csv_file, encoding='latin1')
    return data


class suppress_rdkit_sanity(object):
    """
    Context manager for doing a "deep suppression" of stdout and stderr
    during certain calls to RDKit.
    """
    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


class Benchmark:
    """
    benchmark method used by the unittests
    """
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
    """
    Consists of atomtypes (by periodic number), rdkit molecular
    fragments and custom fragments (that are currently hard
    coded into the engine). These are the building blocks that
    the engine can use to mutate the molecular candidate via the
    _mutate() function
    """
    def __init__(self, atoms, rdkitFrags, customFrags):
        self.Atoms = atoms
        self.RdkitFrags = rdkitFrags
        self.CustomFrags = customFrags


class Chromosome(Chem.rdchem.Mol):
    """
    The main object handled by the engine. The Chromosome object
    inherits the RWMol and Mol attributes from rdkit. Two additional
    attributes are added: genes and fitness. Genes is the SMILES
    encoding of the molecule, fitness is the score (0-1) returned
    by the fitness function
    """
    def __init__(self, genes, fitness):
        Chem.rdchem.Mol.__init__(self)
        self.Genes = genes
        self.Fitness = fitness
        self.Mol = Chem.MolFromSmiles(genes)
        self.RWMol = Chem.MolFromSmiles(genes)
        self.RWMol = Chem.RWMol(Chem.MolFromSmiles(genes))


def generate_geneset():
    """
    Populates the GeneSet class with atoms and fragments to be used
    by the engine. As it stands these are hardcoded into the engine
    but will probably be adapted in future versions

    Parameters
    ----------
    None

    Returns
    ----------
    GeneSet : object
        returns an instance of the GeneSet class containing atoms,
        rdkit fragments, and custom fragments
    """
    atoms = [6, 7, 8, 9, 5, 15, 16, 17]
    fName = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
    rdkitFrags = FragmentCatalog.FragCatParams(1, 5, fName)
    customFrags = FragmentCatalog.FragCatalog(rdkitFrags)
    fcgen = FragmentCatalog.FragCatGenerator()
    m = Chem.MolFromSmiles('CCCC')
    fcgen.AddFragsFromMol(m, customFrags)
    return GeneSet(atoms, rdkitFrags, customFrags)


def _generate_parent(parent_candidates, get_fitness):
    df = parent_candidates
    ohPickMe = random.sample(range(df.shape[0]), 1)
    genes = df[ohPickMe[0]]
    fitness, prediction = get_fitness(genes)
    return Chromosome(genes, fitness)


def _mutate(parent, geneSet, get_fitness, target):
    def replace_atom(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.Atoms
        if childGenes.RWMol.GetAtomWithIdx(oldGene).IsInRing():
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
        childGenes.RWMol.AddBond(newGeneNumber, oldGene, Chem.BondType.SINGLE)
        return childGenes

    def remove_atom(childGenes, GeneSet, oldGene):
        if childGenes.RWMol.GetAtomWithIdx(oldGene).GetExplicitValence() != 1:
            genes = Chem.MolToSmiles(parent.Mol)
            return Chromosome(genes, 0)
        childGenes.RWMol.RemoveAtom(oldGene)
        return childGenes

    def add_custom_fragment(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.CustomFrags
        newGene = Chromosome(geneSet.GetEntryDescription(
            random.sample(range(geneSet.GetNumEntries()), 1)[0]), 0)
        oldGene = oldGene + newGene.Mol.GetNumAtoms()
        combined = Chem.EditableMol(Chem.CombineMols(newGene.Mol,
                                    childGenes.Mol))
        combined.AddBond(0, oldGene, order=Chem.rdchem.BondType.SINGLE)
        childGenes = combined.GetMol()
        try:
            childGenes = Chromosome(Chem.MolToSmiles(childGenes), 0)
            return childGenes
        except BaseException:
            return 0

    def add_rdkit_fragment(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.RdkitFrags
        try:
            newGene = Chromosome(Chem.MolToSmiles(geneSet.GetFuncGroup(
                random.sample(range(geneSet.GetNumFuncGroups()), 1)[0])), 0)
        except BaseException:
            return 0
        oldGene = oldGene + newGene.Mol.GetNumAtoms()
        combined = Chem.EditableMol(Chem.CombineMols(newGene.Mol,
                                    childGenes.Mol))
        combined.AddBond(1, oldGene, order=Chem.rdchem.BondType.SINGLE)
        combined.RemoveAtom(0)
        try:
            childGenes = Chromosome(Chem.MolToSmiles(combined.GetMol()), 0)
            return childGenes
        except BaseException:
            return 0

    def remove_custom_fragment(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.CustomFrags
        newGene = Chromosome(geneSet.GetEntryDescription(
            random.sample(range(geneSet.GetNumEntries()), 1)[0]), 0)
        try:
            truncate = Chem.DeleteSubstructs(childGenes.Mol, newGene.Mol)
            childGenes = truncate
            childGenes = Chromosome(Chem.MolToSmiles(childGenes), 0)
            return childGenes
        except BaseException:
            return 0

    def remove_rdkit_fragment(childGenes, GeneSet, oldGene):
        geneSet = GeneSet.RdkitFrags
        try:
            newGene = Chromosome(Chem.MolToSmiles(geneSet.GetFuncGroup(
                random.sample(range(geneSet.GetNumFuncGroups()), 1)[0])), 0)
        except BaseException:
            return 0
        try:
            truncate = Chem.DeleteSubstructs(childGenes.Mol, newGene.Mol)
            childGenes = truncate
            childGenes = Chromosome(Chem.MolToSmiles(childGenes), 0)
            return childGenes
        except BaseException:
            return 0
    childGenes = Chromosome(parent.Genes, 0)
    mutate_operations = [add_atom, remove_atom, remove_custom_fragment,
                         replace_atom, add_rdkit_fragment, add_custom_fragment,
                         remove_rdkit_fragment]
    i = random.choice(range(len(mutate_operations)))
    mutation = mutate_operations[i].__name__
    try:
        oldGene = random.sample(range(childGenes.RWMol.GetNumAtoms()), 1)[0]
    except BaseException:
        return Chromosome(parent.Genes, 0), mutation
    childGenes = mutate_operations[i](childGenes, geneSet, oldGene)
    try:
        childGenes.RWMol.UpdatePropertyCache(strict=True)
        Chem.SanitizeMol(childGenes.RWMol)
        genes = Chem.MolToSmiles(childGenes.RWMol)
        if "." in genes:
            return Chromosome(parent.Genes, 0), mutation
        fitness, prediction = get_fitness(genes)
        return Chromosome(genes, fitness), mutation
    except BaseException:
        return Chromosome(parent.Genes, 0), mutation
