import numpy as np
from os.path import dirname, join
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
import re
import salty
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator as Calculator
from sklearn.preprocessing import StandardScaler
from math import log


def build_model_from_md(df, property_to_model, temperature=[297, 316], pressure=[99, 102],
                        output_ranges=[[200, 3000]], md_temperature=298.15,
                        md_pressure=101.325):
    """
    creates new qspr models using md data

    Parameters
    ----------
    df : pandas DataFrame
        salt_log data from the genetic algorithm. Contains
        the headers 'Salt Smiles' and 'MD Calculation'. Current
        support is only for cpt and density
    property_to_model : str
        current support is for 'cpt' or 'density'
    temperature : array, optional
        temperature bounds on experimental data to add. Default
        297, 316 K
    pressure : array, optional
        pressure bounds on experimental data to add. Default
        99, 102 kpa
    output_ranges : array, optional
        property bounds on experimental data to add. Default
        200, 3000 (kg/m3 or kj/molK)
    md_temperature : float, optional
        temperature used to generate the md data. Default
        298.15 K
    md_pressure : float, optional
        pressure used to generate the md data. Dfault
        101.325 kPa

    Returns
    -------
    newmodel : salt dev_model object
    new_MD_data_index : int
        start index of the newly incorporated MD data

    Summary
    -------
    Create 4 lists from df: cation/anion smiles, cpt, density
    Nans will be used for cation/anion name in the newmodel
    output
    """

    cpt = []
    density = []
    cation_smi = []
    anion_smi = []
    for i in range(df.shape[0]):
        calculation = df["MD Calculation"][i]
        cpt.append(re.findall("\d+\.\d+", calculation)[0])
        density.append(re.findall("\d+\.\d+", calculation)[1])
        cation_smi.append(df['Salt Smiles'][i].split(".")[0])
        anion_smi.append(df['Salt Smiles'][i].split(".")[1])

    ###based on coco's Deslist
    module_path = dirname(__file__)
    data = df
    n = data.shape[0]
    f = open(join(module_path, 'data', 'Deslist'), 'r')
    Deslist = []
    for line in f:
        Deslist.append(line.strip('\n\t'))
    calc = Calculator(Deslist)
    D = len(Deslist)
    d = len(Deslist) * 2 + 8  # *2 for cat/an desc then T,P,prop1,prop2,name,smi
    X = np.zeros((n, d))
    X[:, -8] = md_temperature
    X[:, -7] = md_pressure
    for i in range(n):
        cation = Chem.MolFromSmiles(cation_smi[i])
        anion = Chem.MolFromSmiles(anion_smi[i])
        X[i][:D] = calc.CalcDescriptors(cation)
        X[i][D:2 * D] = calc.CalcDescriptors(anion)
    X[:, -5] = density
    X[:, -6] = cpt
    cols_cat = [s + "-cation" for s in Deslist]
    cols_ani = [s + "-anion" for s in Deslist]
    cols = cols_cat + cols_ani + ["Temperature, K", "Pressure, kPa",
                                  "Heat capacity at constant pressure, J/K/mol",
                                  "Specific density, kg/m<SUP>3</SUP>", "name-anion",
                                  "smiles-anion", "name-cation", "smiles-cation"]
    X = pd.DataFrame(X, columns=cols)
    X.iloc[:, -4] = np.nan
    X.iloc[:, -2] = np.nan
    X.iloc[:, -3] = anion_smi
    X.iloc[:, -1] = cation_smi  # X is the df with the new simulation data
    new_MD_data_index = X.shape[0]  # this will be used to plot the new data predictions after model re-training
    devmodel = salty.aggregate_data([property_to_model], T=temperature, P=pressure,
                                    data_ranges=output_ranges, scale_center=False)
    cols = devmodel.Data.columns
    new_data = pd.concat([devmodel.Data, X])

    ###Property strings
    # Specific density, kg/m<SUP>3</SUP>
    # Heat capacity at constant pressure, J/K/mol
    if property_to_model == 'density':
        prop = "Specific density, kg/m<SUP>3</SUP>"
        to_drop = "Heat capacity at constant pressure, J/K/mol"
    else:
        to_drop = "Specific density, kg/m<SUP>3</SUP>"
        prop = "Heat capacity at constant pressure, J/K/mol"

    new_data.drop(columns=[to_drop], inplace=True)
    new_data = new_data[cols]
    new_data.reset_index(inplace=True, drop=True)

    # the following is modified
    # stuff from aggregate_data
    # Create summary of dataset
    exp_data = [prop, "Temperature, K", "Pressure, kPa"]
    merged = new_data
    unique_salts = merged["smiles-cation"] + merged["smiles-anion"]
    unique_cations = repr(merged["smiles-cation"].unique())
    unique_anions = repr(merged["smiles-anion"].unique())
    actual_data_ranges = []
    for i in range(len(exp_data)):
        actual_data_ranges.append("{} - {}".format(
            str(merged[exp_data[i]].min()), str(merged[exp_data[i]].max())))
    a = np.array([len(unique_salts.unique()), unique_cations, unique_anions,
                  len(unique_salts)])
    a = np.concatenate((a, actual_data_ranges))
    cols1 = ["Unique salts", "Cations", "Anions", "Total datapoints"]
    cols = cols1 + exp_data
    data_summary = pd.DataFrame(a, cols)
    merged = new_data
    metaDf = merged.select_dtypes(include=["object"])
    dataDf = merged.select_dtypes(include=[np.number])
    cols = dataDf.columns.tolist()
    instance = StandardScaler()
    dataDf.is_copy = False
    dataDf.iloc[:, -1] = dataDf.iloc[:, -1].apply(lambda x: log(float(x)))
    scaled_data = pd.DataFrame(instance.fit_transform(
        dataDf.iloc[:, :-1]), columns=cols[:-1])
    df = pd.concat([scaled_data, dataDf.iloc[:, -1:], metaDf], axis=1)
    mean_std_of_coeffs = pd.DataFrame([instance.mean_, instance.scale_],
                                      columns=cols[:-1])
    new_model = salty.dev_model(mean_std_of_coeffs, data_summary, df)
    print(new_model.Data_summary)
    return new_model, new_MD_data_index
