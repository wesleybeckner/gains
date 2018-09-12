from math import sqrt
from scipy.spatial import ConvexHull
from sklearn.preprocessing import MinMaxScaler
from sklearn.neighbors import KernelDensity
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


def build_model_from_md(df, property_to_model, temperature=[298.1, 299], pressure=[101, 102],
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

    devmodel = salty.aggregate_data(property_to_model, T=temperature, P=pressure,
                                    data_ranges=output_ranges, scale_center=False)
    cols = devmodel.Data.columns
    new_data = pd.concat([devmodel.Data, X])  # may have to sort in future version

    if property_to_model == ['density']:
        prop = "Specific density, kg/m<SUP>3</SUP>"
        to_drop = "Heat capacity at constant pressure, J/K/mol"
    elif property_to_model == ['cpt']:
        to_drop = "Specific density, kg/m<SUP>3</SUP>"
        prop = "Heat capacity at constant pressure, J/K/mol"
    elif property_to_model == ["cpt", "density"]:
        prop = ["Heat capacity at constant pressure, J/K/mol", "Specific density, kg/m<SUP>3</SUP>"]

    if property_to_model != ["cpt", "density"]:
        new_data.drop(columns=[to_drop], inplace=True)

    new_data = new_data[cols]
    new_data.reset_index(inplace=True, drop=True)

    if property_to_model == ["cpt", "density"]:
        exp_data = [prop[0], prop[1], "Temperature, K", "Pressure, kPa"]
    else:
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
    for i in range(1, len(property_to_model) + 1):
        dataDf.iloc[:, -i] = dataDf.iloc[:, -i].apply(lambda x: log(float(x)))

    scaled_data = pd.DataFrame(instance.fit_transform(
        dataDf.iloc[:, :-len(property_to_model)]), columns=cols[:-len(property_to_model)])
    df = pd.concat([scaled_data, dataDf.iloc[:, -len(property_to_model):], metaDf],
                   axis=1)  # may have to sort in future version
    mean_std_of_coeffs = pd.DataFrame([instance.mean_, instance.scale_],
                                      columns=cols[:-len(property_to_model)])
    new_model = salty.dev_model(mean_std_of_coeffs, data_summary, df)
    print(new_model.Data_summary)
    return new_model, new_MD_data_index


def calculate_minimum_distances(data, x3, y3):
    """
    calculates the minimum distance of x3,y3 from any
    boundary of the convex hull

    Parameters
    ----------
    data : pandas DataFrame
        2-column DataFrame comprising the convex hull
    x3 : float
        data point associated with the first column
    y3 : float
        data point associated with the second column

    Returns
    -------
    minimum distance : float
        percent distance from the nearest edge of the convex hull
    """
    instance = MinMaxScaler(feature_range=(0.1, 0.9))
    data = instance.fit_transform(data)

    [[x3, y3]] = instance.transform([[x3, y3]])

    hull = ConvexHull(data)
    distances = []

    for simplex_all in hull.simplices:
        x1_a, x2_a = data[simplex_all, 0]
        y1_a, y2_a = data[simplex_all, 1]
        m_a = (y2_a - y1_a) / (x2_a - x1_a)  # slope
        b_a = y2_a - (x2_a * m_a)  # intercept
        distances.append(
            float(abs(m_a * x3 - y3 + b_a)) / float(sqrt(m_a ** 2 + 1)))

    new_hull = ConvexHull(
        np.append(np.array([[x3, y3]]), data, axis=0))
    if hull.area >= new_hull.area:
        return (-np.min(distances))
    else:
        return (np.min(distances))


def gaussian_pdf(column):
    x = column.values
    x_d = np.linspace(min(x), max(x), 10000)

    # instantiate and fit the KDE model
    kde = KernelDensity(bandwidth=0.01, kernel='gaussian')
    kde.fit(x[:, None])

    # score_samples returns the log of the probability density
    return kde.score_samples(x_d[:, None]), x_d
