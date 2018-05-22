model_ID = "cpt_density_m1"
target = [627,1390]
pdb_files = 9
write_file = False

def get_fitness(anion, genes, target, model_ID):
    cation = Chem.MolFromSmiles(genes)
    model = genetic.load_data("{}.sav".format(model_ID), pickleFile=True)
    deslist = genetic.load_data("{}_descriptors.csv".format(model_ID))
    feature_vector = []

    for item in deslist:

        if "anion" in item:
            with genetic.suppress_stdout_stderr():
                feature_vector.append(calculator([item.partition('-')
                                      [0]]).CalcDescriptors(anion)[0])
        elif "cation" in item:
            with genetic.suppress_stdout_stderr():
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
    prediction = np.round(np.exp(model.predict(np.array(features_normalized).
                          reshape(1, -1))[0]), decimals=2)
    error = abs((prediction - target) / target)
    error = np.average(error)

    return 1 - error, prediction

summary = genetic.load_data("{}_summ.csv".format(model_ID))
parent_candidates = None
if parent_candidates is None:
    parent_candidates = eval(summary.iloc[1][0])
anion_candidates = eval(summary.iloc[2][0])
cols = ["Salt ID", "Salt Smiles", "Cation Heavy Atoms",
        "Tanimoto Similarity Score", "Molecular Relative", "Anion",
        "Model Prediction", "MD Calculation", "Error"]
salts = pd.DataFrame(columns=cols)

for i in range(1, pdb_files + 1):
    anion_smiles = random.sample(list(anion_candidates), 1)[0]
    if i < 10:
        anion = Chem.MolFromPDBFile("A0{}.pdb".format(i))
        cation = Chem.MolFromPDBFile("C0{}.pdb".format(i))
    else:
        anion = Chem.MolFromPDBFile("A{}.pdb".format(i))
        cation = Chem.MolFromPDBFile("C{}.pdb".format(i))
    best = genetic.Chromosome(Chem.MolToSmiles(cation), 0)
    tan_sim_score, sim_index =\
        genetic.molecular_similarity(best, parent_candidates)
    cation_heavy_atoms = best.Mol.GetNumAtoms()
    salt_smiles = best.Genes + "." + Chem.MolToSmiles(anion)

    scr, pre = get_fitness(anion, best.Genes, target, model_ID)
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
    salts = new

if write_file:
    pd.DataFrame.to_csv(new, path_or_buf="salt_log.csv", index=False)    

