from rdkit import Chem
from standardiser import standardise
import requests
import numpy as np
from sklearn.model_selection import train_test_split


def standardise_smiles(smiles):
    mols = []
    for smi in smiles:
        try:
            mol = Chem.MolFromSmiles(smi)
        except:
            mol=np.nan
        mols += [mol]
    st_mols = []
    for mol in mols:
        if mol is not None:
            try:
                st_mol = standardise.run(mol)
            except:
                st_mol = np.nan
        else:
            st_mol = np.nan
        st_mols += [st_mol]

    st_smiles = []
    for st_mol in st_mols:
        if st_mol is not None:
            try:
                st_smi = Chem.MolToSmiles(st_mol)
            except:
                st_smi=np.nan
        else:
            st_smi = np.nan
        st_smiles += [st_smi]
    return st_smiles

def pubchem_smiles(molname):
    try:
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/CanonicalSMILES/JSON".format(molname)
        response = requests.get(url)
        data = response.json()
        canonical_smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except:
        canonical_smiles = np.nan
    return canonical_smiles

def random_split(df, size):
    indices = np.arange(len(df))
    X_train, X_test, y_train, y_test, i_train, i_test = train_test_split(df["st_smiles"], df["bin"], indices, test_size=size, stratify=df["bin"])
    train = df.iloc[i_train]
    test = df.iloc[i_test]
    return train, test


    