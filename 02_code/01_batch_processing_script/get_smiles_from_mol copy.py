import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
import os

# Initialize tqdm for pandas
tqdm.pandas()

df = pd.read_csv('0Rxn_11157.csv')
comp_info_dir = os.path.join(os.getcwd(), 'comp_info')

def mol_to_smiles(mol_path):
    try:
        mol = Chem.MolFromMolFile(mol_path)
        return AllChem.MolToSmiles(mol)
    except Exception as e:
        print(f"Error converting MOL to SMILES for {mol_path}: {e}")
        return None

# Extract SMILES strings and add as a new column with a progress bar
df["SMILES"] = df["mcid"].progress_apply(lambda cid: mol_to_smiles(os.path.join(comp_info_dir, cid, f"{cid}.mol")))
df.to_csv('0Rxn_11157_with_smiles.csv', index=False)
print(df)