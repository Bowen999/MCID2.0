import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pubchempy as pcp

df = pd.read_csv('mcidxkegg_1811882_all.csv')
df = df.head(50000)

def smiles2name(smiles):
    """
    Convert a SMILES string to its common name, IUPAC name, or molecular formula.
    If all conversions fail, return the original SMILES string.
    
    :param smiles: A SMILES representation of a molecule.
    :return: A string representing the molecule's common name, IUPAC name, formula, or original SMILES.
    """
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds:
            match = compounds[0]
            
            # Extract common name
            if match.synonyms:
                return match.synonyms[0]
            
            # Extract IUPAC name if common name is absent
            if match.iupac_name:
                return match.iupac_name
        
        # Compute molecular formula if both names are absent
        mol = Chem.MolFromSmiles(smiles)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        if formula:
            return formula
        
    except Exception as e:
        print(f"Error with SMILES {smiles}: {e}")

    # Fallback to original SMILES if all above stages fail
    return smiles

# Parallel computation with progress bar
num_cores = -3  # Use all available cores
results = Parallel(n_jobs=num_cores)(delayed(smiles2name)(smile) for smile in tqdm(df['smiles'], total=df.shape[0]))

df['Name'] = results
df.to_csv('mcidxkegg_with_name.csv', index=False)
