import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import multiprocessing
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map  # Import this


tqdm.pandas()

df = pd.read_csv('mcidxkegg_1811882_all.csv')


def smiles_to_iupac_name(smiles):
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds:
            match = compounds[0]
            return match.iupac_name
    except Exception as e:
        print(f"Error with SMILES {smiles}: {e}")
    return smiles


def smiles_to_common_name(smiles):
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds:
            match = compounds[0]
            # `synonyms` is a list of alternate names for the compound; the first one is usually the most common name.
            return match.synonyms[0] if match.synonyms else smiles
    except Exception as e:
        print(f"Error with SMILES {smiles}: {e}")
    return smiles

def smiles_to_formula(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.rdMolDescriptors.CalcMolFormula(mol)
    except:
        return None
    
    
def smiles2name(smiles):
    """
    Convert a SMILES string to its common name, IUPAC name, or molecular formula.
    If all conversions fail, return the original SMILES string.
    
    :param smiles: A SMILES representation of a molecule.
    :return: A string representing the molecule's common name, IUPAC name, formula, or original SMILES.
    """
    
    # [Stage 1: Fetch compound details from PubChem]
    try:
        print(smiles)
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds:
            match = compounds[0]
            
            # [Stage 2: Extract common name]
            if match.synonyms:
                return match.synonyms[0]

            
            # [Stage 3: Extract IUPAC name if common name is absent]
            if match.iupac_name:
                return match.iupac_name
        
        # [Stage 4: Compute molecular formula if both names are absent]
        mol = Chem.MolFromSmiles(smiles)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        if formula:
            return formula
        
    except Exception as e:
        print(f"Error with SMILES {smiles}: {e}")

    # [Stage 5: Fallback to original SMILES if all above stages fail]
    return smiles




# Parallelize the `smiles2name` function across all cores
def process_chunk(chunk):
    return chunk['smiles'].apply(smiles2name)

# def parallel_apply(df, func, n_cores=4):
#     # Split the dataframe into chunks
#     chunks = np.array_split(df, n_cores)
    
#     # Use ProcessPoolExecutor to process chunks in parallel
#     with ProcessPoolExecutor(max_workers=n_cores) as executor:
#         results = list(executor.map(func, chunks))
    
#     # Concatenate results back into single dataframe
#     return pd.concat(results)
def parallel_apply(df, func, n_cores=4):
    # Create tqdm progress bar
    pbar = tqdm(total=len(df), desc="Processing", dynamic_ncols=True)
    
    results = []
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        for result in executor.map(func, df['smiles']):
            results.append(result)
            pbar.update(1)  # Update progress bar by one for each row processed
    
    pbar.close()
    
    # Convert results list to a pandas Series
    return pd.Series(results)



if __name__ == '__main__':
    # Set start method to 'fork' if you're on macOS
    if multiprocessing.get_start_method(allow_none=True) != 'fork' and hasattr(multiprocessing, 'set_start_method'):
        multiprocessing.set_start_method('fork')

    df['Name'] = parallel_apply(df, process_chunk)
    df.to_csv('mcidxkegg_with_name.csv', index=False)