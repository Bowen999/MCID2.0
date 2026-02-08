from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
import pubchempy as pcp

# Load your DataFrame
df = pd.read_csv('mcidxkegg_1811882_all.csv')
df = df.head(50000)

# Define SMARTS patterns for different functional groups
functional_groups_smarts = {
    'amine': '[NX3;H2,H1,H0;!$(NC=O)]',
    'phenol': 'c1cc(ccc1O)O',
    'hydroxyl': '[OX2H]',
    'carboxyl': 'C(=O)[OX1H0-,OX2H1]',
    'carbonyl': '[CX3]=[OX1]'
}

def identify_functional_groups(mol):
    groups_present = {}
    for group_name, smarts in functional_groups_smarts.items():
        substructure = Chem.MolFromSmarts(smarts)
        groups_present[group_name] = mol.HasSubstructMatch(substructure)
    return groups_present

def smiles2info(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            formula = rdMolDescriptors.CalcMolFormula(mol)
            groups = identify_functional_groups(mol)
            compounds = pcp.get_compounds(smiles, namespace='smiles')
            if compounds:
                match = compounds[0]
                common_name = match.synonyms[0] if match.synonyms else None
                iupac_name = match.iupac_name if match.iupac_name else None
            else:
                common_name, iupac_name = None, None
            return [common_name, iupac_name, formula] + list(groups.values())
    except Exception as e:
        print(f"Error with SMILES {smiles}: {e}")
        return [smiles, None, None] + [None] * len(functional_groups_smarts)

# Function to process a chunk of the DataFrame and save to CSV
def process_and_save_chunk(chunk, chunk_idx, columns):
    num_cores = -1  # Use all available cores
    results = Parallel(n_jobs=num_cores)(delayed(smiles2info)(smile) for smile in tqdm(chunk['smiles'], total=chunk.shape[0]))
    
    # Create a new DataFrame for the chunk
    chunk_df = pd.DataFrame(results, columns=columns)
    
    # Add any additional columns from the original DataFrame if needed
    for col in chunk.columns.difference(['smiles']):
        chunk_df[col] = chunk[col].values
    
    # Save the chunk to CSV
    chunk_df.to_csv(f'mcidxkegg_with_functional_groups_chunk_{chunk_idx}.csv', index=False)

# Define the columns for the new DataFrame
columns = ['CommonName', 'IUPACName', 'Formula'] + list(functional_groups_smarts.keys())

# Process the DataFrame in chunks and save each chunk to a CSV file
chunk_size = 10000
for start in range(0, df.shape[0], chunk_size):
    end = min(start + chunk_size, df.shape[0])
    chunk = df[start:end]
    process_and_save_chunk(chunk, start // chunk_size, columns)
