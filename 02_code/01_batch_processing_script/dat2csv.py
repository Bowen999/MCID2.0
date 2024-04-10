# Import necessary libraries
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

# Define a function to count the total lines in the file with specified encoding
def count_lines_with_encoding(file_path, encoding='ISO-8859-1'):
    with open(file_path, 'r', encoding=encoding) as file:
        return sum(1 for _ in file)

def parse_compounds_to_dataframe_with_mass(file_path, encoding='ISO-8859-1'):
    # Define the column names
    column_names = ['UNIQUE-ID', 'TYPES', 'COMMON-NAME', 'INCHI', 'INCHI-KEY', 'SMILES', 'Mass']
    
    # Initialize an empty list to hold the compound data
    compounds_data = []
    
    # Initialize an empty dictionary to hold the current compound's data
    current_compound = {key: 'Na' for key in column_names}
    
    # Get the total number of lines in the file for the progress bar with specified encoding
    total_lines = count_lines_with_encoding(file_path, encoding)
    
    # Open the file and read line by line with progress bar and specified encoding
    with open(file_path, 'r', encoding=encoding) as file:
        for line in tqdm(file, total=total_lines, desc="Processing"):
            # Check for the end of a compound block
            if line.strip() == '//':
                # Calculate Monoisotopic Mass if SMILES is not 'Na'
                smiles = current_compound['SMILES']
                if smiles != 'Na':
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:  # Ensure the molecule could be parsed
                        mass = ExactMolWt(mol)
                    else:
                        mass = 'Na'
                else:
                    mass = 'Na'
                current_compound['Mass'] = mass

                # Append the current compound to the compounds data list
                compounds_data.append(current_compound)
                # Reset the current compound dictionary
                current_compound = {key: 'Na' for key in column_names}
                continue
            
            # Split the line into key and value
            parts = line.split(' - ', 1)
            if len(parts) == 2:
                key, value = parts
                key = key.strip()
                value = value.strip()
                # Check if the key is in the column names and update the current compound dictionary
                if key in column_names:
                    current_compound[key] = value
    
    # Create a DataFrame from the compounds data
    df = pd.DataFrame(compounds_data, columns=column_names)
    
    # Count and print the number of rows where SMILES is not 'Na'
    valid_smiles_count = df[df['SMILES'] != 'Na'].shape[0]
    print(f"Number of rows with valid SMILES: {valid_smiles_count}")
    
    return df

    

# Supporting function for counting lines in the file, needed for progress bar
def count_lines_with_encoding(file_path, encoding='ISO-8859-1'):
    with open(file_path, 'r', encoding=encoding) as file:
        return sum(1 for _ in file)

# Note: Before running the above function, ensure RDKit is installed in your Python environment.



df = parse_compounds_to_dataframe_with_mass('/Users/bowen/Desktop/MCID2.0/01_source_data/microbe/Algal/data/compounds.dat')
df = df[df['SMILES'] != 'Na']
df.to_csv('/Users/bowen/Desktop/Algal.csv')
valid_smiles_count = df[df['SMILES'] != 'Na'].shape[0]
print(f"Number of rows with valid SMILES: {valid_smiles_count}")