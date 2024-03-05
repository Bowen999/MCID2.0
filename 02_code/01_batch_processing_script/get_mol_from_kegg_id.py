# from bioservices import KEGG
# import time
# import os
# import pandas as pd
# from rdkit import Chem
# from tqdm import tqdm

# kegg_db = pd.read_csv('0Rxn_11157.csv')
# compound_ids = kegg_db['mcid'].to_list()


# def kegg_id_to_mol(kegg_id, kegg_service):
#     try:
#         mol_data = kegg_service.get(kegg_id, "mol")
#         return mol_data
#     except Exception as e:
#         print(f"Error fetching MOL for {kegg_id}: {e}")
#         return None

# # Initialize KEGG service
# k = KEGG()

# comp_info_dir = os.path.join(os.getcwd(), 'comp_info')
# os.makedirs(comp_info_dir, exist_ok=True)

# # Use tqdm to show progress bar
# for cid in tqdm(compound_ids, desc="Processing compounds", unit="compound"):
#     mol_data = kegg_id_to_mol(cid, k)
#     if mol_data:
#         # Create a directory for the compound within comp_info
#         compound_dir = os.path.join(comp_info_dir, cid)
#         os.makedirs(compound_dir, exist_ok=True)
        
#         # Save the MOL file within the directory
#         mol_path = os.path.join(compound_dir, f"{cid}.mol")
#         with open(mol_path, "w") as f:
#             f.write(mol_data)
#     time.sleep(1)  # Introduce a delay of 1 second between requests

from bioservices import KEGG
import time
import os
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

kegg_db = pd.read_csv('0Rxn_11157.csv')
compound_ids = kegg_db['mcid'].to_list()

def kegg_id_to_mol(kegg_id, kegg_service):
    try:
        mol_data = kegg_service.get(kegg_id, "mol")
        return mol_data
    except Exception as e:
        print(f"Error fetching MOL for {kegg_id}: {e}")
        return None

# Initialize KEGG service
k = KEGG()

comp_info_dir = os.path.join(os.getcwd(), 'comp_info')
os.makedirs(comp_info_dir, exist_ok=True)

# Remove compounds that already have directories in comp_info
existing_compounds = set(os.listdir(comp_info_dir))
compound_ids = [cid for cid in compound_ids if cid not in existing_compounds]

# Use tqdm to show progress bar
for cid in tqdm(compound_ids, desc="Processing compounds", unit="compound"):
    try:
        mol_data = kegg_id_to_mol(cid, k)
        if mol_data:
            # Create a directory for the compound within comp_info
            compound_dir = os.path.join(comp_info_dir, cid)
            os.makedirs(compound_dir, exist_ok=True)
        
            # Save the MOL file within the directory
            mol_path = os.path.join(compound_dir, f"{cid}.mol")
            with open(mol_path, "w") as f:
                f.write(mol_data)
        time.sleep(1)  # Introduce a delay of 1 second between requests
    except Exception as e:
        print(f"Error processing {cid}: {e}")
