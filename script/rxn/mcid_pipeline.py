from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw

from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdchem import Mol

import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
import pandas as pd
import io
from tqdm import tqdm




"""
 ___                   _
|_ _|_ __  _ __  _   _| |_
 | || '_ \| '_ \| | | | __|
 | || | | | |_) | |_| | |_
|___|_| |_| .__/ \__,_|\__|
          |_|
"""
initial_db = pd.read_parquet('/Users/bowen/Desktop/mcid/database/kegg_compounds_demo_100.parquet')
rxn_rules = pd.read_csv('rxn_info.csv', encoding='latin1')
outputfile = '/Users/bowen/Desktop/mcid/database/kegg_demo_1_rxn.parquet'





"""
 _____                 _   _
|  ___|   _ _ __   ___| |_(_) ___  _ __  ___
| |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
|  _|| |_| | | | | (__| |_| | (_) | | | \__ \
|_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
"""
def predict_products(substrate_smiles, smirks_pattern):
    try:
        # Convert the substrate SMILES to a molecule object
        substrate_mol = Chem.MolFromSmiles(substrate_smiles)

        # Check if conversion was successful; if not, return an empty list
        if substrate_mol is None:
            return []

        # Create a reaction object from the SMIRKS pattern
        reaction = AllChem.ReactionFromSmarts(smirks_pattern)

        # Apply the reaction to the substrate molecule
        products_sets = reaction.RunReactants((substrate_mol,))

        # If no products are generated, return an empty list
        if not products_sets:
            return []

        # Initialize a set to store unique product SMILES
        unique_products_smiles = set()

        # Iterate through the product sets and convert each product to a SMILES string
        for product_set in products_sets:
            for product in product_set:
                # Canonicalize the SMILES and add to the set to ensure uniqueness
                product_smiles = Chem.MolToSmiles(product, isomericSmiles=True)
                unique_products_smiles.add(product_smiles)

        # If no unique products were added, return an empty list
        if not unique_products_smiles:
            return []

        # Convert the set to a list and return
        return list(unique_products_smiles)

    except Exception as e:
        # If any error occurs, return an empty list
        return []




def filter_by_mass_difference(product_smiles, substrate_smiles, mass_difference):
    """
    Filters a list of product molecules based on their monoisotopic mass difference relative to a substrate molecule,
    both represented as SMILES strings.
    Products are filtered to either those with the minimum positive mass difference (mass_difference = '+')
    or those with the maximum negative mass difference (mass_difference = '-') from the substrate.
    Invalid SMILES are ignored.

    Args:
    product_smiles (list of str): List of SMILES strings for product molecules.
    substrate_smiles (str): SMILES string for the substrate molecule.
    mass_difference (str): Specifies the direction of mass difference to filter by ('+' for positive, '-' for negative).

    Returns:
    list of str: SMILES strings of products matching the specified mass difference criteria.

    Raises:
    ValueError: If substrate SMILES is invalid or mass_difference is not '+' or '-'"
    """
    try:
        # Convert substrate SMILES to molecule and calculate its monoisotopic mass
        substrate_mol = Chem.MolFromSmiles(substrate_smiles)
        if substrate_mol is None:
            raise ValueError("Invalid substrate SMILES.")

        substrate_mass = Descriptors.ExactMolWt(substrate_mol)

        # Initialize a list to hold products and their masses
        product_masses = []

        # Iterate through the product SMILES
        for smi in product_smiles:
            # Convert product SMILES to molecule and calculate its monoisotopic mass
            product_mol = Chem.MolFromSmiles(smi)
            if product_mol is None:
                continue  # Skip invalid product SMILES

            product_mass = Descriptors.ExactMolWt(product_mol)
            mass_diff = product_mass - substrate_mass

            # Store product masses and SMILES
            if (mass_difference == '+' and mass_diff > 0) or (mass_difference == '-' and mass_diff < 0):
                product_masses.append((mass_diff, smi))

        # Filter based on mass_difference criteria
        if mass_difference == '+':
            # Find the minimum positive mass difference
            min_mass_diff = min([x[0] for x in product_masses] or [0])
            # Filter products with the minimum mass difference
            matching_products = [smi for diff, smi in product_masses if diff == min_mass_diff]
        elif mass_difference == '-':
            # Find the maximum negative mass difference
            max_mass_diff = max([x[0] for x in product_masses] or [0])
            # Filter products with the maximum mass difference
            matching_products = [smi for diff, smi in product_masses if diff == max_mass_diff]
        else:
            # If mass_difference is not '+' or '-', return an empty list (or handle as needed)
            raise ValueError("Invalid mass_difference value. Use '+' or '-'.")

        return matching_products

    except Exception as e:
        # Handle exceptions, possibly logging them or returning an error message
        print(f"An error occurred: {e}")  # Optionally log the error message
        return []


def find_sub_structure(smiles, smarts):
    """
    This function finds and highlights the substructure in a molecule defined by a SMARTS pattern.

    Parameters:
    - smiles (str): The SMILES representation of the molecule.
    - smarts (str): The SMARTS pattern to search for within the molecule.

    Returns:
    - An image of the molecule with the matching substructure highlighted, if any matches are found.
    """
    # Convert the SMILES string to an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "Invalid SMILES string."

    # Convert the SMARTS string to an RDKit molecule object
    query = Chem.MolFromSmarts(smarts)
    if not query:
        return "Invalid SMARTS pattern."

    # Find the atoms in the molecule that match the SMARTS pattern
    matches = mol.GetSubstructMatches(query, uniquify=False)

    print(matches)
    img = Draw.MolToImage(mol, highlightAtoms=sum(matches, ()), subImgSize=(500, 500))
    return img


def smarts_to_formula(smarts):
    # Convert the SMARTS string to a molecule object
    molecule = Chem.MolFromSmarts(smarts)
    if molecule is None:
        return "Invalid SMARTS string"
    try:
        # Attempt to sanitize the molecule to ensure properties can be calculated
        Chem.SanitizeMol(molecule)
        # Calculate the molecular formula
        formula = rdMolDescriptors.CalcMolFormula(molecule)
        return formula + ' (ignore H)'
    except Exception as e:
        return f"Error processing molecule: {str(e)}"


def plot_reaction_scheme(substrate_smiles, product_smiles, smarts_pattern):
    try:
        # Create RDKit molecule objects
        mol_substrate = Chem.MolFromSmiles(substrate_smiles)
        mol_product = Chem.MolFromSmiles(product_smiles)

        # Generate images for the substrate and product
        substrate_img = Draw.MolToImage(mol_substrate, size=(300, 300))
        product_img = Draw.MolToImage(mol_product, size=(300, 300))

        # Find and generate the sub-structure image
        sub_structure_img = find_sub_structure(mol_substrate, smarts_pattern)
    except Exception as e:
        print(f"Error processing molecules: {e}")
        return

    # Create a figure and axes with 1 row and 5 columns
    fig, ax = plt.subplots(1, 5, figsize=(15, 3))

    # Plotting substrate
    try:
        ax[0].imshow(substrate_img)
    except Exception:
        ax[0].text(0.5, 0.5, 'Missing', horizontalalignment='center', verticalalignment='center')
    ax[0].axis("off")

    # First arrow
    ax[1].text(0.5, 0.5, '→', horizontalalignment='center', verticalalignment='center', fontsize=30, color="darkgrey")
    ax[1].axis("off")

    # Plotting sub-structure
    try:
        ax[2].imshow(sub_structure_img)
    except Exception:
        ax[2].text(0.5, 0.5, 'Missing', horizontalalignment='center', verticalalignment='center')
    ax[2].axis("off")

    # Second arrow
    ax[3].text(0.5, 0.5, '→', horizontalalignment='center', verticalalignment='center', fontsize=30, color="darkgrey")
    ax[3].axis("off")

    # Plotting product
    try:
        ax[4].imshow(product_img)
    except Exception:
        ax[4].text(0.5, 0.5, 'Missing', horizontalalignment='center', verticalalignment='center')
    ax[4].axis("off")

    plt.tight_layout()
    plt.show()


def find_sub_structure(smiles, smarts):
    """
    Finds and highlights the substructure in a molecule defined by a SMARTS pattern.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "Invalid SMILES string."

    query = Chem.MolFromSmarts(smarts)
    if not query:
        return "Invalid SMARTS pattern."

    matches = mol.GetSubstructMatches(query, uniquify=True)
    if not matches:
        return None  # Adjusted to return None for no matches

    img = Draw.MolToImage(mol, highlightAtoms=sum(matches, ()), size=(300, 300))
    return img

def plot_reaction_scheme(substrate_smiles, product_smiles, smarts_pattern):
    """
    Plots a reaction scheme showing the substrate, sub-structure, and product.
    """
    try:
        mol_substrate = Chem.MolFromSmiles(substrate_smiles)
        substrate_img = Draw.MolToImage(mol_substrate, size=(300, 300))
    except Exception as e:
        substrate_img = "Error"

    try:
        mol_product = Chem.MolFromSmiles(product_smiles)
        product_img = Draw.MolToImage(mol_product, size=(300, 300))
    except Exception as e:
        product_img = "Error"

    # Use the provided find_sub_structure function
    sub_structure_img = find_sub_structure(substrate_smiles, smarts_pattern)
    if isinstance(sub_structure_img, str):  # Check if the function returned an error message
        sub_structure_img = "Error"

    # Create a figure and axes with 1 row and 5 columns
    fig, ax = plt.subplots(1, 5, figsize=(15, 3))

    # Helper function to plot images or placeholders
    def plot_image_or_placeholder(ax, img, placeholder_text="Missing"):
        if img == "Error":
            ax.text(0.5, 0.5, placeholder_text, horizontalalignment='center', verticalalignment='center', fontsize=12)
            ax.axis("off")
        else:
            ax.imshow(img)
            ax.axis("off")

    # Plotting substrate
    plot_image_or_placeholder(ax[0], substrate_img)

    # First arrow
    ax[1].text(0.5, 0.5, '→', horizontalalignment='center', verticalalignment='center', fontsize=30, color="royalblue")
    ax[1].axis("off")

    # Plotting sub-structure
    plot_image_or_placeholder(ax[2], sub_structure_img, "No Match")

    # Second arrow
    ax[3].text(0.5, 0.5, '→', horizontalalignment='center', verticalalignment='center', fontsize=30, color="royalblue")
    ax[3].axis("off")

    # Plotting product
    plot_image_or_placeholder(ax[4], product_img)

    plt.tight_layout()
    plt.show()


def predict_products_group(substrate):
    subgroup='[#6]1[#6][#6][#6][#6][#6]1[OH]'
    smirks='[OH:1]>>[O:1]-S(=O)(=O)O'
    mol = Chem.MolFromSmiles(substrate)
    pattern = Chem.MolFromSmarts(subgroup)
    matches = mol.GetSubstructMatches(pattern)
    data = []

    reaction = AllChem.ReactionFromSmarts(smirks)

    for match in matches:
        products = reaction.RunReactants((mol,))
        if products:
            product_mol = products[0][0]
            modified_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=True)
            data.append({
                'substrate_smiles': substrate,
                'smiles': modified_smiles,
                'reaction': 'sulfate conjugation',
                'reaction_rules': smirks,
                'reaction_id': 'R26'
            })

    if not data:  # In case no reaction takes place
        data.append({
            'substrate_smiles': substrate,
            'smiles': substrate,
            'reaction': 'sulfate conjugation',
            'reaction_rules': smirks,
            'reaction_id': 'R26'
        })

    return pd.DataFrame(data)

def drop_group(substrate, replacement_group=None):
    subgroup = 'C1C(C(OC1N2C=NC3=C(N=CN=C32)N)CO)O'
    substrate_mol = Chem.MolFromSmiles(substrate)
    subgroup_mol = Chem.MolFromSmiles(subgroup)
    replacement_mol = Chem.MolFromSmiles(replacement_group) if replacement_group else None

    matches = substrate_mol.GetSubstructMatches(subgroup_mol, uniquify=True)
    data = []

    for match in matches:
        editable_mol = Chem.RWMol(Chem.Mol(substrate_mol))
        atoms_to_remove = []

        for atom_idx in match:
            atom = editable_mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in match:
                    break
            else:
                atoms_to_remove.append(atom_idx)

        for atom_idx in sorted(atoms_to_remove, reverse=True):
            editable_mol.RemoveAtom(atom_idx)

        modified_molecule = editable_mol.GetMol()
        modified_smiles = Chem.MolToSmiles(modified_molecule, isomericSmiles=True)
        data.append({
            'substrate_smiles': substrate,
            'smiles': modified_smiles,
            'reaction': 'loss of deoxyadenosine',
            'reaction_rules': subgroup,
            'reaction_id': 'R59'
        })

    if not data:
        data.append({
            'substrate_smiles': substrate,
            'smiles': substrate,
            'reaction': 'loss of deoxyadenosine',
            'reaction_rules': subgroup,
            'reaction_id': 'R59'
        })

    return pd.DataFrame(data)



def expand_reaction_rules(df):
    # Using str.split() to split the 'ReactionRule' entries into lists of rules
    df['ReactionRule'] = df['ReactionRule'].str.split(', ')
    # Exploding the lists into separate rows
    df_expanded = df.explode('ReactionRule')
    df_expanded['mass_difference'] = df_expanded['Reaction'].astype(str).str[1]
    df_expanded = df_expanded[df_expanded['ReactionRule'].str.len() >= 1]

    return df_expanded


def summarize_reaction_data(substrate, rxn_rules_df):
    """
    Generates a summary table of reactions based on given substrate and a DataFrame of reaction rules.

    :param substrate: The substrate molecule for the reactions.
    :param rxn_rules_df: A DataFrame containing the reaction rules and associated metadata.
    :return: A list of dictionaries, each representing the reaction data for reactions that produce valid products.
    """
    summary_table = []

    # Iterate over each row in the DataFrame
    for index, rule in rxn_rules_df.iterrows():
        smirks = rule['ReactionRule']
        mass_difference = rule['mass_difference']
        description = rule['Description']
        reaction_id = rule['ID']

        # Predict products using the SMIRKS pattern
        all_products = predict_products(substrate, smirks)

        # Filter products by the specified mass difference
        products = filter_by_mass_difference(all_products, substrate, mass_difference)

        if products:
            for product in products:
                summary_table.append({
                    'substrate_smiles': substrate,
                    'smiles': product,
                    'reaction': description,
                    'reaction_rules': smirks,
                    'reaction_id': reaction_id
                })

    result_df = pd.DataFrame(summary_table)
    return result_df



def pipeline_function(substrate, rxn_rules_df):
    # Step 1: Execute each function
    try:
        df1 = summarize_reaction_data(substrate, rxn_rules_df)
        df2 = predict_products_group(substrate)
        df3 = drop_group(substrate)
    except Exception as e:
        # Handle potential errors in function execution
        print(f"Error during processing: {e}")
        return None

    # Step 2: Concatenate all DataFrames
    result_df = pd.concat([df1, df2, df3], axis=0, ignore_index=True)

    return result_df


def process_substrate_df(substrate_df, rxn_rules_df):
    # Container for the resulting DataFrames
    results = []

    # Loop through each row in the substrate DataFrame with a progress bar
    for index, row in tqdm(substrate_df.iterrows(), total=substrate_df.shape[0], desc="Processing substrates"):
        substrate = row['smiles']
        substrate_name = row['name']

        # Call the pipeline function
        result_df = pipeline_function(substrate, rxn_rules_df)

        if result_df is not None:
            # Add the substrate_name column
            result_df['substrate_name'] = substrate_name
            results.append(result_df)

    # Concatenate all results into one DataFrame if there are any results
    if results:
        final_df = pd.concat(results, ignore_index=True)
        return final_df
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no results

"""
 ____                      _            _ _
|  _ \ _   _ _ __    _ __ (_)_ __   ___| (_)_ __   ___
| |_) | | | | '_ \  | '_ \| | '_ \ / _ \ | | '_ \ / _ \
|  _ <| |_| | | | | | |_) | | |_) |  __/ | | | | |  __/
|_| \_\\__,_|_| |_| | .__/|_| .__/ \___|_|_|_| |_|\___|
                    |_|     |_|
"""
rxn_rules = expand_reaction_rules(rxn_rules)
result_df = process_substrate_df(initial_db, rxn_rules)


"""
 ____
|  _ \ ___ _ __ ___   _____   _____
| |_) / _ \ '_ ` _ \ / _ \ \ / / _ \
|  _ <  __/ | | | | | (_) \ V /  __/
|_| \_\___|_| |_| |_|\___/ \_/ \___|
"""
def process_smiles(input_df):
    """
    Process a DataFrame to remove invalid SMILES entries and compute chemical properties.

    Parameters:
        input_df (pd.DataFrame): Input DataFrame with a column 'smiles'.

    Returns:
        pd.DataFrame: Enhanced DataFrame with new chemical properties.
        int: Count of unique inchikey_first_block.
    """
    # Ensure SMILES data is treated as a string and handle missing values
    input_df['smiles'] = input_df['smiles'].astype(str).replace('nan', None)

    # Function to check if a SMILES string can be converted to a valid molecule
    def is_valid_smiles(smi):
        if smi is None:
            return False
        return Chem.MolFromSmiles(smi) is not None

    # Remove rows with invalid or missing SMILES strings
    n_before = len(input_df)
    valid_smiles_df = input_df[input_df['smiles'].apply(is_valid_smiles)].copy()
    n_invalid = n_before - len(valid_smiles_df)
    if n_invalid > 0:
        print(f"  Removed {n_invalid:,} rows with invalid product SMILES")

    # Function to compute various chemical properties from SMILES string
    def get_chem_properties(smiles):
        try:
            molecule = Chem.MolFromSmiles(smiles)
            inchi = Chem.MolToInchi(molecule)
            inchi_key = Chem.MolToInchiKey(molecule)
            formula = rdMolDescriptors.CalcMolFormula(molecule)
            exact_mass = Descriptors.ExactMolWt(molecule)
            inchi_key_first_block = inchi_key.split('-')[0]
            return pd.Series([inchi, inchi_key, inchi_key_first_block, formula, exact_mass])
        except Exception as e:
            return pd.Series([None, None, None, None, None])

    # Apply the function and add new columns to the DataFrame with progress bar
    results = []
    for smi in tqdm(valid_smiles_df['smiles'], desc="Generating product properties"):
        results.append(get_chem_properties(smi))
    valid_smiles_df[['inchi', 'inchikey', 'inchikey_first_block', 'formula', 'exact_mass']] = pd.DataFrame(results, index=valid_smiles_df.index)

    # Count unique 'inchikey_first_block' values
    unique_inchi_key_count = valid_smiles_df['inchikey_first_block'].nunique()

    return valid_smiles_df, unique_inchi_key_count


"""
 _____ _                                _
|_   _| |__   ___ _ __ _ __ ___   ___  | |
  | | | '_ \ / _ \ '__| '_ ` _ \ / _ \ | |
  | | | | | |  __/ |  | | | | | | (_) || |
  |_| |_| |_|\___|_|  |_| |_| |_|\___(_)_|
"""
def get_molecular_features(smiles):
    """Extract molecular features relevant for thermodynamic estimation."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        atom_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

        features = {
            'mw': Descriptors.MolWt(mol),
            'n_atoms': mol.GetNumAtoms(),
            'n_heavy': mol.GetNumHeavyAtoms(),
            'n_C': atom_counts.get('C', 0),
            'n_H': atom_counts.get('H', 0),
            'n_O': atom_counts.get('O', 0),
            'n_N': atom_counts.get('N', 0),
            'n_P': atom_counts.get('P', 0),
            'n_S': atom_counts.get('S', 0),
            'n_hdonors': Descriptors.NumHDonors(mol),
            'n_hacceptors': Descriptors.NumHAcceptors(mol),
            'n_rotatable': Descriptors.NumRotatableBonds(mol),
            'n_rings': rdMolDescriptors.CalcNumRings(mol),
            'n_aromatic_rings': rdMolDescriptors.CalcNumAromaticRings(mol),
            'tpsa': Descriptors.TPSA(mol),
            'logp': Descriptors.MolLogP(mol),
            'formal_charge': Chem.GetFormalCharge(mol),
        }
        return features
    except Exception:
        return None


def estimate_reaction_dg(substrate_smiles, product_smiles,
                         favorable_threshold=-30.0, unfavorable_threshold=30.0):
    """
    Estimate dG for a reaction based on structural changes using a group contribution approach.

    Returns: (dG_estimate, uncertainty, feasibility_category, structural_changes)
    """
    sub_feat = get_molecular_features(substrate_smiles)
    prod_feat = get_molecular_features(product_smiles)

    if sub_feat is None:
        return None, None, "Error", "Invalid substrate SMILES"
    if prod_feat is None:
        return None, None, "Error", "Invalid product SMILES"

    delta = {k: prod_feat[k] - sub_feat[k] for k in sub_feat}

    dG = 0.0
    notes = []

    if delta['n_H'] < -1:
        dG += abs(delta['n_H']) * 30
        notes.append(f"Dehydrogenation (-{abs(delta['n_H'])}H)")
    elif delta['n_H'] > 1:
        dG -= delta['n_H'] * 25
        notes.append(f"Hydrogenation (+{delta['n_H']}H)")

    if delta['n_O'] > 0:
        dG -= delta['n_O'] * 50
        notes.append(f"Oxidation (+{delta['n_O']}O)")
    elif delta['n_O'] < 0:
        dG += abs(delta['n_O']) * 40
        notes.append(f"Reduction (-{abs(delta['n_O'])}O)")

    if delta['n_P'] < 0:
        dG -= abs(delta['n_P']) * 30
        notes.append(f"Dephosphorylation (-{abs(delta['n_P'])}P)")
    elif delta['n_P'] > 0:
        dG += delta['n_P'] * 30
        notes.append(f"Phosphorylation (+{delta['n_P']}P)")

    if delta['n_C'] < 0:
        dG -= abs(delta['n_C']) * 20
        notes.append(f"Carbon loss (-{abs(delta['n_C'])}C)")
    elif delta['n_C'] > 0:
        dG += delta['n_C'] * 25
        notes.append(f"Carbon addition (+{delta['n_C']}C)")

    if delta['n_N'] != 0:
        dG += delta['n_N'] * 15
        notes.append(f"Nitrogen change ({delta['n_N']:+d}N)")

    if delta['formal_charge'] != 0:
        dG += abs(delta['formal_charge']) * 10
        notes.append(f"Charge change ({delta['formal_charge']:+d})")

    if delta['n_rings'] > 0:
        dG -= delta['n_rings'] * 15
        notes.append(f"Ring formation (+{delta['n_rings']} rings)")
    elif delta['n_rings'] < 0:
        dG += abs(delta['n_rings']) * 20
        notes.append(f"Ring breaking (-{abs(delta['n_rings'])} rings)")

    uncertainty = max(abs(dG) * 0.3, 20.0)

    if dG < favorable_threshold:
        feasibility = "Favorable"
    elif dG < 0:
        feasibility = "Slightly favorable"
    elif dG < unfavorable_threshold:
        feasibility = "Near equilibrium"
    elif dG < unfavorable_threshold * 2:
        feasibility = "Slightly unfavorable"
    else:
        feasibility = "Unfavorable"

    return dG, uncertainty, feasibility, "; ".join(notes) if notes else "Minimal structural change"


def add_thermodynamics(df):
    """Apply thermodynamic estimation to each row of the DataFrame."""
    results = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Estimating thermodynamics"):
        dg, unc, feas, changes = estimate_reaction_dg(row['substrate_smiles'], row['smiles'])
        results.append({
            'dg_estimate_kj_mol': dg,
            'dg_uncertainty_kj_mol': unc,
            'thermodynamic_feasibility': feas,
            'structural_changes': changes
        })
    thermo_df = pd.DataFrame(results)
    return pd.concat([df.reset_index(drop=True), thermo_df], axis=1)


"""
 ____            _                                        _
|  _ \ ___  ___| |_      _ __  _ __ ___   ___ ___  ___ ___(_)_ __   __ _
| |_) / _ \/ __| __|____| '_ \| '__/ _ \ / __/ _ \/ __/ __| | '_ \ / _` |
|  __/ (_) \__ \ ||_____| |_) | | | (_) | (_|  __/\__ \__ \ | | | | (_| |
|_|   \___/|___/\__|    | .__/|_|  \___/ \___\___||___/___/_|_| |_|\__, |
                         |_|                                        |___/
"""
# Post-process pipeline results to filter invalid entries and compute chemical properties
n_raw = len(result_df)
processed_data, unique_count = process_smiles(result_df)
n_valid = len(processed_data)

# Remove products already in substrate database (by inchikey_main_block)
substrate_inchikeys = set(initial_db['inchikey_main_block'].dropna().tolist())
before_substrate_filter = len(processed_data)
filtered_data = processed_data[~processed_data['inchikey_first_block'].isin(substrate_inchikeys)]
n_removed_existing = before_substrate_filter - len(filtered_data)
n_after_substrate_filter = len(filtered_data)
print(f"  Removed {n_removed_existing:,} products already in substrate database")

# Add thermodynamic estimation columns
filtered_data = add_thermodynamics(filtered_data)

# Remove rows where thermodynamic feasibility is Unfavorable
n_unfavorable = len(filtered_data[filtered_data['thermodynamic_feasibility'] == 'Unfavorable'])
filtered_data = filtered_data[filtered_data['thermodynamic_feasibility'] != 'Unfavorable']
n_final = len(filtered_data)
print(f"  Removed {n_unfavorable:,} thermodynamically unfavorable products")

# Save all unique products to parquet
filtered_data.to_parquet(outputfile, index=False)


"""
 ____  _        _   _     _   _
/ ___|| |_ __ _| |_(_)___| |_(_) ___ ___
\___ \| __/ _` | __| / __| __| |/ __/ __|
 ___) | || (_| | |_| \__ \ |_| | (__\__ \\
|____/ \__\__,_|\__|_|___/\__|_|\___|___/
"""
print("\n" + "=" * 70)
print("PIPELINE STATISTICS")
print("=" * 70)

# Basic counts
n_substrates = result_df['substrate_smiles'].nunique() if not result_df.empty else 0
n_reactions = result_df['reaction'].nunique() if not result_df.empty else 0

print(f"\nInput substrates:              {n_substrates:,}")
print(f"Distinct reaction types:       {n_reactions:,}")
print(f"\nTotal reactions generated:     {n_raw:,}")
print(f"After removing invalid SMILES: {n_valid:,}")
print(f"After removing existing in DB: {n_after_substrate_filter:,}")
print(f"After removing unfavorable:    {n_final:,}")

# Reaction type breakdown
if not result_df.empty and 'reaction' in result_df.columns:
    print(f"\nReaction type breakdown (before filtering):")
    rxn_counts = result_df['reaction'].value_counts()
    for rxn_name, count in rxn_counts.items():
        print(f"  {rxn_name}: {count:,}")

# Thermodynamic feasibility distribution (of final output)
if 'thermodynamic_feasibility' in filtered_data.columns and len(filtered_data) > 0:
    print(f"\nThermodynamic feasibility distribution (final output):")
    feas_counts = filtered_data['thermodynamic_feasibility'].value_counts()
    order = ['Favorable', 'Slightly favorable', 'Near equilibrium',
             'Slightly unfavorable', 'Error']
    for cat in order:
        if cat in feas_counts.index:
            count = feas_counts[cat]
            pct = count / len(filtered_data) * 100
            print(f"  {cat}: {count:,} ({pct:.1f}%)")

    # dG summary stats
    valid_dg = filtered_data['dg_estimate_kj_mol'].dropna()
    if len(valid_dg) > 0:
        print(f"\ndG distribution (kJ/mol):")
        print(f"  Mean:   {valid_dg.mean():.1f}")
        print(f"  Median: {valid_dg.median():.1f}")
        print(f"  Std:    {valid_dg.std():.1f}")
        print(f"  Min:    {valid_dg.min():.1f}")
        print(f"  Max:    {valid_dg.max():.1f}")

print(f"\nOutput saved to: {outputfile}")
print("=" * 70)
