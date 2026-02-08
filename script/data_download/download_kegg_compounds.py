#!/usr/bin/env python3
"""
Download all compounds from KEGG database and save to CSV.

Output columns:
- compound_id: KEGG compound ID (e.g., C00001)
- name: Primary compound name
- name_all: All other names (semicolon-separated)
- exact_mass: Exact molecular mass
- smiles: SMILES string from MOL structure
- inchikey: InChIKey derived from SMILES
- inchikey_main_block: First block of InChIKey (connectivity layer)
- pathway: Associated pathways (format: "map00052: Galactose metabolism, ...")

Usage:
    python download_kegg_compounds.py

Notes:
    - Downloads ~19,571 compounds (as of 2024)
    - Takes approximately 7-8 hours due to API rate limits
    - Saves progress every 500 compounds to allow resumption
    - Requires: requests, pandas, rdkit, tqdm
"""

import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey
import time
import os
from tqdm import tqdm


def get_compound_list():
    """Fetch list of all compound IDs from KEGG."""
    print("Fetching compound list from KEGG...")
    response = requests.get("https://rest.kegg.jp/list/compound")
    response.raise_for_status()
    
    compound_ids = []
    for line in response.text.strip().split('\n'):
        compound_id = line.split('\t')[0].replace('cpd:', '')
        compound_ids.append(compound_id)
    
    print(f"Found {len(compound_ids)} compounds")
    return compound_ids


def get_compound_details(compound_id):
    """Fetch compound details from KEGG API."""
    result = {
        'compound_id': compound_id,
        'name': None,
        'name_all': None,
        'exact_mass': None,
        'smiles': None,
        'inchikey': None,
        'inchikey_main_block': None,
        'pathway': None
    }
    
    try:
        # Get compound info
        response = requests.get(f"https://rest.kegg.jp/get/{compound_id}")
        if response.status_code != 200:
            return result
        
        content = response.text
        lines = content.split('\n')
        
        # Parse the flat file format
        current_field = None
        names = []
        pathways = []
        
        for line in lines:
            if line.startswith('NAME'):
                current_field = 'NAME'
                name_part = line[12:].strip().rstrip(';')
                if name_part:
                    names.append(name_part)
            elif line.startswith('EXACT_MASS'):
                current_field = 'EXACT_MASS'
                result['exact_mass'] = line[12:].strip()
            elif line.startswith('PATHWAY'):
                current_field = 'PATHWAY'
                pathway_part = line[12:].strip()
                if pathway_part:
                    parts = pathway_part.split(maxsplit=1)
                    if len(parts) == 2:
                        pathways.append(f"{parts[0]}: {parts[1]}")
            elif line.startswith(' ') and current_field:
                # Continuation of previous field
                value = line.strip().rstrip(';')
                if current_field == 'NAME' and value:
                    names.append(value)
                elif current_field == 'PATHWAY' and value:
                    parts = value.split(maxsplit=1)
                    if len(parts) == 2:
                        pathways.append(f"{parts[0]}: {parts[1]}")
            elif line and not line.startswith(' '):
                current_field = None
        
        # Set names
        if names:
            result['name'] = names[0]
            if len(names) > 1:
                result['name_all'] = '; '.join(names[1:])
        
        # Set pathways
        if pathways:
            result['pathway'] = ', '.join(pathways)
        
        # Get MOL structure and convert to SMILES
        mol_response = requests.get(f"https://rest.kegg.jp/get/{compound_id}/mol")
        if mol_response.status_code == 200 and mol_response.text.strip():
            mol = Chem.MolFromMolBlock(mol_response.text)
            if mol:
                result['smiles'] = Chem.MolToSmiles(mol)
                inchikey = MolToInchiKey(mol)
                if inchikey:
                    result['inchikey'] = inchikey
                    result['inchikey_main_block'] = inchikey.split('-')[0]
        
    except Exception as e:
        print(f"Error processing {compound_id}: {e}")
    
    return result


def load_progress(checkpoint_file):
    """Load progress from checkpoint file if exists."""
    if os.path.exists(checkpoint_file):
        df = pd.read_csv(checkpoint_file)
        completed_ids = set(df['compound_id'].tolist())
        print(f"Loaded {len(completed_ids)} completed compounds from checkpoint")
        return df.to_dict('records'), completed_ids
    return [], set()


def save_progress(results, checkpoint_file):
    """Save current progress to checkpoint file."""
    df = pd.DataFrame(results)
    df.to_csv(checkpoint_file, index=False)


def main():
    # Configuration
    output_file = "/Users/bowen/Desktop/mcid/database/kegg_download.csv"
    checkpoint_file = "/Users/bowen/Desktop/mcid/database/kegg_compounds_checkpoint.csv"
    save_interval = 500  # Save progress every N compounds
    api_delay = 0.1  # Delay between API calls (seconds)
    
    # Get all compound IDs
    compound_ids = get_compound_list()
    
    # Load any existing progress
    results, completed_ids = load_progress(checkpoint_file)
    
    # Filter out already completed compounds
    remaining_ids = [cid for cid in compound_ids if cid not in completed_ids]
    print(f"Remaining compounds to download: {len(remaining_ids)}")
    
    # Download compounds
    for i, compound_id in enumerate(tqdm(remaining_ids, desc="Downloading compounds")):
        result = get_compound_details(compound_id)
        results.append(result)
        
        # Rate limiting
        time.sleep(api_delay)
        
        # Save checkpoint periodically
        if (i + 1) % save_interval == 0:
            save_progress(results, checkpoint_file)
            print(f"\nCheckpoint saved: {len(results)} compounds")
    
    # Save final results
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    print(f"\nDownload complete! Saved {len(df)} compounds to {output_file}")
    
    # Clean up checkpoint file
    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)
        print("Checkpoint file removed")
    
    # Print summary
    print("\n=== Summary ===")
    print(f"Total compounds: {len(df)}")
    print(f"With SMILES: {df['smiles'].notna().sum()}")
    print(f"With pathways: {df['pathway'].notna().sum()}")


if __name__ == "__main__":
    main()
