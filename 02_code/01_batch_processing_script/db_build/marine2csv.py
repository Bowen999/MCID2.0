import pandas as pd
import requests

# Load the dataset
file_path = '/Users/bowen/Desktop/MCID2.0/01_source_data/marine/CMNPD_1.0_calc_prop.tsv'
data = pd.read_csv(file_path, delimiter='\t')

# Adding 'InChIKey first block' column
data['InChIKey first block'] = data['STANDARD_INCHI_KEY'].apply(lambda x: x.split('-')[0])

# Function to fetch compound name and synonyms from PubChem with error handling
def fetch_pubchem_data_safe(inchikey):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    properties_url = f"{base_url}/compound/inchikey/{inchikey}/property/IUPACName/JSON"
    synonyms_url = f"{base_url}/compound/inchikey/{inchikey}/synonyms/JSON"
    name, synonyms = None, []

    try:
        # Fetching the name
        name_response = requests.get(properties_url)
        if name_response.status_code == 200:
            name_data = name_response.json()
            name = name_data['PropertyTable']['Properties'][0]['IUPACName']
    except Exception as e:
        print(f"Error fetching name for {inchikey}: {e}")

    try:
        # Fetching the synonyms
        synonyms_response = requests.get(synonyms_url)
        if synonyms_response.status_code == 200:
            synonyms_data = synonyms_response.json()
            synonyms = synonyms_data['InformationList']['Information'][0]['Synonym']
    except Exception as e:
        print(f"Error fetching synonyms for {inchikey}: {e}")   

    return name, synonyms

# Adding data to dataframe
data['Name'] = None
data['Synonyms'] = None

# Uncomment the following block to fetch data from PubChem
for index, row in data.iterrows():
    name, synonyms = fetch_pubchem_data_safe(row['STANDARD_INCHI_KEY'])
    data.at[index, 'Name'] = name
    data.at[index, 'Synonyms'] = synonyms

# Deduplicate based on the 'InChIKey first block'
deduplicated_data = data.drop_duplicates(subset=['InChIKey first block'], keep='first')

# Save or process the deduplicated data
deduplicated_data.to_csv('/Users/bowen/Desktop/MCID2.0/03_initial_result/marine/cmnpd_0_rxn.csv', index=False)

# Optional: display a few rows to check the final output
print(deduplicated_data.head())
