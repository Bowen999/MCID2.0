import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
from tqdm import tqdm


hmdb_class_file = '/Users/bowen/Desktop/MCID2.0/01_source_data/lipid/HMDB_36_classyfire_21_annotations.csv'
class_id = 'CHEMONTID:0000001'
result_path = '/Users/bowen/Desktop/inorganic.csv'



# Read HMDB Class, and only keep the row that in the given class
human_class = pd.read_csv(hmdb_class_file)
hmdb_inorganic = human_class[human_class['Smiles']==class_id]['CompoundID'].to_list()
print(hmdb_inorganic)


# Get the information (SMIELS, InChIkey ....)
def fetch_hmdb_data(hmdb_ids, retries=3, delay=5, timeout=20):
    base_url = "https://hmdb.ca/metabolites/"
    data = []
    failed_ids = []

    for hmdb_id in tqdm(hmdb_ids, desc="Fetching HMDB Data"):
        attempt = 0
        success = False
        while attempt < retries and not success:
            try:
                url = f"{base_url}{hmdb_id}"
                response = requests.get(url, timeout=timeout)
                if response.status_code == 200:
                    soup = BeautifulSoup(response.content, 'html.parser')
                    common_name = soup.find('th', text='Common Name').find_next_sibling('td').text.strip()
                    smiles = soup.find('th', text='SMILES').find_next_sibling('td').text.strip()
                    inchi = soup.find('th', text='InChI Identifier').find_next_sibling('td').text.strip()
                    inchi_key = soup.find('th', text='InChI Key').find_next_sibling('td').text.strip()
                    
                    data.append({
                        'HMDB ID': hmdb_id,
                        'URL': url,
                        'Common Name': common_name,
                        'SMILES': smiles,
                        'InChI Identifier': inchi,
                        'InChI Key': inchi_key
                    })
                    success = True
                else:
                    print(f"Attempt {attempt + 1}: Failed to retrieve data for {hmdb_id}. Status code: {response.status_code}")
                    attempt += 1
                    time.sleep(delay)
            except requests.exceptions.Timeout:
                print(f"Attempt {attempt + 1}: Timeout while fetching data for {hmdb_id}. Retrying...")
                attempt += 1
                time.sleep(delay)

        if not success:
            print(f"Failed to retrieve data after {retries} attempts for {hmdb_id}.")
            failed_ids.append(hmdb_id)

    return pd.DataFrame(data), failed_ids



# Result
df, failed_ids = fetch_hmdb_data(hmdb_inorganic)
df.to_csv(result_path, index=False)
print(df)
print("Failed IDs:", failed_ids)


