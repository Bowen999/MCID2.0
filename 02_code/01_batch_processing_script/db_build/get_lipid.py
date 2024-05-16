import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
from tqdm import tqdm

hmdb_class_file = '/Users/bowen/Desktop/MCID2.0/01_source_data/lipid/HMDB_36_classyfire_21_annotations.csv'
class_id = 'CHEMONTID:0000012' # lipid
result_path = '/Users/bowen/Desktop/lipid.csv'



# Read HMDB Class, and only keep the row that in the given class
human_class = pd.read_csv(hmdb_class_file)
hmdb_lipid = human_class[human_class['Smiles']==class_id]['CompoundID'].to_list()
print(hmdb_lipid)


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
                    def find_value(tag, text):
                        element = soup.find(tag, text=text)
                        if element and element.find_next_sibling('td'):
                            return element.find_next_sibling('td').text.strip()
                        return None

                    common_name = find_value('th', 'Common Name')
                    smiles = find_value('th', 'SMILES')
                    inchi = find_value('th', 'InChI Identifier')
                    inchi_key = find_value('th', 'InChI Key')

                    if None not in [common_name, smiles, inchi, inchi_key]:
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
                        print(f"Data missing for {hmdb_id}, skipping...")
                        attempt += 1
                        time.sleep(delay)
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
df, failed_ids = fetch_hmdb_data(hmdb_lipid)
df.to_csc(result_path, index=False)
print(df)
print("Failed IDs:", failed_ids)
