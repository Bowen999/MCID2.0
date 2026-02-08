import xml.etree.ElementTree as ET
import pandas as pd
from tqdm import tqdm

# Define the path to your XML file
xml_file_path = '/Users/bowen/Documents/个人资料&申请/完结项目/LIME/LIME1/DB_Parsing/HMDB/hmdb_metabolites.xml'

# Parse the XML from the file
tree = ET.parse(xml_file_path)
root = tree.getroot()

# Create an empty list to hold each metabolite's data
data = []

# Get up to 100 metabolite elements
metabolites = root.findall('{http://www.hmdb.ca}metabolite')

# Iterate through each metabolite in the XML with a progress bar
for metabolite in tqdm(metabolites, desc="Processing Metabolites"):
    accession = metabolite.find('{http://www.hmdb.ca}accession').text if metabolite.find('{http://www.hmdb.ca}accession') is not None else ''
    name = metabolite.find('{http://www.hmdb.ca}name').text if metabolite.find('{http://www.hmdb.ca}name') is not None else ''
    mono_mw = metabolite.find('{http://www.hmdb.ca}monisotopic_molecular_weight').text if metabolite.find('{http://www.hmdb.ca}monisotopic_molecular_weight') is not None else ''
    smiles = metabolite.find('{http://www.hmdb.ca}smiles').text if metabolite.find('{http://www.hmdb.ca}smiles') is not None else ''
    inchikey = metabolite.find('{http://www.hmdb.ca}inchikey').text if metabolite.find('{http://www.hmdb.ca}inchikey') is not None else ''
    
    # Handle potentially missing taxonomy details with try-except or conditional check
    taxonomy_element = metabolite.find('{http://www.hmdb.ca}taxonomy')
    super_class = ''
    _class = ''
    if taxonomy_element is not None:
        super_class = taxonomy_element.find('{http://www.hmdb.ca}super_class').text if taxonomy_element.find('{http://www.hmdb.ca}super_class') is not None else ''
        _class = taxonomy_element.find('{http://www.hmdb.ca}class').text if taxonomy_element.find('{http://www.hmdb.ca}class') is not None else ''
    
    # Append this metabolite's data to the list
    data.append({
        'Accession': accession,
        'Name': name,
        'Monisotopic Molecular Weight': mono_mw,
        'SMILES': smiles,
        'InChIKey': inchikey,
        'Super Class': super_class,
        'Class': _class
    })

# Convert the list to a DataFrame
df = pd.DataFrame(data)

# Show the DataFrame
print(df)
df.to_csv('lipid_in_hmdb.csv')
