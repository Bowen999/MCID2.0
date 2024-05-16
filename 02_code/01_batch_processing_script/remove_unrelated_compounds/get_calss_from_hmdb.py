import xml.etree.ElementTree as ET
import pandas as pd
from tqdm import tqdm

xml_file_path = '/Users/bowen/Documents/个人资料&申请/完结项目/LIME/LIME1/DB_Parsing/HMDB/hmdb_metabolites.xml'


# Parse the XML from the file
tree = ET.parse(xml_file_path)
root = tree.getroot()
print(tree)

# Create an empty list to hold each metabolite's data
data = []

# Get up to 100 metabolite elements
metabolites = root.findall('{http://www.hmdb.ca}metabolite')

# Iterate through each metabolite in the XML with a progress bar
for metabolite in tqdm(metabolites, desc="Processing Metabolites"):
    accession = metabolite.find('{http://www.hmdb.ca}accession').text
    name = metabolite.find('{http://www.hmdb.ca}name').text
    mono_mw = metabolite.find('{http://www.hmdb.ca}monisotopic_molecular_weight').text
    smiles = metabolite.find('{http://www.hmdb.ca}smiles').text
    inchikey = metabolite.find('{http://www.hmdb.ca}inchikey').text
    super_class = metabolite.find('{http://www.hmdb.ca}taxonomy/{http://www.hmdb.ca}super_class').text
    _class = metabolite.find('{http://www.hmdb.ca}taxonomy/{http://www.hmdb.ca}class').text
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
df.to_csv('lipid_in_hmdb.csv')