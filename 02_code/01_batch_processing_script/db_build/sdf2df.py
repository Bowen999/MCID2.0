# import pandas as pd

# # Path to the data file
file_path = '/Users/bowen/Desktop/MCID2.0/01_source_data/human/structures-2.sdf'
result_path = '/Users/bowen/Desktop/MCID2.0/03_initial_result/human/human2.csv'


import pandas as pd

def parse_and_process_sdf(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    
    compounds = content.split('M  END')
    
    data = []
    
    for compound in compounds:
        if compound.strip() == '':
            continue
        compound_data = {}
        lines = compound.split('\n')
        
        for i, line in enumerate(lines):
            if line.startswith('> <DATABASE_ID>'):
                compound_data['Entry'] = lines[i + 1].strip()
            elif line.startswith('> <GENERIC_NAME>'):
                compound_data['Name'] = lines[i + 1].strip()
            elif line.startswith('> <JCHEM_TRADITIONAL_IUPAC>'):
                synonym = lines[i + 1].strip()
                # Check if there is another synonym on the next tag
                if '> <JCHEM_IUPAC>' in lines[i + 2]:
                    synonym += ', ' + lines[i + 3].strip()
                compound_data['Synonyms'] = synonym if synonym else ''
            elif line.startswith('> <EXACT_MASS>'):
                try:
                    compound_data['Exact mass'] = float(lines[i + 1].strip())
                except ValueError:
                    compound_data['Exact mass'] = 0
                    # continue
            elif line.startswith('> <FORMULA>'):
                compound_data['Formula'] = lines[i + 1].strip()
            elif line.startswith('> <SMILES>'):
                compound_data['SMILES'] = lines[i + 1].strip()
            elif line.startswith('> <INCHI_IDENTIFIER>'):
                compound_data['InChI'] = lines[i + 1].strip()
            elif line.startswith('> <INCHI_KEY>'):
                inchi_key = lines[i + 1].strip()
                compound_data['InChIKey'] = inchi_key
                compound_data['InChIKey first block'] = inchi_key.split('-')[0]
        
        data.append(compound_data)
    
    df = pd.DataFrame(data)
    df = df.dropna(subset=['Exact mass'])
    print(df.shape[0])

    # Deduplicate based on 'InChIKey first block' by merging 'Name' into 'Synonyms' of the first occurrence
    def aggregate_names(group):
        group['Synonyms'] = group['Synonyms'].fillna('')  # Replace NaN with empty strings
        if len(group) > 1:
            first_row = group.iloc[0]
            subsequent_names = group['Name'].iloc[1:].tolist()
            first_row['Synonyms'] += ', ' + ', '.join(subsequent_names)
            return first_row
        return group.iloc[0]

    df = df.groupby('InChIKey first block').apply(aggregate_names).reset_index(drop=True)

    # Count unique 'InChIKey first block'
    unique_inchi_first_block = df['InChIKey first block'].nunique()
    
    return df, unique_inchi_first_block



result, unique_inchi_first_block = parse_and_process_sdf(file_path)
print("Number of unique 'InChIKey first block':", unique_inchi_first_block)
print()
print(result.head())
result.to_csv(result_path, index=False)
