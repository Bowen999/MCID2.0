import sqlite3
import pandas as pd

def mz_list(input_string):
    # Split by newline and filter out empty strings
    str_values = filter(None, input_string.split('\n'))
    
    # Convert to integers or floats and return as a list
    numbers = []
    for value in str_values:
        try:
            # Try converting to float first
            num = float(value)
            # If it's an integer (i.e., has no fractional part), convert to int
            if num.is_integer():
                num = int(num)
            numbers.append(num)
        except ValueError:
            # If conversion fails, skip this value
            continue
            
    return numbers


def cal_adduct(mz_values, adduct_type):
    # Define mass adjustments for each adduct type
    adduct_masses = {
        "neutral": 0,
        "m-h": -1.007276,
        "m+h": 1.007276,
        "m+na": 22.989218,
        "m+k": 38.963158,
        "m+nh4": 18.033823
    }
    
    try:
        delta_mz = adduct_masses[adduct_type]
        exact_masses = [mz - delta_mz for mz in mz_values]
        return exact_masses
    except KeyError:
        return f"Error: {adduct_type} is not a valid adduct type."
    except Exception as e:
        return f"Error: {e}"



def search_in_db(exact_mass_list, tol, tol_unit, db):
    # Connect to the database
    conn = sqlite3.connect('mcid_1.db')
    cursor = conn.cursor()

    data = {"m/z": [], "Num. of Hits": [], "Hits": []}
    
    for exact_mass in exact_mass_list:
        # Calculate error based on tolerance unit
        if tol_unit == "Da":
            error = tol
        elif tol_unit == "ppm":
            error = (tol / 1e6) * exact_mass
        else:
            raise ValueError("Invalid tol_unit. It should be either 'Da' or 'ppm'.")
        
        # Set lower and upper bounds
        lower_bound = exact_mass - error
        upper_bound = exact_mass + error
        
        # Query the database, using double quotes to escape the table name
        query = f'SELECT * FROM "{db}" WHERE exact_mass BETWEEN ? AND ?'
        cursor.execute(query, (lower_bound, upper_bound))
        matches = cursor.fetchall()
        
        hits = []
        for match in matches:
            mcid = match[0]
            if db == "0_rxn":
                # Assuming the name is the second column in the match
                name = match[1]
                # Hyperlink only the mcid, not the name
                hit = f"<a href='https://www.kegg.jp/entry/{mcid}'>{mcid}</a> ({name})"
            else:
                # Hyperlink to a local resource using mcid
                hit = f"<a href='http://127.0.0.1:5000/compounds/{mcid}'>{mcid}</a>"
            hits.append(hit)
        
        # If there are no hits, append "No Match"
        if not hits:
            hits = ["No Match"]
        
        # Add to data
        data["m/z"].append(exact_mass)
        data["Num. of Hits"].append(len(matches))  # Use matches instead of hits for the count
        data["Hits"].append(", ".join(hits))
    
    # Close the database connection
    conn.close()    
    
    # Convert data to a DataFrame
    df = pd.DataFrame(data)
    
    return df



def compound_id(input, adduct_type, tol, tol_unit, db):
    mz_values = mz_list(input)
    exact_masses = cal_adduct(mz_values, adduct_type)
    
    result = search_in_db(exact_masses, tol, tol_unit, db)
    return result