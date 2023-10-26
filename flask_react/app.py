from flask import Flask, render_template, send_from_directory, request, jsonify
import os
import pandas as pd
import numpy as np
from flask_cors import CORS
import mycompoundid as mcid

app = Flask(__name__, static_folder='frontend/build/static', template_folder="frontend/build")
CORS(app)


@app.route('/')
def index():
    return render_template('index.html')



# Search function
def search(input, tol, tol_unit, num_rxn, adduct):
    print('input:', input)
    print('tol:', tol)
    print('tol_unit:', tol_unit)
    print('num_rxn:', num_rxn)
    print('adduct:', adduct)
    return pd.DataFrame({'input': input, 'tol': tol, 'tol_unit': tol_unit,'num_rxn': num_rxn, 'adduct': adduct}, index=[0])



@app.route('/api/search', methods=['POST'])
def search_api():
    data = request.get_json()

    # Access keys with default values
    input_val = data.get('input', '')  # default to an empty string
    tol = data.get('tol', '1')  # default to '1'
    tol = float(tol)
    tol_unit = data.get('tol_unit', 'Da')  # default to 'Da'
    num_rxn = data.get('num_rxn', 'no_reaction')  # default to 'no_reaction'
    adduct = data.get('adduct', 'neutral')  # default to 'neutral'

    print('input:', input_val, type(input_val))
    print('tol:', tol, type(tol))
    print('tol_unit:', tol_unit, type(tol_unit))
    print('num_rxn:', num_rxn, type(num_rxn))
    print('adduct:', adduct, type(adduct))
    
    # Calling the mcid.compound_id function with the gathered data
    result = mcid.compound_id(input_val, adduct, tol, tol_unit, num_rxn)
    result['Hits'] = result['Hits'].str.replace('\n', '<br></br>')
    table_html = result.to_html(classes='data', index=False, escape=False)

    # table_html = result.to_html(classes='data', index=False)
    return jsonify({'tableHTML': table_html})





# Test function
def test(a, b):
    data = np.full((a, b), 'X')  # Create a matrix with `a` rows and `b` columns filled with 'X'
    return pd.DataFrame(data)

@app.route('/api/test', methods=['POST'])
def compute_test():
    data = request.get_json()
    a = int(data['firstValue'])
    b = int(data['secondValue'])
    df = test(a, b)
    # Convert dataframe to HTML
    return df.to_html(classes="table table-bordered")




@app.route('/<path:path>') # For any other routes, serve the index.html as a catch-all
def serve_file(path):
    return send_from_directory('frontend/build', path)

@app.route('/api/sum', methods=['POST'])
def compute_sum():
    data = request.get_json()
    result = data['firstValue'] + data['secondValue']
    return jsonify({'result': result})

if __name__ == '__main__':
    app.run(debug=True)
