from flask import Flask, render_template, send_from_directory, request, jsonify, abort
import os
import pandas as pd
import numpy as np
from flask_cors import CORS
import mycompoundid as mcid
from flask import Flask, send_from_directory


# *************** Part1: Router 路由 **********************#
# Initialize Flask app with specific folders for static files and templates
# 初始化 Flask 应用，指定静态文件和模板的文件夹
app = Flask(__name__, static_folder='frontend/build/static', template_folder="frontend/build")
CORS(app)

@app.route('/test') #Test page
def hello():
    return "test"


# Route for serving the main index page
# 用于提供主索引页面的路由
@app.route('/')
def index():
    return render_template('index.html')



# Catch-all route for serving files from the frontend build directory
# 从前端构建目录提供文件的全捕获路由
@app.route('/<path:path>')  # Catch-all route for any other paths
def serve(path):
    # First check if the path corresponds to a file in the frontend/build directory
    if path != "" and os.path.exists(os.path.join('frontend/build', path)):
        return send_from_directory('frontend/build', path)
    else:
        # If the file is not found in either location, return the index.html or a 404 page
        return render_template('index.html')  # or return "Page not found", 404
  
  
    
# Route for serving custom HTML pages for compounds information
# 用于提供化合物信息自定义 HTML 页面的路由
@app.route('/mcid<page>')
def custom_html(page):
    filepath = os.path.join('compounds_html', f'mcid{page}.html')
    if os.path.exists(filepath):
        return send_from_directory('compounds_html', f'mcid{page}.html')
    else:
        abort(404)  # Page not found




# Function for handling search queries
# 用于处理搜索查询的函数和路由
def search(input, tol, tol_unit, num_rxn, adduct):
    print('input:', input)
    print('tol:', tol)
    print('tol_unit:', tol_unit)
    print('num_rxn:', num_rxn)
    print('adduct:', adduct)
    return pd.DataFrame({'input': input, 'tol': tol, 'tol_unit': tol_unit,'num_rxn': num_rxn, 'adduct': adduct}, index=[0])

# @app.route('/api/search', methods=['POST'])
# def search_api():
#     data = request.get_json()

#     # Access keys with default values
#     input_val = data.get('input', '')  # default to an empty string
#     tol = data.get('tol', '1')  # default to '1'
#     tol = float(tol)
#     tol_unit = data.get('tol_unit', 'Da')  # default to 'Da'
#     num_rxn = data.get('num_rxn', 'no_reaction')  # default to 'no_reaction'
#     adduct = data.get('adduct', 'neutral')  # default to 'neutral'

#     print('input:', input_val, type(input_val))
#     print('tol:', tol, type(tol))
#     print('tol_unit:', tol_unit, type(tol_unit))
#     print('num_rxn:', num_rxn, type(num_rxn))
#     print('adduct:', adduct, type(adduct))
    
#     # Calling the mcid.compound_id function with the gathered data
#     result = mcid.compound_id(input_val, adduct, tol, tol_unit, num_rxn)
#     result['Hits'] = result['Hits'].str.replace('\n', '<br></br>')
#     table_html = result.to_html(classes='data', index=False, escape=False)
    
#     # print(table_html)
#     # table_html = result.to_html(classes='data', index=False)
#     return jsonify({'tableHTML': table_html})
@app.route('/api/search', methods=['POST'])
def search_api():
    try:
        data = request.get_json()

        # Access keys with default values
        input_val = data.get('input', '')  # default to an empty string
        tol = data.get('tol', '1')  # default to '1'
        tol = float(tol)
        tol_unit = data.get('tol_unit', 'Da')  # default to 'Da'
        num_rxn = data.get('num_rxn', 'no_reaction')  # default to 'no_reaction'
        adduct = data.get('adduct', 'neutral')  # default to 'neutral'

        # Calling the mcid.compound_id function with the gathered data
        result = mcid.compound_id(input_val, adduct, tol, tol_unit, num_rxn)
        result['Hits'] = result['Hits'].str.replace('\n', '<br></br>')

        if result.empty:
            raise ValueError("No data found for the provided inputs.")

        table_html = result.to_html(classes='data', index=False, escape=False)

        return jsonify({'tableHTML': table_html})
    except Exception as e:
        # Log the exception if needed
        print(f"An error occurred: {e}")
        return jsonify({'error': str(e)}), 400


if __name__ == '__main__':
    app.run(host="0.0.0.0", port=5000, debug=True)