# How to Create a Flask based Web APP


## 1.Create a flask project

### Step 1: Install Flask

```
pip install Flask
```
### Step 2: Creating the App

**app.js:**

```
from flask import Flask

app = Flask(__name__)

@app.route('/')
def hello_world():
    return 'Hello, World!'

if __name__ == '__main__':
    app.run(debug=True)
```
This code creates a new Flask web server. When you navigate to the root URL (/), it will return "Hello, World!"

### Step 3: Running the App

```
python app.py
```


### Step 4: Add More route:

**app.py:**

```
from flask import Flask
from flask import Flask, render_template

app = Flask(__name__)

@app.route('/')
def hello_world():
    return 'Hello, World!'


@app.route('/about')
def about():
    return 'About Page'

@app.route('/html')
def contact():
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
```

Create a new folder named **templates** in your project directory.  
Inside the templates folder, create a new HTML file named **index.html:**

```
<!DOCTYPE html>
<html>
<head>
    <title>My Flask App</title>
</head>
<body>
    <h1>Welcome to My Flask App!</h1>
    <p>This is a simple Flask app.</p>
</body>
</html>
```
## 2.Connect Flask (backend) with React (frontend)
#### 目标：新建一个Contact页面（使用React构建）

2.1 新建React项目
```
npx create-react-app frontend
```

2.2 **Contact.js:**

```
// frontend/src/Contact.js

import React from 'react';

function Contact() {
    return (
        <div>
            <h2>Contact Page</h2>
            <p>This is the contact page.</p>
        </div>
    );
}

export default Contact;
```

2.3 integrate **Contact** component in **App.js**:

```
// frontend/src/App.js

import React from 'react';
import './App.css';
import Contact from './Contact';

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <Contact />
      </header>
    </div>
  );
}

export default App;
```

2.4 Build the React Frontend

```
npm run build
```

2.5 **app.py**:

```
from flask import Flask, render_template, send_from_directory
import os

app = Flask(__name__, static_folder='frontend/build/static', template_folder="frontend/build") #连接React

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/<path:path>') #Router连接React
def serve_file(path):
    return send_from_directory('frontend/build', path)

if __name__ == '__main__':
    app.run(debug=True)
```

## 3.Add one more page in React
#### 目标：新加一个页面Input，可以输入两个值（A，B）到表单，并且提交
3.1 Setting up React Router  
如果需要在不同的页面之间切换，需要一个React Route

```
npm install react-router-dom
```

3.2 新建一个**Input.js**

```
// frontend/src/Input.js

import React, { useState } from 'react';

function Input() {
    const [firstValue, setFirstValue] = useState('');
    const [secondValue, setSecondValue] = useState('');

    const handleSubmit = (e) => {
        e.preventDefault(); // Prevent the default form submission behavior（这行代码的目的是阻止表单的默认提交行为，这样页面就不会重新加载，也不会发送任何数据到服务器。）
        alert(`First value: ${firstValue}, Second value: ${secondValue}`);
    };

    return (
        <div>
            <h2>Input Page</h2>
            <form onSubmit={handleSubmit}>
                <input 
                    type="number" 
                    placeholder="Enter first value" 
                    value={firstValue} 
                    onChange={e => setFirstValue(e.target.value)}
                />
                <input 
                    type="number" 
                    placeholder="Enter second value" 
                    value={secondValue} 
                    onChange={e => setSecondValue(e.target.value)}
                />
                <button type="submit">Submit</button>
            </form>
        </div>
    );
}

export default Input;
```
*useState是React的一个Hook写法，const [count, setCount] = useState(0); 来设置变量*

form的写法

```
<h2>HTML Forms</h2>

<form action="/action_page.php">
  <label for="fname">First name:</label><br>
  <input type="text" id="fname" name="fname" value=" "><br><br>
  
  <label for="lname">Last name:</label><br>
  <input type="text" id="lname" name="lname" value="Doe"><br><br>

  <label for="age">Age:</label><br>
  <input type="number" id="age" name="age" value="18"><br><br>
  
  <label for="gender">Genger:</label><br>
  <input type="radio" id="male" name="gender" value="male">
  <label for="male">Male</label><br>
  <input type="radio" id="female" name="gender" value="female">
  <label for="female">Female</label><br><br>
  
  
  
  <input type="submit" value="Submit">
</form> 
</html>
```
<h2>HTML Forms</h2>

<form action="/action_page.php">
  <label for="fname">First name:</label><br>
  <input type="text" id="fname" name="fname" value=" "><br><br>
  
  <label for="lname">Last name:</label><br>
  <input type="text" id="lname" name="lname" value="Doe"><br><br>

  <label for="age">Age:</label><br>
  <input type="number" id="age" name="age" value="18"><br><br>
  
  <label for="gender">Genger:</label><br>
  <input type="radio" id="male" name="gender" value="male">
  <label for="male">Male</label><br>
  <input type="radio" id="female" name="gender" value="female">
  <label for="female">Female</label><br><br>
  
  
  
  <input type="submit" value="Submit">
</form> 
</html>


## 4.Connnect front end with back end
#### 目标：在input界面输入a，b。点击submit。跳转到结果页。结果页是一个a*b的dataframe。dataframe的生成在后端（Python）中完成

4.1 write the function to generate df (**app.py**):

```
from flask import Flask, render_template, send_from_directory, request, jsonify
import os
import pandas as pd
import numpy as np

app = Flask(__name__, static_folder='frontend/build/static', template_folder="frontend/build")

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/test', methods=['POST'])
def compute_test():
    data = request.get_json()
    a = int(data['firstValue'])
    b = int(data['secondValue'])
    df = test(a, b)
    # Convert dataframe to HTML
    return df.to_html(classes="table table-bordered")

def test(a, b):
    data = np.full((a, b), 'X')  # Create a matrix with `a` rows and `b` columns filled with 'X'
    return pd.DataFrame(data)

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
```
 

4.2 **Inpute.js**:

```
// frontend/src/Input.js

import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';

function Input() {
    const [firstValue, setFirstValue] = useState('');
    const [secondValue, setSecondValue] = useState('');
    const navigate = useNavigate();

    const handleSubmit = async (e) => {
        e.preventDefault();

        try {
            const response = await fetch('/api/test', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    firstValue: parseInt(firstValue),
                    secondValue: parseInt(secondValue)
                })
            });

            if (!response.ok) {
                throw new Error('Network response was not ok');
            }

            const data = await response.text();
            
            // Navigate to the result page and pass the data along
            navigate('/result', { state: { tableHTML: data } });

        } catch (error) {
            console.error("There was a problem with the fetch operation:", error.message);
        }
    };

    return (
        <div>
            <h2>Input Page</h2>
            <form onSubmit={handleSubmit}>
                <input 
                    type="number" 
                    placeholder="Enter first value" 
                    value={firstValue} 
                    onChange={e => setFirstValue(e.target.value)}
                />
                <input 
                    type="number" 
                    placeholder="Enter second value" 
                    value={secondValue} 
                    onChange={e => setSecondValue(e.target.value)}
                />
                <button type="submit">Submit</button>
            </form>
        </div>
    );
}

export default Input;
```

4.3 构建一个**Result.js**来展示结果：

```
import React from 'react';

function Result(props) {
    return (
        <div>
            <h2>Result Page</h2>
            <div dangerouslySetInnerHTML={{ __html: props.tableHTML }} />
        </div>
    );
}

export default Result;
```


4.4 更新**App.js**里的Router

```

import React from 'react';
import './App.css';
import Contact from './Contact';
import Input from './Input';
import Result from './Result';
import { BrowserRouter, Routes, Route, Link, Outlet, useLocation } from 'react-router-dom';

function App() {
  return (
    <BrowserRouter>
      <div className="App">
        <nav>
          <Link to="/">Home</Link>
          <Link to="/contact">Contact</Link>
          <Link to="/input">Input</Link>
        </nav>

        <Routes>
          <Route path="/" element={<Home />} />
          <Route path="/contact" element={<Contact />} />
          <Route path="/input" element={<Input />} />
          <Route path="/result" element={<ResultWrapper />} />
        </Routes>
      </div>
    </BrowserRouter>
  );
}

function Home() {
  return <div><h2>Welcome to the app!</h2></div>;
}

function ResultWrapper() {
  const location = useLocation();
  const tableHTML = location.state?.tableHTML;
  return <Result tableHTML={tableHTML} />;
}

export default App;
```

## 5. DO more
* WSGI server 

```
pip install gunicorn
gunicorn app:app
gunicorn app:app -b 0.0.0.0:8000
```

* CORS

```
from flask_cors import CORS

app = Flask(__name__, static_folder='frontend/build/static', template_folder="frontend/build")
CORS(app)

```

* Webpack

