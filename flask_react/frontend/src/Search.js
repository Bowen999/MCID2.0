import React, { useState } from 'react';
import axios from 'axios';
import { useNavigate } from 'react-router-dom';
import { Input, Radio, Button, Modal } from 'antd';
import { Tooltip } from 'antd';
import { QuestionCircleOutlined } from '@ant-design/icons';
const { TextArea } = Input;
// import { Button } from 'antd';


// ***** 1. FUNCTION 功能 ********
// 1.1 INPUT 获取变量
function Search() {
  const [inputValue, setInputValue] = useState("28.0\n181.07389\n114.0667");
  const [tol, setTol] = useState("30");
  const [tol_unit, setTolUnit] = useState("ppm"); 
  const [database, setDatabase] = useState("kegg");
  const [num_rxn, setNumRxn] = useState("0_rxn");
  const [adduct, setAdduct] = useState("neutral");
  const navigate = useNavigate();


  // 1.2 OTHER FUNCTIONS 其余功能
  const [hasMassChanged, setHasMassChanged] = useState(false);

  const [loading, setLoading] = useState(false);


  // 1.3 ERROR HANDLE 输入错误检查
  const validateInput = () => {
    // Check if inputValue is empty
    if (!inputValue.trim()) {
      Modal.error({
        title: 'No Input',
        content: 'Please enter m/z values.',
      });
      return false;
    }

    // Check if tolerance is empty or not a number
    if (!tol.trim() || isNaN(parseFloat(tol))) {
      Modal.error({
        title: 'Invalid Tolerance',
        content: 'Tolerance must be a number.',
      });
      return false;
    }

    // Validate each m/z value
    const values = inputValue.split('\n').map(v => v.trim()).filter(v => v);
    for (let value of values) {
      if (isNaN(parseFloat(value))) {
        Modal.error({
          title: 'Invalid Input',
          content: 'Please check your m/z and tolerance, it should be float or int.',
        });
        return false;
      }
    }

    return true;
  };



  // 1.4 SUMBIT HANLE 提交
  const handleSubmit = async (e) => {
    e.preventDefault();
    // Call the validation function
    if (!validateInput()) {
      return; // Stop the submission if validation fails
    }
    setLoading(true);

    try {
      const response = await axios.post("http://111.229.151.228:5000/api/search", {
          input: inputValue,
          tol,
          database,
          tol_unit, 
          num_rxn,
          adduct
      });


      if (response.data) {
        navigate('/result', { state: { tableHTML: response.data.tableHTML } });
      }
    } catch (error) {
      console.error("There was an error sending the data!", error);
    }
  };

  


// ***** 2. LAYOUT 布局 ********
  return (
    <div className="search-container">
  
      {/* Left Part: Mass Label and Text Area */}
        <div className="form-section">
            {/* Mass */}
            <div className="mass-input-container">
                <div className="mass-label">Query m/z</div>
                <textarea 
                    rows={14}
                    placeholder="one feature (m/z) per line" 
                    className="mass-textarea"
                    style={{ color: hasMassChanged ? 'black' : 'gray' }} // Color changes based on hasMassChanged
                    defaultValue="one feature (m/z) per line
28.0
181.07389
114.0667" 
                    onChange={(e) => {
                        setInputValue(e.target.value);
                        setHasMassChanged(true);
                    }} 
                />
            </div>
        </div>
  


      {/* Right Part: Other Form Fields */}
      <div className="info-section">
          {/* Tolerance */}
          <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
            <br></br><br></br><br></br>
              <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Tolerance</div>
              <div style={{ flex: 2 }}>
                  <input type="text" placeholder="10" style={{ fontSize: '16px', width: '40%' }} onChange={(e) => setTol(e.target.value)} />
                  <select style={{ fontSize: '16px', marginLeft: '10px' }} onChange={(e) => setTolUnit(e.target.value)}>
                    <option value="ppm">ppm</option>
                    <option value="Da">Da</option>
                  </select>
              </div>
          </div>




          {/* Database Section */}
          <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
            <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Database</div>
            <div style={{ flex: 2 }}>
              <Radio.Group defaultValue="kegg" onChange={(e) => setDatabase(e.target.value)}>
                <Radio.Button value="kegg">MCIDxKEGG</Radio.Button>
                <Radio.Button value="plant">MCIDxPlant</Radio.Button>
              </Radio.Group>
            </div>
          </div>

          

          {/* Num. of Reactions */}
          <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
            <br></br>

              <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>
                  Num. of Reactions
                  <Tooltip title="The 1 rxn database is based on 76 metabolic reactions, while the 2 rxn extends predictions through two successive reaction steps.">
                      <QuestionCircleOutlined style={{ marginLeft: '5px', cursor: 'pointer' }} />
                  </Tooltip>
              </div>

              <div style={{ flex: 2 }}>
                  <Radio.Group className="custom-radio-color" defaultValue="0_rxn" onChange={(e) => setNumRxn(e.target.value)}>
                      <Radio.Button value="0_rxn">0 rxn</Radio.Button>
                      <Radio.Button value="1_rxn">1 rxn</Radio.Button>
                      <Radio.Button value="2_rxn">2 rxn</Radio.Button>
                  </Radio.Group>
              </div>

          </div>
  

          {/* Adduct Type */}
          <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
              <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Adduct Typ
                <Tooltip title="An adduct is a compound formed by the addition of an ion to a molecule and altering its mass.">
                  <QuestionCircleOutlined style={{ marginLeft: '5px', cursor: 'pointer' }} />
                </Tooltip>
              </div>
              <div style={{ flex: 2 }}>
                  <table>
                      <tbody>
                          <tr>
                              <td><input type="radio" name="adduct" value="neutral" defaultChecked onChange={(e) => setAdduct(e.target.value)} /> Neutral</td>
                              <td><input type="radio" name="adduct" value="m-h" onChange={(e) => setAdduct(e.target.value)} /> [M-H]-</td>
                          </tr>
                          <tr>
                              <td><input type="radio" name="adduct" value="m+h" onChange={(e) => setAdduct(e.target.value)} /> [M+H]+</td>
                              <td><input type="radio" name="adduct" value="m+na" onChange={(e) => setAdduct(e.target.value)} /> [M+Na]+</td>
                          </tr>
                          <tr>
                              <td><input type="radio" name="adduct" value="m+k" onChange={(e) => setAdduct(e.target.value)} /> [M+K]+</td>
                              <td><input type="radio" name="adduct" value="m+nh4" onChange={(e) => setAdduct(e.target.value)} /> [M+NH4]+</td>
                          </tr>
                      </tbody>
                  </table>
              </div>
          </div>


  
          {/* Submit Button */}
          <div style={{ display: 'flex', justifyContent: 'center', marginTop: '20px' }}>
          <br></br>
          <br></br>
          <br></br>
          <br></br>
              <Button 
                  type="primary" 
                  onClick={handleSubmit} 
                  style={{ 
                      backgroundColor: loading ? '#04AA6D' : '#00720a', 
                      borderColor: '#00720a', 
                      width: '30%', 
                      alignSelf: 'center'
                  }}
                  className={loading ? 'breathing' : ''}
              >
                  {loading ? <><i className="fa fa-spinner fa-spin"></i> Loading</> : 'Submit'}
              </Button>
          </div>
      </div>
    </div>
  );  
  }

  
export default Search;