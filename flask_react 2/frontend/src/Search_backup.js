import React, { useState } from 'react';
import axios from 'axios';
import { useNavigate } from 'react-router-dom';
import { Input, Radio, Button } from 'antd';
const { TextArea } = Input;
// import { Button } from 'antd';



function Search() {
  const [inputValue, setInputValue] = useState("504\n534.99\n28.0\n549");
  const [tol, setTol] = useState("0.01");
  const [tol_unit, setTolUnit] = useState("Da"); // Step 1: Initialize the tol_unit state
  const [num_rxn, setNumRxn] = useState("1_rxn");
  const [adduct, setAdduct] = useState("neutral");
  const navigate = useNavigate();

  const [hasMassChanged, setHasMassChanged] = useState(false);

  const [loading, setLoading] = useState(false);


  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);

    try {
      const response = await axios.post("http://111.229.151.228:5000/api/search", {
          input: inputValue,
          tol,
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

  



  return (
    <div className="search-container" style={{ display: 'flex', justifyContent: 'space-between', padding: '40px', textAlign: 'left' }}>
    
        
        {/* Left Part: Form */}
        {/* <div style={{ flex: 2 }}> */}
        <div className="form-section">
            <form onSubmit={handleSubmit}>
                

                {/* Mass */}
                <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
                    <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Mass</div>
                    <div style={{ flex: 2 }}>
                        <textarea 
                            rows={4} 
                            placeholder="one feature (m/z) per line" 
                            style={{ 
                                fontSize: '16px', 
                                width: '60%',
                                color: hasMassChanged ? 'black' : 'gray' 
                            }} 
                            defaultValue="504
534.99
28.0
549" 
                            onChange={(e) => {
                                setInputValue(e.target.value);
                                setHasMassChanged(true);
                            }} 
                        />
                    </div>
                </div>


                
                {/* Tolerance */}
                <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
                    <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Tolerance</div>
                    <div style={{ flex: 2 }}>
                        {/* default is 0.01 Da */}
                        <input type="text" placeholder="0.01" style={{ fontSize: '16px', width: '40%' }} onChange={(e) => setTol(e.target.value)} />
                        <select style={{ fontSize: '16px', marginLeft: '10px' }} onChange={(e) => setTolUnit(e.target.value)}>
                            <option value="Da">Da</option>
                            <option value="ppm">ppm</option>
                        </select>
                    </div>
                </div>
                
                {/* Num. of Reactions */}
                <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
                    <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Num. of Reactions</div>
                    <div style={{ flex: 2 }}>
                        {/* default is 1 reaction */}
                        <Radio.Group className="custom-radio-color" defaultValue="1_rxn" onChange={(e) => setNumRxn(e.target.value)}>
                            <Radio.Button value="0_rxn">0 reaction</Radio.Button>
                            <Radio.Button value="1_rxn">1 reaction</Radio.Button>
                            <Radio.Button value="2_rxn">2 reactions</Radio.Button>
                        </Radio.Group>
                    </div>
                </div>

                {/* Adduct Type */}
                <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
                    <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Adduct Type</div>
                    <div style={{ flex: 2 }}>
                        <br />
                        {/* <input type="radio" name="adduct" value="neutral" onChange={(e) => setAdduct(e.target.value)} /> Neutral */}
                        <input type="radio" name="adduct" value="neutral" defaultChecked onChange={(e) => setAdduct(e.target.value)} /> Neutral
                        <input type="radio" name="adduct" value="m-h" onChange={(e) => setAdduct(e.target.value)} /> [M-H]-
                        <br />
                        <input type="radio" name="adduct" value="m+h" onChange={(e) => setAdduct(e.target.value)} /> [M+H]+
                        <input type="radio" name="adduct" value="m+na" onChange={(e) => setAdduct(e.target.value)} /> [M+Na]+
                        <input type="radio" name="adduct" value="m+k" onChange={(e) => setAdduct(e.target.value)} /> [M+K]+
                        <input type="radio" name="adduct" value="m+nh4" onChange={(e) => setAdduct(e.target.value)} /> [M+NH4]+
                    </div>
                </div>
                


                {/* Submit Button */}
                <div style={{ display: 'flex', justifyContent: 'center', marginTop: '20px' }}>
                    <br></br>
                    {/* <Button 
                        type="primary" 
                        onClick={handleSubmit} 
                        style={{ backgroundColor: '#00720a', borderColor: '#00720a', width: '30%', alignSelf: 'center' }}>
                        Submit
                    </Button> */}
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
            </form>
        </div>
        
        {/* Right Part: MCID 2.0 Information */}
        {/* <div style={{ flex: 2, paddingLeft: '20px', textAlign: 'left' }}> */}
        <div className="info-section">
            <h2 id="mcid-2-0">MCID 2.0</h2>
            <p>MCID 2.0 (MCIDxKEGG) is an evidence-based metabolome library, by merging metabolites from the <a href="https://www.genome.jp/kegg/">KEGG</a>, and generating predicted theoretical metabolites derived from biological reactions.</p>
            <p>We&#39;ve harnessed 11,164 distinct metabolites to create a zero-reaction database. Furthermore, with the incorporation of <a href="./about/reaction">76 biological reactions</a>, we&#39;ve formed both a one-reaction and a two-reaction database based on the number of reactions</p>
            <p>Queries can be conducted using m/z values from LC-MS or <a href="./about/CIL">CIL LC-MS</a>, or directly through the MCID.</p>
        </div>

    </div>
  );
  }

  
export default Search;