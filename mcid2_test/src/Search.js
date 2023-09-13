import React from 'react';
import { Radio, Button } from 'antd';
import type { RadioChangeEvent } from 'antd';

const onChange = (e: RadioChangeEvent) => {
  console.log(`radio checked:${e.target.value}`);
};

const SearchContent = () => {
  return (
    <div className="site-layout-content" style={{ display: 'flex' }}>
      
      {/* Left Part */}
      <div style={{ flex: 3, paddingRight: '20px', display: 'flex', flexDirection: 'column', justifyContent: 'space-between' }}>
        
        {/* Mass */}
        <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
          <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Mass</div>
          <div style={{ flex: 2 }}>
            <input type="text" placeholder="Enter m/z" style={{ fontSize: '16px', width: '60%' }} />
          </div>
        </div>
        
        {/* Tolerance */}
        <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
          <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Tolerance</div>
          <div style={{ flex: 2 }}>
            <input type="text" placeholder="Enter tolerance" style={{ fontSize: '16px', width: '40%' }} />
            <select style={{ fontSize: '16px', marginLeft: '10px' }}>
              <option value="Da">Da</option>
              <option value="ppm">ppm</option>
            </select>
          </div>
        </div>
        
        {/* Num. of Reactions */}
        <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
          <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Num. of Reactions</div>
          <div style={{ flex: 2 }}>
            <Radio.Group className="custom-radio-color" onChange={onChange} defaultValue="no_reaction">
              <Radio.Button value="no_reaction">No reaction</Radio.Button>
              <Radio.Button value="1_reaction">1 reaction</Radio.Button>
              <Radio.Button value="2_reaction">2 reaction</Radio.Button>
            </Radio.Group>
          </div>
        </div>
        
        {/* Adduct Type */}
        <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
          <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Adduct Type</div>
          <div style={{ flex: 2 }}>
            <input type="radio" name="adduct" value="neutral" /> Neutral
            <br />
            <input type="radio" name="adduct" value="m+h" /> [M+H]+
            <input type="radio" name="adduct" value="m+na" /> [M+Na]+
            <br />
            <input type="radio" name="adduct" value="m-h" /> [M-H]-
            <br />
            <br />
            <br />
          </div>
        </div>
        
        {/* Button */}
        {/* put search button at center */}
        <Button type="primary" style={{ backgroundColor: '#00720a', borderColor: '#00720a', width: '30%', alignSelf: 'center' }}>
          Search
        </Button>
        <br />
        <br />
        <br />            
      </div>
      
      {/* Right Part */}
      <div style={{ flex: 2, paddingLeft: '20px' }}>
        <h2>LC-MS Metabolomics</h2>
        <p>The MS Search program allows a user to search a query mass to generate a list of possible matches with the metabolites in an evidence-based metabolome library (EML). This library is composed of 8,021 known human endogenous metabolites and their predicted metabolic products (375,809 compounds from one metabolic reaction and 10,583,901 from two reactions). If the MS/MS spectrum of the query mass ion is available, this program also allows the user to interpret the MS/MS spectral pattern against the list of metabolite candidates generated from the mass search to narrow down the list into one or a few unique structures.</p>
      </div>
      
    </div>
  );
};

export default SearchContent;
