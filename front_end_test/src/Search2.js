import React from 'react';
import { Radio, Button, Input } from 'antd';
import type { RadioChangeEvent } from 'antd';

const { TextArea } = Input; // Extract TextArea from Input

const onChange = (e: RadioChangeEvent) => {
  console.log(`radio checked:${e.target.value}`);
};


const SearchContent2 = () => {
  return (
    <div className="site-layout-content" style={{ display: 'flex' }}>
      
      {/* Left Part */}
      <div style={{ flex: 3, paddingRight: '20px', display: 'flex', flexDirection: 'column', justifyContent: 'space-between' }}>
        
        {/* Mass */}
         {/* <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
          <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Mass</div>
          <div style={{ flex: 2 }}>
            <input type="text" placeholder="Enter m/z" style={{ fontSize: '16px', width: '60%' }} />
          </div>
        </div> */}

        {/* Mass */}
        <div style={{ display: 'flex', alignItems: 'center', marginBottom: '20px' }}>
          <div style={{ flex: 1, fontWeight: 'bold', color: '#00720a', fontSize: '16px' }}>Mass</div>
          <div style={{ flex: 2 }}>
            <TextArea rows={4} placeholder="one feature (m/z) per line" style={{ fontSize: '16px', width: '60%' }} /> {/* Use TextArea here */}
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
            <br />
            <input type="radio" name="adduct" value="neutral" /> Neutral
            <input type="radio" name="adduct" value="m-h" /> [M-H]-
            <br />
            <input type="radio" name="adduct" value="m+h" /> [M+H]+
            <input type="radio" name="adduct" value="m+na" /> [M+Na]+
            <input type="radio" name="adduct" value="m+k" /> [M+K]+
            <input type="radio" name="adduct" value="m+nh4" /> [M+NH4]+
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
      <h2 id="mcid-2-0">MCIDxPlant</h2>
      <p>MCID 2.0 (MCIDxKEGG) is an evidence-based metabolome library, by merging metabolites from the <a href="https://www.genome.jp/kegg/">KEGG</a>, and generating predicted theoretical metabolites derived from biological reactions.</p>
      <p>We&#39;ve harnessed 11,164 distinct metabolites to create zero-reaction database. Furthermore, with the incorporation of <a href="./about/reaction">76 biological reactions</a>, we&#39;ve formed both a one-reaction and a two-reaction database based on the number of reactions  </p>
      <p>Queries can be conducted using m/z values from LC-MS or <a href="./about/CIL">CIL LC-MS</a>, or directly through the MCID.</p>
      </div>
      
    </div>
  );
};

export default SearchContent2;
