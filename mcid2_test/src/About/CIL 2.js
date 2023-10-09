import React from 'react';
import Sidebar from './Sidebar';

const CIL = () => {
  return (
    <div style={{ display: 'flex' }}>
      
      {/* Sidebar */}

      <div style={{ flex: '1', fontSize: '18px' }}>
      <br />
        <Sidebar />
      </div>
      
      {/* Content */}
      <div style={{ flex: '4', padding: '0 20px' }}>
        <h1>CIL Metabolomics</h1>
        <p>To increase metabolome coverage and achieve accurate quantification for all detectable metabolites, we use a Chemical Isotope Labeling (CIL) Metabolomics Platform (Zhao et al., 2019, Zhao et al., 2020) for metabolome analysis. The whole metabolome is analyzed by combining the analysis of four submetabolomes: amine/phenol, carboxyl, carbonyl and hydroxyl submetabolome. The combined results from four channels are able to cover 85% to 95% of the entire chemical space of the metabolome. However, a user can choose the labelling method(s) to analyze one or more chemical-group-based submetabolomes, depending on budget and desired coverage. Using one submetabolome analysis, we detect and perform accurate relative quantification of up to 2,500 metabolites (not features) in many different types of samples, including urine, blood, plant extract and cell extracts. Using a combination of four submetabolome analyses, we can cover over 7,000 to 10,000 metabolites (sample dependent) per sample.</p>
      </div>

    </div>
  );
};

export default CIL;
