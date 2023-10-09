import React from 'react';
import Sidebar from './Sidebar';

const reactions = [
    { id: 1, reaction: '-H2', massDifference: -2.015650, description: 'dehydrogenation' },
    { id: 2, reaction: '+H2', massDifference: 2.015650, description: 'hydrogenation' },
    { id: 3, reaction: '-CH2', massDifference: -14.015650, description: 'demethylation' },
    { id: 4, reaction: '+CH2', massDifference: 14.015650, description: 'methylation' },
    { id: 5, reaction: '-NH', massDifference: -15.010899, description: 'loss of NH' },
    { id: 6, reaction: '+NH', massDifference: 15.010899, description: 'addition of NH' },
    { id: 7, reaction: '-O', massDifference: -15.994915, description: 'loss of oxygen' },
    { id: 8, reaction: '+O', massDifference: 15.994915, description: 'oxidation' },
    { id: 9, reaction: '-NH3', massDifference: -17.026549, description: 'loss of ammonia' },
    { id: 10, reaction: '+NH3', massDifference: 17.026549, description: 'addition of ammonia' },
    { id: 11, reaction: '-H2O', massDifference: -18.010565, description: 'loss of water' },
    { id: 12, reaction: '+H2O', massDifference: 18.010565, description: 'addition of water' },
    { id: 13, reaction: '-CO', massDifference: -27.994915, description: 'loss of CO' },
    { id: 14, reaction: '+CO', massDifference: 27.994915, description: 'addition of CO' },
    { id: 15, reaction: '-C2H4', massDifference: -28.031300, description: 'loss of C2H4' },
    { id: 16, reaction: '+C2H4', massDifference: 28.031300, description: 'addition of C2H4' },
    { id: 17, reaction: '-C2H2O', massDifference: -42.010565, description: 'deacetylation' },
    { id: 18, reaction: '+C2H2O', massDifference: 42.010565, description: 'acetylation' },
    { id: 19, reaction: '-CO2', massDifference: -43.989830, description: 'loss of CO2' },
    { id: 20, reaction: '+CO2', massDifference: 43.989830, description: 'addition of CO2' },
    { id: 21, reaction: 'SO3H->SH', massDifference: -47.984745, description: 'sulfonic acid to thiol' },
    { id: 22, reaction: 'SH->SO3H', massDifference: 47.984745, description: 'thiol to sulfonic acid' },
    { id: 23, reaction: '-C2H3NO', massDifference: -57.021464, description: 'loss of glycine' },
    { id: 24, reaction: '+C2H3NO', massDifference: 57.021464, description: 'glycine conjugation' },
    { id: 25, reaction: '-SO3', massDifference: -79.956817, description: 'loss of sulfate' },
    { id: 26, reaction: '+SO3', massDifference: 79.956817, description: 'sulfate conjugation' },
    { id: 27, reaction: '-HPO3', massDifference: -79.966333, description: 'loss of phosphate' },
    { id: 28, reaction: '+HPO3', massDifference: 79.966333, description: 'addition of phosphate' },
    { id: 29, reaction: '-C4H3N3', massDifference: -93.032697, description: 'loss of cytosine' },
    { id: 30, reaction: '+C4H3N3', massDifference: 93.032697, description: 'addition of cytosine' },
    { id: 31, reaction: '-C4H2N2O', massDifference: -94.016713, description: 'loss of uracil' },
    { id: 32, reaction: '+C4H2N2O', massDifference: 94.016713, description: 'addition of uracil' },
    { id: 33, reaction: '-C3H5NOS', massDifference: -103.009186, description: 'loss of cysteine' },
    { id: 34, reaction: '+C3H5NOS', massDifference: 103.009186, description: 'cysteine conjugation' },
    { id: 35, reaction: '-C2H5NO2S', massDifference: -107.004101, description: 'loss of taurine' },
    { id: 36, reaction: '+C2H5NO2S', massDifference: 107.004101, description: 'taurine conjugation' },
    { id: 37, reaction: '-C5H4N2O', massDifference: -108.032363, description: 'loss of thymine' },
    { id: 38, reaction: '+C5H4N2O', massDifference: 108.032363, description: 'addition of thymine' },
    { id: 39, reaction: '- (C5H5N5 - H2O)', massDifference: -117.043930, description: 'loss of adenine' },
    { id: 40, reaction: '+ (C5H5N5 - H2O)', massDifference: 117.043930, description: 'addition of adenine' },
    { id: 41, reaction: '-C3H5NO2S', massDifference: -119.004101, description: 'loss of S-cysteine' },
    { id: 42, reaction: '+C3H5NO2S', massDifference: 119.004101, description: 'S-cysteine conjugation' },
    { id: 43, reaction: '-C5H8O4', massDifference: -132.042260, description: 'loss of D-ribose' },
    { id: 44, reaction: '+C5H8O4', massDifference: 132.042260, description: 'addition of D-ribose' },
    { id: 45, reaction: '-C5H3N5', massDifference: -133.038845, description: 'loss of guanine' },
    { id: 46, reaction: '+C5H3N5', massDifference: 133.038845, description: 'addition of guanine' },
    { id: 47, reaction: '-C5H8O6', massDifference: -176.032089, description: 'loss of glucuronic acid' },
    { id: 48, reaction: '+C5H8O6', massDifference: 176.032089, description: 'addition of glucuronic acid' },
    { id: 49, reaction: '-C6H10O6', massDifference: -182.052430, description: 'loss of glucose' },
    { id: 50, reaction: '+C6H10O6', massDifference: 182.052430, description: 'addition of glucose' },
    { id: 51, reaction: '-C5H4N2O3', massDifference: -156.032363, description: 'loss of uridine' },
    { id: 52, reaction: '+C5H4N2O3', massDifference: 156.032363, description: 'addition of uridine' },
    { id: 53, reaction: '-C5H7N5O3', massDifference: -195.075365, description: 'loss of adenosine' },
    { id: 54, reaction: '+C5H7N5O3', massDifference: 195.075365, description: 'addition of adenosine' },
    { id: 55, reaction: '-C5H5N5O3', massDifference: -197.071006, description: 'loss of guanosine' },
    { id: 56, reaction: '+C5H5N5O3', massDifference: 197.071006, description: 'addition of guanosine' },
    { id: 57, reaction: '-C5H6N2O2', massDifference: -126.032716, description: 'loss of cytidine' },
    { id: 58, reaction: '+C5H6N2O2', massDifference: 126.032716, description: 'addition of cytidine' },
    { id: 59, reaction: '-C5H5N5O', massDifference: -151.039595, description: 'loss of deoxyadenosine' },
    { id: 60, reaction: '+C5H5N5O', massDifference: 151.039595, description: 'addition of deoxyadenosine' },
    { id: 61, reaction: '-C5H3N5O', massDifference: -149.035236, description: 'loss of deoxyguanosine' },
    { id: 62, reaction: '+C5H3N5O', massDifference: 149.035236, description: 'addition of deoxyguanosine' },
    { id: 63, reaction: '-C5H4N2O', massDifference: -110.032364, description: 'loss of deoxyuridine' },
    { id: 64, reaction: '+C5H4N2O', massDifference: 110.032364, description: 'addition of deoxyuridine' },
    { id: 65, reaction: '-C5H5N2O', massDifference: 111.036728, description: 'loss of deoxycytidine' },
    { id: 66, reaction: '+C5H5N2O', massDifference: 111.036728, description: 'addition of deoxycytidine' },
    { id: 67, reaction: '-C10H15N5O10P2', massDifference: -427.052500, description: 'loss of ADP' },
    { id: 68, reaction: '+C10H15N5O10P2', massDifference: 427.052500, description: 'addition of ADP' },
    { id: 69, reaction: '-C10H16N5O13P3', massDifference: -507.018150, description: 'loss of ATP' },
    { id: 70, reaction: '+C10H16N5O13P3', massDifference: 507.018150, description: 'addition of ATP' },
    { id: 71, reaction: '-C10H14N5O7P', massDifference: -347.057910, description: 'loss of AMP' },
    { id: 72, reaction: '+C10H14N5O7P', massDifference: 347.057910, description: 'addition of AMP' },
    { id: 73, reaction: '-C5H9NO4', massDifference: -147.053160, description: 'loss of glutamate' },
    { id: 74, reaction: '+C5H9NO4', massDifference: 147.053160, description: 'addition of glutamate' },
    { id: 75, reaction: '-C5H9NO3', massDifference: -131.037480, description: 'loss of glutamine' },
    { id: 76, reaction: '+C5H9NO3', massDifference: 131.037480, description: 'addition of glutamine' }



    
  // ... (Continue for all other rows)
];

const Rxn = () => {
  return (
    <table style={{ width: '100%', borderCollapse: 'collapse', margin: '20px 0', fontSize: '16px', textAlign: 'left' }}>
      <thead>
        <tr>
          <th style={{ padding: '10px', borderBottom: '1px solid #ddd', backgroundColor: '#f2f2f2' }}> </th>
          <th style={{ padding: '10px', borderBottom: '1px solid #ddd', backgroundColor: '#f2f2f2' }}>Reaction</th>
          <th style={{ padding: '10px', borderBottom: '1px solid #ddd', backgroundColor: '#f2f2f2' }}>Mass Difference (Da)</th>
          <th style={{ padding: '10px', borderBottom: '1px solid #ddd', backgroundColor: '#f2f2f2' }}>Description</th>
        </tr>
      </thead>
      <tbody>
        {reactions.map(rxn => (
          <tr key={rxn.id} style={{ borderBottom: '1px solid #ddd', hover: { backgroundColor: '#f5f5f5' } }}>
            <td style={{ padding: '10px' }}>{rxn.id}</td>
            <td style={{ padding: '10px' }}>{rxn.reaction}</td>
            <td style={{ padding: '10px' }}>{rxn.massDifference}</td>
            <td style={{ padding: '10px' }}>{rxn.description}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
};

const Rxn_page = () => {
  return (
    <div style={{ display: 'flex' }}>
      
      {/* Sidebar */}
      <div style={{ flex: '1', fontSize: '18px' }}>
        <Sidebar />
      </div>
      
      {/* Content */}
      <div style={{ flex: '4', padding: '0 20px' }}>
        <h1>CIL Metabolomics</h1>

        {/* Table */}
        <Rxn />
      </div>
    </div>
  );
};

export default Rxn_page;

