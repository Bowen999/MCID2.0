import React from 'react';

const Intro = () => {
  return (
    <div>
      <h1>Introduction</h1>
          <p>"MyCompoundID 2.0 (MCIDxKEGG) is an evidence-based metabolome library.</p>
          <p>As an upgraded version of previous MyCompoundID, MCIDxKEGG predited theoretical metabolites based on metabolites from KEGG compound database via 76 biological reactions.</p>
          <p>More metabolites and their predicted products are included. And the 76 reactions have been optimized for more accurate prediction.</p>
          <p>In the current version, the zero-reaction database contains 11,164 metabolites. And the 1,811,882 predicted metabolites are enrolled in one-reaction database.</p>
          <p>The structure information is included to assist users in determining the unknown metabolites.</p>
          <p>Besides, the new MCID ID system can facilitate users to save and re-search predicted metabolites.</p>
          <p>Compared to using the standard library alone or the previous MCID, MyCompoundID 2.0 can provide better coverage.</p>
          <p>In the future, MCID 2.0 can be further expanded by recruiting more metabolites and more reactions.</p>
    </div>
  );
};

export default Intro;
