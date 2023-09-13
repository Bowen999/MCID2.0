import React from 'react';
import Sidebar from './Sidebar';

const FAQ = () => {
  return (
    <div style={{ display: 'flex' }}>
      
      {/* Sidebar */}

      <div style={{ flex: '1', fontSize: '18px' }}>
      <br />
        <Sidebar />
      </div>
      
      {/* Content */}
      <div style={{ flex: '4', padding: '0 20px' }}>
        <h1>FAQ</h1>
        <p>....</p>
      </div>

    </div>
  );
};

export default FAQ;
