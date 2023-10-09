// Sidebar.js
import React from 'react';
import { Link } from 'react-router-dom'; // Assuming you are using react-router

const Sidebar = () => {
  return (
    <div className="sidebar">
      <ul>
        <li><Link to="/about">Introduction</Link></li>
        <li><Link to="/about/CIL">CIL Metabolomics</Link></li>
        <li><Link to="/about/reaction">Possible Reactions</Link></li>
        <li><Link to="/about/services">Services</Link></li>
        <li><Link to="/about/FAQ">FAQ</Link></li>
        <li><Link to="/about/cite">How to cite</Link></li>
      </ul>
    </div>
  );
}

export default Sidebar;
