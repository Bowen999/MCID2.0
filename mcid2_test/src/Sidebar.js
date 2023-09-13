// Sidebar.js
import React from 'react';
import { Link } from 'react-router-dom'; // Assuming you are using react-router

const Sidebar = () => {
  return (
    <div className="sidebar">
      <ul>
        <li><Link to="/Intro">Introduction</Link></li>
        <li><Link to="/About/CIL">CIL</Link></li>
        <li><Link to="/About/reaction">Possible Reactions</Link></li>
        <li><Link to="/About/serverse">Serverse</Link></li>
        <li><Link to="/About/FAQ">FAQ</Link></li>
        <li><Link to="/About/cite">How to cite</Link></li>
      </ul>
    </div>
  );
}

export default Sidebar;
