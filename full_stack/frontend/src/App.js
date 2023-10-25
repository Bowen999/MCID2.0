import React from 'react';
import './App.css';
import Contact from './Contact';
import Input from './Input';
import Result from './Result';
import Search from './Search';
import { BrowserRouter, Routes, Route, Link, Outlet, useLocation } from 'react-router-dom';

function App() {
  return (
    <BrowserRouter>
      <div className="App">
        <nav>
          <Link to="/">Home</Link>
          {/* <Link to="/contact">Contact</Link> */}
          {/* <Link to="/input">Input</Link> */}
          <Link to="/search">Search</Link>
        </nav>

        <Routes>
          <Route path="/" element={<Home />} />
          {/* <Route path="/contact" element={<Contact />} /> */}
          {/* <Route path="/input" element={<Input />} /> */}
          <Route path="/result" element={<ResultWrapper />} />
          <Route path="/search" element={<Search />} />
        </Routes>
      </div>
    </BrowserRouter>
  );
}

function Home() {
  return <div><br></br><br></br><br></br><h2>Test Demo of MCID 2.0</h2><h4>Click Search to try</h4></div>;
}

function ResultWrapper() {
  const location = useLocation();
  const tableHTML = location.state?.tableHTML;
  return <Result tableHTML={tableHTML} />;
}

export default App;

