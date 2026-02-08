import React from 'react';
import { BrowserRouter, Routes, Route, NavLink, useLocation, Navigate } from 'react-router-dom';
import './App.css';
import Home from './Home'; // Import Home component
import Result from './Result';
import Search from './Search';
import FAQ from './FAQ';
import Contact from './Contact';

function App() {
  return (
    <BrowserRouter>
      <div className="App">
        <Header />
        <Content />
        <Footer />
      </div>
    </BrowserRouter>
  );
}

function Header() {
  return (
    <header style={{ backgroundColor: '#ffffff' }}>
      <ul>
        <li>
          <NavLink exact to="/home" activeClassName="active">
            Home
          </NavLink>
        </li>
        <li>
          <NavLink to="/search" activeClassName="active">
            Search
          </NavLink>
        </li>
        <li>
          <NavLink to="/faq" activeClassName="active">
            FAQ
          </NavLink>
        </li>
        <li>
          <NavLink to="/contact" activeClassName="active">
            Contact us
          </NavLink>
        </li>
      </ul>
    </header>
  );
}

function Content() {
  return (
    <main className="main-content">
      <Routes>
        <Route path="/" element={<Navigate to="/home" />} /> {/* Redirect to /home */}
        <Route path="/home" element={<Home />} /> {/* Home Route */}
        <Route path="/result" element={<ResultWrapper />} />
        <Route path="/search" element={<Search />} />
        <Route path="/faq" element={<FAQ />} />
        <Route path="/contact" element={<Contact />} />
      </Routes>
    </main>
  );
}

function Footer() {
  return (
    <footer style={{ color: 'gray' }}>
      <p>Copyright &copy; 2023 University of Alberta</p>
    </footer>
  );
}

function ResultWrapper() {
  const location = useLocation();
  const tableHTML = location.state?.tableHTML;
  return <Result tableHTML={tableHTML} />;
}


export default App;
