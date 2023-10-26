import React from 'react';
import './App.css';
import Contact from './Contact';
import Input from './Input';
import Result from './Result';
import Search from './Search';
import { BrowserRouter, Routes, Route, NavLink, useLocation } from 'react-router-dom';

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
          <NavLink exact to="/" activeClassName="active">
            Home
          </NavLink>
        </li>
        {/* Uncomment and modify as per your routes */}
        {/* 
        <li>
          <NavLink to="/news" activeClassName="active">
            News
          </NavLink>
        </li>
        <li>
          <NavLink to="/contact" activeClassName="active">
            Contact
          </NavLink>
        </li>
        */}
        <li>
          <NavLink to="/search" activeClassName="active">
            Search
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
        <Route path="/" element={<Home />} />
        {/* <Route path="/contact" element={<Contact />} /> */}
        {/* <Route path="/input" element={<Input />} /> */}
        <Route path="/result" element={<ResultWrapper />} />
        <Route path="/search" element={<Search />} />
      </Routes>
    </main>
  );
}

function Footer() {
  return (
    <footer>
      <p>Copyright &copy; 2013 University of Alberta</p>
    </footer>
  );
}

// function Home() {
//   return <div><br></br><br></br><br></br><br></br><br></br><br></br><h2>Test Demo of MCID 2.0</h2><h4>Click Search to try</h4><br></br><br></br><br></br></div>;
// }

// function Home() {
//   return (
//     <div>
//       <br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />
//       <h2>Test Demo of MCID 2.0</h2>
//       <h4>Click <span style={{ color: '#00720a', fontSize: '1.2em' }}>Search</span> to try</h4>
//       <br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />
//     </div>
//   );
// }
// function Home() {
//   return (
//     <div>
//       <br /><br /><br /><br /><br /><br />
//       <img src={`${process.env.PUBLIC_URL}/mcid.png`} alt="MCID Logo" style={{ display: 'block', margin: '0 auto' }} />
//       <h2>Test Demo of MCID 2.0</h2>
//       <h4>Click <span style={{ color: '#00720a', fontSize: '1.2em' }}>Search</span> to try</h4>
//       <br /><br /><br />
//     </div>
//   );
// }
function Home() {
  return (
    <div>
      <br /><br /><br /><br /><br /><br />
      <figure style={{ textAlign: 'center' }}>
        <img 
          src={`${process.env.PUBLIC_URL}/mcid.png`} 
          alt="MCID Logo" 
          style={{ display: 'block', margin: '0 auto', width: '250px' }} 
        />
      </figure>
      <br /><br />
      <h2>Test Demo of MCID 2.0</h2>
      <h4>Click <span style={{ color: '#00720a', fontSize: '1.2em' }}>Search</span> to try</h4>
      <br /><br /><br />
    </div>
  );
}




function ResultWrapper() {
  const location = useLocation();
  const tableHTML = location.state?.tableHTML;
  return <Result tableHTML={tableHTML} />;
}

export default App;

