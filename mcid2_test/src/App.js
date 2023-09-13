import React from 'react';
import { Layout, Menu, Breadcrumb } from 'antd';
import {
  BrowserRouter as Router,
  Routes,
  Route,
  Link,
} from 'react-router-dom';
import Footer from './Footer';
import Intro from './Intro';
import Search from './Search';
import Search2 from './Search2';  // Assuming MCIDxPlant is mapped to Search2.js
import Tutorial from './Tutorial';
import Contact from './Contact';
import Home from './Home';

const { Header, Content } = Layout;

const App: React.FC = () => {
  return (
    <Router>
      <Layout>
        <Header
          style={{
            position: 'sticky',
            top: 0,
            zIndex: 1,
            width: '100%',
            display: 'flex',
            alignItems: 'center',
          }}
        >
          <div className="demo-logo" />
          <Menu
            theme="light"
            mode="horizontal"
            defaultSelectedKeys={['1']}
          >
            <Menu.Item key="1">
              <Link to="/">Home</Link>
            </Menu.Item>
            <Menu.Item key="2">
              <Link to="/MCIDxKEGG">MCIDxKEGG</Link>
            </Menu.Item>
            <Menu.Item key="3">
              <Link to="/MCIDxPlant">MCIDxPlant</Link>
            </Menu.Item>
            <Menu.Item key="4">
              <Link to="/about">About</Link>
            </Menu.Item>
            <Menu.Item key="5">
              <Link to="/tutorial">Tutorial</Link>
            </Menu.Item>
            <Menu.Item key="6">
              <Link to="/contact">Contact Us</Link>
            </Menu.Item>
          </Menu>
        </Header>
        <Content className="site-layout" style={{ padding: '0 50px' }}>
          <Breadcrumb style={{ margin: '16px 0' }}>
            <Breadcrumb.Item>
              <Link to="/">Home</Link>
            </Breadcrumb.Item>
            {/* ... other breadcrumb items as needed ... */}
          </Breadcrumb>
          <div style={{ padding: 24, minHeight: 380 }}>
            <Routes>
              <Route path="/" element={<Home />} />
              <Route path="/MCIDxKEGG" element={<Search />} />
              <Route path="/MCIDxPlant" element={<Search2 />} />
              <Route path="/about" element={<Intro />} />
              <Route path="/tutorial" element={<Tutorial />} />
              <Route path="/contact" element={<Contact />} />
              {/* ... other routes as needed ... */}
            </Routes>
          </div>
        </Content>
        <Footer />
      </Layout>
    </Router>
  );
};

export default App;
