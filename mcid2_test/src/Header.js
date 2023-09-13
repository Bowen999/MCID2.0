import React from 'react';
import { Link } from 'react-router-dom';
import { Menu } from 'antd';  // Assuming you're using Ant Design

const HeaderContent = () => {

  const menuItems = ['Home', 'MCIDxKEGG', 'MCIDxPlant', 'About', 'Tutorial', 'Contact Us'];
  const menuPaths = ['/', '/MCIDxKEGG', '/MCIDxPlant', '/about', '/tutorial', '/contact'];

  return (
    <header>
      <Menu mode="horizontal">
        {menuItems.map((item, index) => (
          <Menu.Item key={index}>
            <Link to={menuPaths[index]}>{item}</Link>
          </Menu.Item>
        ))}
      </Menu>
    </header>
  );
};

export default HeaderContent;
