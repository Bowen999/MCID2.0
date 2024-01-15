// frontend/src/Contact.js

// Contact.js

import React, { useState, useEffect } from 'react';
import './App.css';

function Contact() {
    const [isSmallScreen, setIsSmallScreen] = useState(window.innerWidth < 600);

    useEffect(() => {
        const checkScreenSize = () => {
            setIsSmallScreen(window.innerWidth < 600);
        };

        window.addEventListener('resize', checkScreenSize);

        return () => {
            window.removeEventListener('resize', checkScreenSize);
        };
    }, []);

    const containerClass = isSmallScreen ? 'contact-container-small' : 'contact-container';
    const imgClass = isSmallScreen ? 'small-screen' : 'large-screen';


    return (
        <div className={containerClass}>
            <div className="contact-info">
                <h2 id="contact-us" className="InfoTitle">Contact us</h2>
                <p><strong>Prof. Liang Li</strong></p>
                <p><strong>Department</strong>: Department of Chemistry, University of Alberta</p>
                <p><strong>Address</strong>: 11227 Saskatchewan Drive, Edmonton, Alberta T6G 2G2, Canada</p>
                <p><strong>E-mail</strong>: liang.li@ualberta.ca</p>
                <br></br>
                <h2 id="feedback" className="InfoTitle">Feedback</h2>
                <p>Lorem ipsum dolor sit amet, ut vix exerci pertinacia...</p>
            </div>
            <div className="contact-image">
                <img 
                    src="https://liweb.chem.ualberta.ca/wp-content/uploads/sites/73/2019/05/Chemistry-2019.jpg" 
                    alt="Chemistry Department" 
                    className={imgClass}
                />
            </div>
        </div>
    );    
}

export default Contact;


