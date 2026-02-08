// frontend/src/Input.js

import React, { useState } from 'react';
import { useNavigate } from 'react-router-dom';

function Input() {
    const [firstValue, setFirstValue] = useState('');
    const [secondValue, setSecondValue] = useState('');
    const navigate = useNavigate();

    const handleSubmit = async (e) => {
        e.preventDefault();

        try {
            const response = await fetch('/api/test', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    firstValue: parseInt(firstValue),
                    secondValue: parseInt(secondValue)
                })
            });

            if (!response.ok) {
                throw new Error('Network response was not ok');
            }

            const data = await response.text();
            
            // Navigate to the result page and pass the data along
            navigate('/result', { state: { tableHTML: data } });

        } catch (error) {
            console.error("There was a problem with the fetch operation:", error.message);
        }
    };

    return (
        <div>
            <h2>Input Page</h2>
            <form onSubmit={handleSubmit}>
                <input 
                    type="number" 
                    placeholder="Enter first value" 
                    value={firstValue} 
                    onChange={e => setFirstValue(e.target.value)}
                />
                <input 
                    type="number" 
                    placeholder="Enter second value" 
                    value={secondValue} 
                    onChange={e => setSecondValue(e.target.value)}
                />
                <button type="submit">Submit</button>
            </form>
        </div>
    );
}

export default Input;