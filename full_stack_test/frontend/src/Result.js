import React from 'react';
function Result(props) {
    return (
        <div>
            <h2>Result Page</h2>
            <div dangerouslySetInnerHTML={{ __html: props.tableHTML }} />
        </div>
    );
}

export default Result;



