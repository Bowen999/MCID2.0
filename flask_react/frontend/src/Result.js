// import React, { useState } from 'react';
// import { Table } from 'antd';
// import { DownloadOutlined } from '@ant-design/icons';
// import { Button } from 'antd';

// const baseWidth = 300;
// const columns = [
//   {
//       title: 'm/z',
//       dataIndex: 'mz',
//       width: baseWidth,
//   },
//   {
//       title: 'Num. of Hits',
//       dataIndex: 'numOfHits',
//       width: baseWidth*2,
//   },
//   {
//       title: 'Hits',
//       dataIndex: 'hits',
//       width: baseWidth*7,
//       render: (hits) => hits.map((hit, index) => (
//           <div key={index}>
//               <FoldableContent content={hit} />
//           </div>
//       ))
//   },
// ];


// function parseHtmlToData(html) {
//   const parser = new DOMParser();
//   const doc = parser.parseFromString(html, 'text/html');
//   const rows = Array.from(doc.querySelectorAll('table tr')).slice(1); // skipping header

//   return rows.map((row, index) => {
//     const cells = Array.from(row.querySelectorAll('td'));
//     return {
//       key: index, 
//       mz: cells[0].textContent,
//       numOfHits: cells[1].textContent,
//       hits: cells[2].innerHTML // Use innerHTML instead of textContent
//     };
//   });
// }

// function FoldableContent({ content }) {
//   const [isExpanded, setIsExpanded] = useState(false);

//   if (content.length <= 100) return <div dangerouslySetInnerHTML={{ __html: content }} />;

//   const displayContent = isExpanded 
//       ? content 
//       : content.slice(0, 100).split(' ').slice(0, -1).join(' ') + '...';

//   return (
//       <div>
//           <div dangerouslySetInnerHTML={{ __html: displayContent }} />
//           <a onClick={() => setIsExpanded(!isExpanded)}>
//               {isExpanded ? 'Show less' : 'Click to show more'}
//           </a>
//       </div>
//   );
// }




// function convertToCSV(data) {
//     const header = ['m/z', 'Num. of Hits', 'Hits'];
//     const csvRows = data.map(row => {
//       return `${row.mz},"${row.numOfHits}","${row.hits.join(', ')}"`;
//     });
//     return `${header.join(',')}\n${csvRows.join('\n')}`;
// }

  
// function downloadCSV(data) {
//     const csvContent = convertToCSV(data);
//     const blob = new Blob([csvContent], { type: 'text/csv' });
//     const url = URL.createObjectURL(blob);
//     const a = document.createElement('a');
//     a.setAttribute('hidden', '');
//     a.setAttribute('href', url);
//     a.setAttribute('download', 'selected_data.csv');
//     document.body.appendChild(a);
//     a.click();
//     document.body.removeChild(a);
// }

  

// const Result = ({ tableHTML }) => {
//     const [selectedRowKeys, setSelectedRowKeys] = useState([]);
  
//     // Parse the HTML into a suitable JSON structure
//     const dataSource = parseHtmlToData(tableHTML); 
//     const selectedRows = selectedRowKeys.map(key => dataSource[key]);
  
//     const onSelectChange = (newSelectedRowKeys) => {
//       setSelectedRowKeys(newSelectedRowKeys);
//     };
  
//     const rowSelection = {
//       selectedRowKeys,
//       onChange: onSelectChange,
//     };
  
 
//     return (
//       <div style={{ display: 'flex', flexDirection: 'column', height: '100%', backgroundColor: 'white' }}>
//         <Table rowSelection={rowSelection} columns={columns} dataSource={dataSource} />
//         <div style={{ marginTop: '20px', marginBottom: '20px', display: 'flex', justifyContent: 'center' }}>
//           <Button 
//             type="primary" 
//             icon={<DownloadOutlined />} 
//             onClick={() => downloadCSV(selectedRows)}
//             style={{ backgroundColor: '#00720a', borderColor: '#00720a' }} // Added style for the button color
//           >
//             Download Selected
//           </Button>
//         </div>
//       </div>
//     );
//   };

  
  

// // Parse HTML to an array of objects compatible with antd's Table
// function parseHtmlToData(html) {
//   const parser = new DOMParser();
//   const doc = parser.parseFromString(html, 'text/html');
//   const rows = Array.from(doc.querySelectorAll('table tr')).slice(1); // skipping header

//   return rows.map((row, index) => {
//     const cells = Array.from(row.querySelectorAll('td'));
//     return {
//       key: index, 
//       mz: cells[0].textContent,
//       numOfHits: cells[1].textContent,
//       hits: cells[2].textContent.split('\n') // split by \n to get a list of hits
//     };
//   });
// }

// export default Result;



// V 2
// import React, { useState } from 'react';
// import { Table, Button } from 'antd';
// import { DownloadOutlined } from '@ant-design/icons';

// const baseWidth = 300;
// const columns = [
//   {
//     title: 'm/z',
//     dataIndex: 'mz',
//     width: baseWidth,
//   },
//   {
//     title: 'Num. of Hits',
//     dataIndex: 'numOfHits',
//     width: baseWidth * 2,
//   },
//   {
//     title: 'Hits',
//     dataIndex: 'hits',
//     width: baseWidth * 7,
//     render: (hits) => <FoldableContent content={hits} />,
//   },
// ];


// function FoldableContent({ content }) {
//   const [isExpanded, setIsExpanded] = useState(false);
//   const hitElements = content.split(', ');
//   const shouldFold = hitElements.length > 10;
//   const displayedContent = shouldFold && !isExpanded ? hitElements.slice(0, 10).join(', ') + '...' : content;

//   return (
//     <div>
//       <div dangerouslySetInnerHTML={{ __html: displayedContent }} />
//       {shouldFold && (
//         <a onClick={() => setIsExpanded(!isExpanded)} style={{ cursor: 'pointer', color: '#00720a' }}>
//           {isExpanded ? 'Show less' : 'Show more'}
//         </a>
//       )}
//     </div>
//   );
// }

// function parseHtmlToData(html) {
//   const parser = new DOMParser();
//   const doc = parser.parseFromString(html, 'text/html');
//   const rows = Array.from(doc.querySelectorAll('table tr')).slice(1); // skipping header

//   return rows.map((row, index) => {
//     const cells = Array.from(row.querySelectorAll('td'));
//     return {
//       key: index,
//       mz: cells[0].textContent,
//       numOfHits: cells[1].textContent,
//       hits: cells[2].innerHTML, // Use innerHTML to preserve HTML tags
//     };
//   });
// }

// function convertToCSV(data) {
//   const header = ['m/z', 'Num. of Hits', 'Hits'];
//   const csvRows = data.map(row => {
//     return `${row.mz},"${row.numOfHits}","${row.hits}"`; // Keep the HTML tags in "Hits" column for CSV
//   });
//   return `${header.join(',')}\n${csvRows.join('\n')}`;
// }

// function downloadCSV(data) {
//   const csvContent = convertToCSV(data);
//   const blob = new Blob([csvContent], { type: 'text/csv' });
//   const url = URL.createObjectURL(blob);
//   const a = document.createElement('a');
//   a.setAttribute('hidden', '');
//   a.setAttribute('href', url);
//   a.setAttribute('download', 'selected_data.csv');
//   document.body.appendChild(a);
//   a.click();
//   document.body.removeChild(a);
// }

// const Result = ({ tableHTML }) => {
//   const [selectedRowKeys, setSelectedRowKeys] = useState([]);
  
//   // Parse the HTML into a suitable JSON structure
//   const dataSource = parseHtmlToData(tableHTML);
//   const selectedRows = selectedRowKeys.map(key => dataSource[key]);

//   const onSelectChange = (newSelectedRowKeys) => {
//     setSelectedRowKeys(newSelectedRowKeys);
//   };

//   const rowSelection = {
//     selectedRowKeys,
//     onChange: onSelectChange,
//   };

//   return (
//     <div style={{ display: 'flex', flexDirection: 'column', height: '100%', backgroundColor: 'white' }}>
//       <Table rowSelection={rowSelection} columns={columns} dataSource={dataSource} />
//       <div style={{ marginTop: '20px', marginBottom: '20px', display: 'flex', justifyContent: 'center' }}>
//         <Button
//           type="primary"
//           icon={<DownloadOutlined />}
//           onClick={() => downloadCSV(selectedRows)}
//           style={{ backgroundColor: '#00720a', borderColor: '#00720a' }} // Added style for the button color
//         >
//           Download Selected
//         </Button>
//       </div>
//     </div>
//   );
// };

// export default Result;





import React, { useState } from 'react';
import { Table, Button } from 'antd';
import { DownloadOutlined } from '@ant-design/icons';

const baseWidth = 300;

const FoldableContent = ({ content }) => {
  const [isExpanded, setIsExpanded] = useState(false);
  const hitElements = content.split(', ');
  const shouldFold = hitElements.length > 10;
  const displayContent = isExpanded || !shouldFold ? content : hitElements.slice(0, 10).join(', ') + '...';

  return (
    <div>
      <div dangerouslySetInnerHTML={{ __html: displayContent }} />
      {shouldFold && (
        <a href="#" onClick={() => window.open(content, '_blank')} style={{ cursor: 'pointer', color: '#00720a' }}>
          {isExpanded ? 'Show less' : 'Click to show more'}
        </a>
      )}
    </div>
  );
};


const columns = [
  {
    title: 'm/z',
    dataIndex: 'mz',
    width: baseWidth,
  },
  {
    title: 'Num. of Hits',
    dataIndex: 'numOfHits',
    width: baseWidth * 2,
  },
  {
    title: 'Hits',
    dataIndex: 'hits',
    width: baseWidth * 7,
    render: hits => <FoldableContent content={hits} />,
  },
];

function parseHtmlToData(html) {
  const parser = new DOMParser();
  const doc = parser.parseFromString(html, 'text/html');
  const rows = Array.from(doc.querySelectorAll('table tr')).slice(1); // skipping header

  return rows.map((row, index) => {
    const cells = Array.from(row.querySelectorAll('td'));
    return {
      key: index,
      mz: cells[0].textContent,
      numOfHits: cells[1].textContent,
      hits: cells[2].innerHTML, // Use innerHTML to preserve HTML tags
    };
  });
}

function convertToCSV(data) {
  const header = ['m/z', 'Num. of Hits', 'Hits'];
  const csvRows = data.map(row => {
    return `${row.mz},"${row.numOfHits}","${row.hits}"`;
  });
  return `${header.join(',')}\n${csvRows.join('\n')}`;
}

function downloadCSV(data) {
  const csvContent = convertToCSV(data);
  const blob = new Blob([csvContent], { type: 'text/csv' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.setAttribute('hidden', '');
  a.setAttribute('href', url);
  a.setAttribute('download', 'selected_data.csv');
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}

const Result = ({ tableHTML }) => {
  const [selectedRowKeys, setSelectedRowKeys] = useState([]);

  // Parse the HTML into a suitable JSON structure
  const dataSource = parseHtmlToData(tableHTML);
  const selectedRows = selectedRowKeys.map(key => dataSource.find(row => row.key === key));

  const onSelectChange = (newSelectedRowKeys) => {
    setSelectedRowKeys(newSelectedRowKeys);
  };

  const rowSelection = {
    selectedRowKeys,
    onChange: onSelectChange,
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', height: '100%', backgroundColor: 'white' }}>
      <Table rowSelection={rowSelection} columns={columns} dataSource={dataSource} />
      <div style={{ marginTop: '20px', marginBottom: '20px', display: 'flex', justifyContent: 'center' }}>
        <Button
          type="primary"
          icon={<DownloadOutlined />}
          onClick={() => downloadCSV(selectedRows)}
          style={{ backgroundColor: '#00720a', borderColor: '#00720a' }}
        >
          Download Selected
        </Button>
      </div>
    </div>
  );
};

export default Result;
