import React, { useState } from 'react';
import { Table } from 'antd';
import { DownloadOutlined } from '@ant-design/icons';
import { Button } from 'antd';

const baseWidth = 300;
// const columns = [
//   {
//     title: 'm/z',
//     dataIndex: 'mz',
//     width: baseWidth,
//   },
//   {
//     title: 'Num. of Hits',
//     dataIndex: 'numOfHits',
//     width: baseWidth*2,
//   },
//   {
//     title: 'Hits',
//     dataIndex: 'hits',
//     width: baseWidth*7,
//     render: (hits) => hits.map((hit, index) => <div key={index} dangerouslySetInnerHTML={{ __html: hit }} />)

//   },
// ];
const columns = [
  {
      title: 'm/z',
      dataIndex: 'mz',
      width: baseWidth,
  },
  {
      title: 'Num. of Hits',
      dataIndex: 'numOfHits',
      width: baseWidth*2,
  },
  {
      title: 'Hits',
      dataIndex: 'hits',
      width: baseWidth*7,
      render: (hits) => hits.map((hit, index) => (
          <div key={index}>
              <FoldableContent content={hit} />
          </div>
      ))
  },
];

function FoldableContent({ content }) {
  const [isExpanded, setIsExpanded] = useState(false);

  if (content.length <= 100) return <div dangerouslySetInnerHTML={{ __html: content }} />;

  const displayContent = isExpanded 
      ? content 
      : content.slice(0, 100).split(' ').slice(0, -1).join(' ') + '...'; // This ensures that the content doesn't break at a single word.

  return (
      <div>
          <div dangerouslySetInnerHTML={{ __html: displayContent }} />
          <a onClick={() => setIsExpanded(!isExpanded)}>
              {isExpanded ? 'Show less' : 'Click to show more'}
          </a>
      </div>
  );
}


function convertToCSV(data) {
    const header = ['m/z', 'Num. of Hits', 'Hits'];
    const csvRows = data.map(row => {
      return `${row.mz},"${row.numOfHits}","${row.hits.join(', ')}"`;
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
    const selectedRows = selectedRowKeys.map(key => dataSource[key]);
  
    const onSelectChange = (newSelectedRowKeys) => {
      setSelectedRowKeys(newSelectedRowKeys);
    };
  
    const rowSelection = {
      selectedRowKeys,
      onChange: onSelectChange,
    };
  
    return (
      <div style={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
        <Table rowSelection={rowSelection} columns={columns} dataSource={dataSource} />
        <div style={{ marginTop: '20px', display: 'flex', justifyContent: 'center' }}>
          <Button 
            type="primary" 
            icon={<DownloadOutlined />} 
            onClick={() => downloadCSV(selectedRows)}
          >
            Download Selected
          </Button>
        </div>
      </div>
    );
  };
  
  

// Parse HTML to an array of objects compatible with antd's Table
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
      hits: cells[2].textContent.split('\n') // split by \n to get a list of hits
    };
  });
}

export default Result;
