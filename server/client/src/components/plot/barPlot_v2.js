import React from 'react';
import Plot from 'plotly.js-dist';

class ExamplePlot extends React.Component {
  render() {

    const trace1 = {
        x: ['PanCancer','SKCM', 'Liver', 'BRCA', 'BLCA'],
        y: [3, 0, 1, 1, 0],
        name: 'promoter',
        type: 'bar'
      };
      
      const trace2 = {
        x: ['PanCancer','SKCM', 'Liver', 'BRCA', 'BLCA'],
        y: [3,1, 0, 0, 1],
        name: '5UTR',
        type: 'bar'
      };  



    const data1 = [trace1, trace2];
      
    const layout = {
        title: 'Evidence support per cancer type',
        xaxis: {
            tickangle: -45
        },
        yaxis: {
            title: 'Number of unique methods'
        },
        barmode: 'group',
        bargap: 0.15,
        bargroupgap: 0.1
      };
      
    //Plotly.newPlot('myDiv', data, layout);  






    return (
      <Plot
        data={data1}
        layout={ layout }
      />
    );
  }
}

export default ExamplePlot;