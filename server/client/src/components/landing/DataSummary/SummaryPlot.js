import React, { Component } from 'react';
import Plot from 'react-plotly.js';

class SummaryPlot extends Component {
  render() {

    const trace1 = {
        x: ['PanCancer','SKCM', 'Liver', 'BRCA', 'BLCA'],
        y: [3, 0, 1, 1, 0],
        width: [0.2,0.2,0.2,0.2,0.2],
        name: 'promoter',
        type: 'bar'
    };
      
    const trace2 = {
        x: ['PanCancer','SKCM', 'Liver', 'BRCA', 'BLCA'],
        y: [3,1, 0, 0, 1],
        width: [0.2,0.2,0.2,0.2,0.2],
        name: '5UTR',
        type: 'bar',
        hovertemplate: '<i>Price</i>: $%{y:.2f}' +
                    '<br><b>X</b>: %{x}<br>' +
                    '<b>%{text}</b>',
        text: ["Text A", "Text B", "Text C", "Text D", "Text E"],
        hoverinfo: 'skip'
    };  



    const data1 = [trace1, trace2];
      
    const layout = {
        title: 'Predicted non-coding cancer drivers per cancer type',
        xaxis: {
            tickangle: -45
        },
        yaxis: {
            title: 'Number of genes'
        },
        barmode: 'stack',
        hovermode: 'closest',
        hoverlabel: { bgcolor: "#FFF" },
        bargap: 0.15,
        bargroupgap: 0.1
    };




    return (
        <Plot
          data={ data1 }
          layout={ layout }
        />
    );
  }
}

export default SummaryPlot;