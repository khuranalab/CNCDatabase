import React from 'react';


//import { Link } from 'react-router-dom';
import { connect } from 'react-redux';
import { fetchCancerDriverList, fetchBygeneSymbol } from '../../actions';
//import { CSVLink } from 'react-csv';
import * as d3 from 'd3';
import Plot from 'react-plotly.js';

class CancerDriverListPlot extends React.Component {
     
    componentDidMount() {
        // action
        //this.props.dispatch(fetchCancerDriverList());
    }

    renderDataBygeneSmbol2(dataBygeneSymbol) {

        const data = dataBygeneSymbol

        console.log("plot data", data);
        console.log("plot data length", data.length) 
        console.log("plot data element", data.map(d => d.element)) 
        const totalDataLength = data.length

        const numOfpublications = d3.map(data, function(d) {return d.pmid; }).keys().length
        console.log("unique papers", numOfpublications);

         //pmid, cancertype, gene, element, mutatedsamplesize, cohortsize, evidencetype, evidencemethod

        
        let dataArrayCancerType = data.map(d => d.cancertype) || 1;
        let dataArrayElement = data.map(d => d.element) || 1;
        let dataArrayEvidenceType = data.map(d => d.evidencetype) || 1;
        let dataArrayEvidenceMethod = data.map(d => d.evidencemethod) || 1;

        function count(arr) {
            return arr.reduce((prev, curr) => (prev[curr] = ++prev[curr] || 1, prev), {})
          }
        
        let countsCancerType = count(dataArrayCancerType); 
        let countsElement = count(dataArrayElement); 
        let countsEvidenceType = count(dataArrayEvidenceType); 
        let countsEvidenceMethod = count(dataArrayEvidenceMethod); 


        let countsCancerTypeKeys = Object.keys(countsCancerType);  
        let countsCancerTypeValues = Object.values(countsCancerType);


        let countsElementKeys = Object.keys(countsElement);  
        let countsElementValues = Object.values(countsElement);

        let countsEvidenceTypeKeys = Object.keys(countsEvidenceType);  
        let countsEvidenceTypeValues = Object.values(countsEvidenceType);

        let countsEvidenceMethodKeys = Object.keys(countsEvidenceMethod);  
        let countsEvidenceMethodValues = Object.values(countsEvidenceMethod);

        //console.log("element counts", count(arr));        
        
        //console.log("element counts key", c1);
        //console.log("element counts value", c2); 
        
            //const elementTypePie = [16, 15, 12, 6, 5, 4, 50];
            //const elementTypePieLabel = ['US', 'China', 'European Union', 'Russian Federation', 'Brazil', 'India', 'Rest of World' ];

            //const elementTypePie = countsElementValues;
            //const elementTypePieLabel = countsElementKeys;

            const dataPieCancerType = [{
              type: "pie",
              hole: .4,
              values: countsCancerTypeValues,
              labels: countsCancerTypeKeys,
              domain: {column: 0},
              name: 'Cancer Type',
              textinfo: "label+value",
              textposition: "outside",
              hoverinfo: 'label+percent',
              automargin: true
            }];

            const layoutPieCancerType = {
              title: 'Cancer Type',
              height: 300,
              width: 300,
              showlegend: false
            };


            const dataPieElement = [{
                type: "pie",
                hole: .4,
                values: countsElementValues,
                labels: countsElementKeys,
                domain: {column: 0},
                name: 'Element Type',
                textinfo: "label+value",
                textposition: "outside",
                hoverinfo: 'label+percent',
                automargin: true
              }];

            const layoutPieElement = {
                title: 'Element Type',
                height: 300,
                width: 400,
                showlegend: false
              };

              const dataPieEvidenceType = [{
                type: "pie",
                hole: .4,
                values: countsEvidenceTypeValues,
                labels: countsEvidenceTypeKeys,
                domain: {column: 0},
                name: 'Evidence Type',
                textinfo: "label+value",
                textposition: "outside",
                hoverinfo: 'label+percent',
                margin: {
                  l: 50,
                  r: 50,
                  b: 50,
                  t: 50,
                  pad: 4
                },
                automargin: false
              }];

              const layoutPieEvidenceType = {
                title: 'Evidence Type',
                height: 300,
                width: 400,
                cliponaxis: false,
                showlegend: false
              };

              const dataPieEvidenceMethod = [{
                type: "pie",
                hole: .4,
                values: countsEvidenceMethodValues,
                labels: countsEvidenceMethodKeys,
                domain: {column: 0},
                name: 'Evidence Method',
                textinfo: "label+value",
                textposition: "outside",
                hoverinfo: 'label+percent',
                automargin: true
              }];

            const layoutPieEvidenceMethod = {
                title: 'Evidence Method',
                height: 300,
                width: 500,
                showlegend: false
              };


        return (
            <div className="container-fluid" style={{ marginTop: 20 }}>

              <div className="d-flex justify-content-around flex-row">
               <Plot
                   data={dataPieCancerType}
                   layout={layoutPieCancerType} 
               />
             

              
               <Plot
                   data={dataPieElement}
                   layout={layoutPieElement} 
               />

               <Plot
                   data={dataPieEvidenceType}
                   layout={layoutPieEvidenceType} 
               />

               <Plot
                   data={dataPieEvidenceMethod}
                   layout={layoutPieEvidenceMethod} 
               />

              </div>

              



            </div>

        );

    }    

    render (){
      
        return(
            <div>{ this.renderDataBygeneSmbol2(this.props.dataBygeneSymbol) }</div>
        );


    }

}

//export default ReactBootstrapTable;


const mapStateToProps = ( state ) => {
    return { dataBygeneSymbol: state.dataCancerDriverList };
}
  
export default connect(mapStateToProps, { fetchCancerDriverList })(CancerDriverListPlot);

//export default ReactTabulatorTable;