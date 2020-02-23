import React from 'react';
//import 'react-tabulator/lib/styles.css';
//import 'react-tabulator/lib/css/bootstrap/tabulator_bootstrap4.min.css';
//import 'react-tabulator/lib/css/tabulator.min.css';
import { ReactTabulator} from 'react-tabulator';
//import { Link } from 'react-router-dom';
import { connect } from 'react-redux';
import { fetchBygeneSymbol } from '../../actions';


class CancerDriverByGeneTable extends React.Component {
     
    componentDidMount() {
        // action
        this.props.fetchBygeneSymbol();
    }

    renderDataBygeneSmbol2(dataBygeneSymbol) {

        const data = dataBygeneSymbol

        console.log("table data", data); 

        const pmidFormatter = function(cell, formatterParams){
            //cell - the cell component
            //formatterParams - parameters set for the column
    
            const data = cell.getData();
            console.log("check",data);
            
            return "<a target='_blank' href='https://www.ncbi.nlm.nih.gov/pubmed/" + data.pmid + "'>" + cell.getValue() + "</a>";
            //return "<Link to='https://mapquest.com/latlng/" + data.col + "'>" + cell.getValue() + "</Link>";
    }


    //pmid, cancertype, gene, element, mutatedsamplesize, cohortsize, evidencetype, evidencemethod

        const columns = [
            {
              field: 'pmid',
              title: 'PMID',
              formatter: pmidFormatter,
              headerFilter: true,
              align: 'center',
              headerTooltip:true
            }, 
            {
              field: 'cancertype',
              title: 'Cancer Type',
              headerFilter: true,
              align: 'center'
              
            }, 
            {
              field: 'gene',
              title: 'Gene Name',
              headerFilter: true,
              align: 'center'
              
            }, 
            {
              field: 'element',
              title: 'Element',
              headerFilter: true,
              align: 'center'
              
            }, 
            {
              field: 'mutatedsamplesize',
              title: 'Mutated Sample Size',
              headerFilter: true,
              align: 'center'
              
            }, 
            {
              field: 'cohortsize',
              title: 'Cohort Size',
              headerFilter: true,
              align: 'center'
              
            },
            {
              field: 'evidencetype',
              title: 'Evidence Type',
              headerFilter: true,
              align: 'center'
              
            }, 
            {
              field: 'evidencemethod',
              title: 'Evidence Method',
              headerFilter: true,
              align: 'center'
            }
          ];

        const options = {
            pagination: "local",
            paginationSize: 20,
            paginationButtonCount:3,
            paginationSizeSelector:[5,10,20,50,100],
            footerElement: "<button>Custom Button</button>",
          };

        return (
            <div className="container-fluid" style={{ marginTop: 20 }}>
               <ReactTabulator
                      data={data}
                      columns={columns}
                      tooltips={true}
                      headerTooltip={true}
                      layout={"fitColumns"}
                      autoColumns={true}
                      options={options}
                  />
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
    return { dataBygeneSymbol: state.dataBygeneSymbol };
}
  
export default connect(mapStateToProps, { fetchBygeneSymbol })(CancerDriverByGeneTable);

//export default ReactTabulatorTable;