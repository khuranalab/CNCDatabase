import React from 'react';
//import 'react-tabulator/lib/styles.css';
import { ReactTabulator} from 'react-tabulator';
//import { Link } from 'react-router-dom';
import { connect } from 'react-redux';
import { fetchBygeneSymbol } from '../../actions';


class ReactTabulatorTable extends React.Component {
     
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




        const columns = [
            {
              field: 'pmid',
              title: 'PMID',
              formatter: pmidFormatter
            }, 
            {
              field: 'cancertype',
              title: 'Cancer Type'
              
            }, 
            {
              field: 'cancertype2',
              title: 'Type'
              
            },
            {
              field: 'gene',
              title: 'Gene Name'
              
            }, 
            {
              field: 'ensemblID',
              title: 'Ensembl ID'
              
            }, 
            {
              field: 'mutationType',
              title: 'Mutation Type'
              
            },
            {
              field: 'element',
              title: 'Element'
              
            }, 
            {
              field: 'mutatedsamplesize',
              title: 'Mutated Sample Size'
              
            }, 
            {
              field: 'cohortsize',
              title: 'Cohort Size'
              
            },
            {
              field: 'mutationrate',
              title: 'Mutation Rate'
              
            }, 
            {
              field: 'evidencetype',
              title: 'Evidence Type'
              
            }, 
            {
              field: 'evidencemethod',
              title: 'Evidence Method'
              
            }
          ];

        const options = {
            pagination: "local",
            paginationSize: 50,
            paginationButtonCount:3,
            footerElement: "<button>Custom Button</button>",
          };

        return (
            <div className="container-fluid" style={{ marginTop: 20 }}>
               <ReactTabulator
                      data={data}
                      columns={columns}
                      tooltips={true}
                      layout={"fitData"}
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
  
export default connect(mapStateToProps, { fetchBygeneSymbol })(ReactTabulatorTable);

//export default ReactTabulatorTable;