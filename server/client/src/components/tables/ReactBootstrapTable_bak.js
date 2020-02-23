import React from 'react';
import BootstrapTable from 'react-bootstrap-table-next';
import paginationFactory from 'react-bootstrap-table2-paginator';
import filterFactory, { textFilter } from 'react-bootstrap-table2-filter';
import { connect } from 'react-redux';
import { fetchBygeneSymbo } from '../../actions';


class ReactBootstrapTable extends React.Component {
     
    componentDidMount() {
        // action
        //this.props.fetchBygeneSymbol();
    }

    renderDataBygeneSmbol2(dataBygeneSymbol) {

        const data = dataBygeneSymbol

        const columns = [
          {
            dataField: 'pmid',
            text: 'PMID',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'cancertype',
            text: 'Cancer Type',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'cancertype2',
            text: 'Type',
            sort: true,
            filter: textFilter()
          },
          {
            dataField: 'gene',
            text: 'Gene Name',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'ensemblID',
            text: 'Ensembl ID',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'mutationType',
            text: 'Mutation Type',
            sort: true,
            filter: textFilter()
          },
          {
            dataField: 'element',
            text: 'Element',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'mutatedsamplesize',
            text: 'Mutated Sample Size',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'cohortsize',
            text: 'Cohort Size',
            sort: true,
            filter: textFilter()
          },
          {
            dataField: 'mutationrate',
            text: 'Mutation Rate',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'evidencetype',
            text: 'Evidence Type',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'evidencemethod',
            text: 'Evidence Method',
            sort: true,
            filter: textFilter()
          }
        ];
          
        const defaultSorted = [{
            dataField: 'gene',
            order: 'desc'
          }
        ];
    
          const customTotal = (from, to, size) => (
            <span className="react-bootstrap-table-pagination-total">
              Showing { from } to { to } of { size } Results
            </span>
          );
          
          const options = {
            paginationSize: 5,
            pageStartIndex: 0,
            // alwaysShowAllBtns: true, // Always show next and previous button
            // withFirstAndLast: false, // Hide the going to First and Last page button
            // hideSizePerPage: true, // Hide the sizePerPage dropdown always
            // hidePageListOnlyOnePage: true, // Hide the pagination list when only one page
            firstPageText: 'First',
            prePageText: 'Back',
            nextPageText: 'Next',
            lastPageText: 'Last',
            nextPageTitle: 'First page',
            prePageTitle: 'Pre page',
            firstPageTitle: 'Next page',
            lastPageTitle: 'Last page',
            showTotal: true,
            paginationTotalRenderer: customTotal,
            sizePerPageList: [{
              text: '5', value: 5
            }, {
              text: '10', value: 10
            }, {
              text: 'All', value: data.length
            }] // A numeric array is also available. the purpose of above example is custom the text
          };

        return (
            <div className="container-fluid" style={{ marginTop: 20 }}>
               <BootstrapTable
               bootstrap4
               keyField="pmid"
               data={ data }
               columns={ columns }
               striped
               hover
               condensed
               filter={ filterFactory() }
               pagination ={ paginationFactory(options) }
               />
           </div>

        );

    }    

    render (){
      
      const products = [
        {   "id": "1",
            "name": "Book",
            "price": "18"
        },
        {
            "id": "2",
            "name": "Mobile",
            "price": "400"
        },
        {
            "id": "3",
            "name": "PC",
            "price": "1000"
        },
        {
            "id": "4",
            "name": "PS4",
            "price": "500"
        },
        {
            "id": "5",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "6",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "7",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "8",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "9",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "10",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "11",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "12",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "13",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "14",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "15",
            "name": "Chromebook",
            "price": "500"
        },
        {
            "id": "16",
            "name": "Chromebook",
            "price": "500"
        }

      ];

      


         

        return(
            <div>{ this.renderDataBygeneSmbol2(this.props.dataBygeneSymbol) }</div>
        );


    }

}

//export default ReactBootstrapTable;

const mapStateToProps = ( state ) => {
    return { dataBygeneSymbol: state.dataBygeneSymbol };
}
  
export default connect(mapStateToProps, { fetchBygeneSymbol })(ReactBootstrapTable);