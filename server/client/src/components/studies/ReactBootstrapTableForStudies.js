import React from 'react';
import BootstrapTable from 'react-bootstrap-table-next';
import paginationFactory from 'react-bootstrap-table2-paginator';
import filterFactory, { textFilter } from 'react-bootstrap-table2-filter';
import { connect } from 'react-redux';
import { fetchAllStudies } from '../../actions';


class ReactBootstrapTableForStudies extends React.Component {
     
    componentDidMount() {
        // action
        this.props.fetchAllStudies();
    }

    renderDataByStudies(dataByStudies) {

        const data = dataByStudies

        const columns = [
          {
            dataField: 'pmid',
            text: 'PMID',
            sort: true,
            filter: textFilter(),
            headerStyle: {
              width: '100px' 
            },
            
          }, 
          {
            dataField: 'year',
            text: 'Year',
            sort: true,
            filter: textFilter(),
            headerStyle: {
                width: '80px' 
            }
          }, 
          {
            dataField: 'first_author_last_name',
            text: 'Last Name',
            sort: true,
            filter: textFilter(),
            headerStyle: {
              width: '100px' 
          }
          },
          {
            dataField: 'journal_name',
            text: 'Journal',
            sort: true,
            filter: textFilter(),
            headerStyle: {
              width: '150px' 
          }
          }, 
          {
            dataField: 'title',
            text: 'Title',
            sort: true,
            filter: textFilter()
          }, 
          {
            dataField: 'link',
            text: 'URL',
            sort: true,
            filter: textFilter()
          },
          
        ];
          
        const defaultSorted = [{
            dataField: 'year',
            order: 'desc'
          }
        ];
    
          const customTotal = (from, to, size) => (
            <span className="react-bootstrap-table-pagination-total">
              Showing { from } to { to } of { size } Results
            </span>
          );
          
          const options = {
            paginationSize: 20,
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
      
      
        return(
            <div>{ this.renderDataByStudies(this.props.dataByStudies) }</div>
        );


    }

}

//export default ReactBootstrapTable;

const mapStateToProps = ( state ) => {
    return { dataByStudies: state.dataByStudies };
}
  
export default connect(mapStateToProps, { fetchAllStudies })(ReactBootstrapTableForStudies);