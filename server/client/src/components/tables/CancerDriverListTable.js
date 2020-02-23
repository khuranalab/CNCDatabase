import React from 'react';
//import 'react-tabulator/lib/styles.css';
//import 'react-tabulator/lib/css/bootstrap/tabulator_bootstrap4.min.css';
//import 'react-tabulator/lib/css/bulma/tabulator_bulma.min.css';
//import 'react-tabulator/lib/css/tabulator_modern.min.css';
//import 'react-tabulator/lib/css/tabulator.min.css';
import { ReactTabulator, reactFormatter } from 'react-tabulator';
//import '../css/tabulator.css';
import '../css/tabulator_setup.css';

//import { Link } from 'react-router-dom';
import { connect } from 'react-redux';
import { fetchCancerDriverList, fetchBygeneSymbol } from '../../actions';
import { CSVLink } from 'react-csv';
import ReactTooltip from 'react-tooltip';


class CancerDriverListTable extends React.Component {
     
    componentDidMount() {
        // action
        //this.props.fetchCancerDriverList();
    }

    renderDataBygeneSmbol2(dataBygeneSymbol) {

        const data = dataBygeneSymbol

        console.log("table data", data);
        console.log("table data", data.length) 
        
        const totalDataLength = data.length

        const pmidFormatter = function(cell, formatterParams){
            //cell - the cell component
            //formatterParams - parameters set for the column
    
            const data = cell.getData();
            console.log("check",data);
            
            return "<a target='_blank' href='https://www.ncbi.nlm.nih.gov/pubmed/" + data.pmid + "'>" + cell.getValue() + "</a>";
            //return "<Link to='https://mapquest.com/latlng/" + data.col + "'>" + cell.getValue() + "</Link>";
        }

        const elementTypeFormatter = function(cell, formatterParams){
          //cell - the cell component
          //formatterParams - parameters set for the column
  
          const data = cell.getData();
          console.log("check",data);
          
          return "<button type=\"button\" class=\"btn btn-primary btn-sm pt-0\">" + cell.getValue() + "</button>";
          //return "<Link to='https://mapquest.com/latlng/" + data.col + "'>" + cell.getValue() + "</Link>";
         }

         const evidenceMethodFormatter = function(cell, formatterParams){
          //cell - the cell component
          //formatterParams - parameters set for the column
  
          const data = cell.getData();
          console.log("check evidence method",data);
          
          return "<p> test1 </p>";
              
          
          
          
          //"<button type=\"button\" class=\"btn btn-outline-primary btn-sm\">" + cell.getValue() + "</button>";
          //return "<Link to='https://mapquest.com/latlng/" + data.col + "'>" + cell.getValue() + "</Link>";
         }

        const TestFormatter = function(props: any) {
          console.log("react cell",props)
          const cellData = props.cell._cell.row.data;
          console.log("cellData",cellData)

          

          const name = cellData.evidencemethod
          const description = "description"

          return (
              <div>
                <a data-tip={`${name}|${description}`} data-for='evidencemethod'> {name} </a>
                
                
                
                
                <ReactTooltip 
                   id='evidencemethod' 
                   place="left" 
                   type="light" 
                   effect="solid"
                   border={true}
                   getContent={(dataTip) => (
                          <div>
                            This is {dataTip}
                            <p>You can put every thing here</p>
                            <ul>
                              <li>Word</li>
                              <li>Chart</li>
                              <li>Else</li>
                            </ul>
                          </div>
                   )}
                />
                   
              </div>
          );          
        }

        const HeaderFormatter1 = function(props: any) {
          console.log("react header cell",props)
          //const cellData = props.cell._cell.row.data;
          //console.log("header cellData",cellData)

          

          //const name = cellData.evidencemethod
          const name = "test"
          const description = "description"

          return (
              <div>
                
                   
              </div>
          );          
        }

        const tooltipsFormatter = function(column, formatterParams){
            const data = column.getDefinition();
            console.log("tooltips check",data)

            return (
              <div>
                test header
              </div> 
            )
        } 


    //pmid, cancertype, gene, element, mutatedsamplesize, cohortsize, evidencetype, evidencemethod

        const columns = [
            {
              field: 'pmid',
              title: 'PMID',
              formatter: pmidFormatter,
              //headerFilter: true,
              align: 'center',
              width: 140,
              headerTooltip:true
            }, 
            {
              field: 'cancertype',
              title: 'Cancer Type',
              //headerFilter: true,
              width: 160,
              align: 'center'
              
            }, 
            {
              field: 'gene',
              title: 'Gene Name',
              //headerFilter: true,
              width: 140,
              align: 'center'
              
            }, 
            {
              field: 'element',
              title: 'Element',
              //headerFilter: true,
              //formatter: elementTypeFormatter,
              width: 180,
              align: 'center'
              
            }, 
            {
              field: 'mutatedsamplesize',
              title: 'Mutated Samples',
              //headerFilter: true,
              align: 'center',
              width: 180,
              headerTooltip: tooltipsFormatter
            }, 
            {
              field: 'cohortsize',
              title: 'Cohort Samples',
              //headerFilter: true,
              width: 180,
              align: 'center'
              
            },
            {
              field: 'evidencetype',
              title: 'Evidence Type',
              //headerFilter: true,
              width: 220,
              align: 'center'
              
            }, 
            {
              field: 'evidencemethod',
              title: 'Evidence Method',
              //headerFilter: true,
             //width: 250,
              tooltips: true,
              align: 'center',

              formatter: reactFormatter(

                <TestFormatter/>
                
              )
            }
          ];

        const options = {
            pagination: "local",
            paginationSize: 20,
            paginationButtonCount:3,
            paginationSizeSelector:[5,10,20,30,50,100],
            selectable:false,
            //footerElement: "<button>Custom Button</button>",
          };

        return (
            <div className="container-fluid">

              

              <div className="d-flex justify-content-between"> 
                  <div className="p-2">
                    Returned {totalDataLength} entries
                  </div>
                  
                  <div className="p-2">  
                    <CSVLink
                      data={data} 
                      filename={"SearchResults.csv"}
                      className="btn btn-outline-primary btn-sm"
                      target="_blank"
                    >
                      Download CSV 

                    </CSVLink>
                  </div>  
               
              </div>

              <div className="row">
                <div className="col-lg-12">    
                    <ReactTabulator
                            data={data}
                            initialSort={[{column:"element",dir:"asc"}]}
                            columns={columns}
                            tooltips={true}
                            headerTooltip={true}
                            layout= {"fitColumns"}
                            autoColumns={true}
                            options={options}
                    />
                </div>


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
  
export default connect(mapStateToProps, { fetchCancerDriverList })(CancerDriverListTable);

//export default ReactTabulatorTable;