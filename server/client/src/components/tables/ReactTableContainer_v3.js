import React, { useMemo } from "react";
import { connect } from "react-redux";
import { fetchCancerDriverList } from "../../actions";

import Table from "./Table_v5";
import "./css/Table.css";
import ReactTooltip from 'react-tooltip';

const ReactTableContainer = ({dataCancerDriverList}) => {

  const data = dataCancerDriverList;

  const columns = useMemo(
    () => [
      
          {
            Header: <span title="PMID">PMID</span>,
            accessor: "pmid",
            Cell: ({ cell: { value } }) => {
                 return (
                   <div> 
                     <a target='_blank' href={"https://www.ncbi.nlm.nih.gov/pubmed/" + value }>{value}</a>
                  </div>
                 ); 
            }
          },
          {
            Header: "Cancer Type",
            accessor: "cancertype",
            Cell: ({ cell: { value } }) => {
              return (
                <div>
                    <a data-tip={`${value}`} data-for='cancertype'> {value} </a>
                    
                    <ReactTooltip 
                    id='cancertype' 
                    place="right" 
                    type="light" 
                    effect="solid"
                    border={true}
                    getContent={(value) => (
                            <div>
                                <p><strong>Cancer Type</strong></p>
                                <hr/>
                                <ul>
                                  <li>Short name: {value}</li>
                                </ul>
                            </div>
                    )}
                    />
                   
                </div>
              ); 
         }
          },
        
          {
            Header: "Gene Name",
            accessor: "gene"
          },
          {
            Header: "COSMIC(Cancer Gene)",
            accessor: "cosmic"
          },
          {
            Header: "Element",
            accessor: "element",
            //Cell: ({ cell: { value } }) => <Genres values={value} />
          },
          {
            Header: "Cohort Size",
            accessor: "cohortsize"
          },
          {
            Header: "Mutated Samples",
            accessor: "mutatedsamplesize",

          },
          {
            Header: "Evidence Type",
            accessor: "evidencetype"
          },
          {
            Header: "Evidence Method",
            accessor: "evidencemethod"
          }
        
      
    ],
    []
  );

  //const [data, setData] = useState([]);

  //useEffect(() => {
  //  (async () => {
  //    const result = await axios("https://api.tvmaze.com/search/shows?q=snow");
  //    setData(result.data);
  //  })();
  //}, []);

  console.log("react table data:",data)

  return (
    <div className="react-table-main">
      <Table columns={columns} data={data} />
    </div>
  );
}

//export default ReactTableContainer;

const mapStateToProps = ( state ) => {
    return { dataCancerDriverList: state.dataCancerDriverList };
}
export default connect( mapStateToProps, { fetchCancerDriverList })(ReactTableContainer);