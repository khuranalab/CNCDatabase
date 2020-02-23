import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import SearchBar from '../searches/Search_bar';
//import QueryBarByGene from '../searches/Query_bar_by_gene';
//import QueryForm2 from '../searches/QueryForm_v3';
//import QueryForm from '../searches/QueryForm';

//import QueryForm from '../searches/QueryFormV7';

import QueryForm from '../searches/QueryFormikForm';
//import MyForm from '../searches/QueryHookFormTest_v2';

//import QueryHookForm from '../searches/QueryHookFormTest';

//import TestTable from './tables/TestTable';
//import TestTable2 from './tables/TestTable2';
//import CancerDriverListTable from '../tables/AntdTableContainer3';

//import InfoTable from './InfoTable';
//import CancerDriverListTable from '../tables/ReactBootstrapTable';

//import CancerDriverByGeneTable from '../tables/CancerDriverByGeneTable';

//import CancerDriverListTable from '../tables/CancerDriverListTable';

import ReactTableContainer from '../tables/ReactTableContainer_v3';

//import CancerDriverListTable from '../tables/CancerDriverListTableMUI';

//import GeneSummary from './GeneSummary';
//import TabContent from './TabContent';
//import { Tabs, Tab } from 'react-bootstrap-tabs';
//import TestPlot from '../studies/TestPlot';
//import ExamplePlot from '../plot/barPlot';
//import CancerDriverListPlot from '../plot/CancerDriverListPlot';

import CancerDriverListPlot from '../plot/TestPlot_v12';

//import QueryBarByGene from '../searches/Query_bar_by_gene';

class Dashboard extends Component {
  render() {
    return (
      <div className="container-fluid pl-5 py-5">
        <div className="row">
          <div className="col-lg-12">
             <h1><strong>Search cancer driver</strong></h1>
             <hr className="style1"/>
       
             <div className="card">
                    <h5 className="card-header"><strong>Search cancer driver list</strong></h5>
                    <div className="card-body">
                       
                        
                          <QueryForm />
                        
                    </div>
             </div>
           
             <div className="card mt-3">
                     <h5 className="card-header"><strong>Summary</strong> </h5>
                     <div className="card-body">
                        
                          <CancerDriverListPlot />
                        
                     </div>
             </div>
  
             <div className="card mt-3">
                     <h5 className="card-header"><strong>Results</strong> </h5>
                      <div className="card-body">
                        
                          <ReactTableContainer />
                       
                      </div>
             </div>


           </div>
        </div>
      </div>
    );
  }
}

export default Dashboard;
