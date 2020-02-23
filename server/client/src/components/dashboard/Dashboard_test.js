import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import SearchBar from '../searches/Search_bar';
import QueryBarByGene from '../searches/Query_bar_by_gene';
//import TestTable from './tables/TestTable';
//import TestTable2 from './tables/TestTable2';
//import AndTableContainer3 from './tables/AntdTableContainer3';

//import InfoTable from './InfoTable';
//import ReactBootstrapTable from './tables/ReactBootstrapTable';
import CancerDriverByGeneTable from '../tables/CancerDriverByGeneTable';
import GeneSummary from './GeneSummary';
//import TabContent from './TabContent';
import { Tabs, Tab } from 'react-bootstrap-tabs';
import TestPlot from '../studies/TestPlot';


class Dashboard_test extends Component {
  render() {
    return (
      <div className="container col-10">
        <br/> 
        <h2>Search cancer driver</h2>
        <br/>
        {/* <!-- Nav tabs --> */}

        <Tabs onSelect={(index, label) => console.log(label + ' selected')}>
           <Tab label="By Gene">
                    <br />
                    <div className="card">
                    <div className="card-body">
                        <QueryBarByGene />
                    </div>
                    </div>
                    <br /> 

                    <div className="card mt-3">
                    <h5 className="card-header">Gene overview (provided by HGNC)</h5>
                    <div className="card-body">
                        <GeneSummary />
                    </div>
                    </div>
                    <br />

                    <div className="card mt-3">
                    <h5 className="card-header">Graphical Summary</h5>
                    <div className="card-body">
                        <TestPlot />
                    </div>
                    </div>
                    <br />


                    <div className="card mt-3">
                    <h5 className="card-header">Search Results</h5>
                    <div className="card-body">
                        <CancerDriverByGeneTable />
                    </div>
                    </div>
                    <br />
           </Tab>

           <Tab label="By Cancer Type">
                    <br />
                    <div className="card">
                    <div className="card-body">
                        <SearchBarA />
                    </div>
                    </div>
                    <br /> 
           </Tab>
           
           <Tab label="By Element Type">
                    Tab 2 content
           </Tab>
        
        
        </Tabs>
        
        
      </div>
    );
  }
}

//export default Dashboard_test;
