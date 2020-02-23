import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
import SearchBar from '../searches/Search_bar';
//import TestTable from './tables/TestTable';
//import TestTable2 from './tables/TestTable2';
//import AndTableContainer3 from './tables/AntdTableContainer3';

//import InfoTable from './InfoTable';
//import ReactBootstrapTable from './tables/ReactBootstrapTable';
import CancerDriverByGeneTable from '../tables/CancerDriverByGeneTable';
import GeneSummary from './GeneSummary';

class Dashboard extends Component {
  render() {
    return (
      <div className="container col-10">
        <div style={{ textAlign: 'center' }}>
          <h1>CNCDatabase</h1>
          Cornell Non-coding database
        </div>
        <br />
        <div className="card">
          <div className="card-body">
            <SearchBar />
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
          <h5 className="card-header">Search Results</h5>
          <div className="card-body">
            <CancerDriverByGeneTable />
          </div>
        </div>
        <br />
        
       
        
        <br />
      </div>
    );
  }
}

export default Dashboard;
