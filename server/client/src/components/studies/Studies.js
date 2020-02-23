import React, { Component } from 'react';
import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import SearchBar from '../searches/Search_bar';
//import TestTable from './tables/TestTable';
//import TestTable2 from './tables/TestTable2';
//import AndTableContainer3 from './tables/AntdTableContainer3';

//import InfoTable from './InfoTable';
import ReactBootstrapTableForStudies from './ReactBootstrapTableForStudies';
import TestPlot from './TestPlot';

import { fetchAllStudies } from '../../actions';

class Studies extends Component {

  componentDidMount() {
    // action
    this.props.fetchAllStudies();
  }



  render() {
    return (
      <div className="container col-10">
        <div className="pt-4">
          <h1>Overview</h1>
          <hr/>
          Cornell Non-coding database
        </div>
        <br />
        <div>
          <TestPlot/>
        </div>  
        <br/>
        
        <br /> 
        <div className="mt-3">
          <h5 >References</h5>
          <hr/>
          <div>
            
          </div>  
          <div >
          <ReactBootstrapTableForStudies />
          </div>
        </div>
        <br />
        
       
        
        <br />
      </div>
    );
  }
}

//export default Studies;

const mapStateToProps = ( state ) => {
  return { dataByStudies: state.dataByStudies };
}

export default connect(mapStateToProps, { fetchAllStudies })(Studies);