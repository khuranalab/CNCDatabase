import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
//import { UndrawBrowserStats } from 'react-undraw-illustrations';
import SummaryPlot from './DataSummary/SummaryPlot';
//import SummaryPiePlot from './DataSummary/SummaryPiePlot';

class DataSummary extends Component {
  render() {
    return (
        <div className="container-fluid pt-3 pb-3 pl-5">
          <div className="row">
            <div className="col-lg-10">
                <h3><strong>Data summary</strong></h3>
                <br/>
                <div className="container d-flex justify-content-center px-5">
                  
                    

                    <SummaryPlot/>

                    


                    
                </div>
            </div>
          </div>
        </div>  
    );
  }
}

export default DataSummary;
