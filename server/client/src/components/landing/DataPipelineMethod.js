import React, { Component } from 'react';
//import Collapsible from 'react-collapsible';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
//import { UndrawBrowserStats } from 'react-undraw-illustrations';
//import Panel from './UseCases/CollapsibleComponent/Panel';
import DataPipelineFigure from './img/dataPipelineFigure.png'

class DataPipelineMethod extends Component {
  render() {
    return (
      <div className="container-fluid pl-5">
         <div className="row">
              <div className="col-lg-10">
                   
                   <h3><strong> Data procecssing pipeline </strong></h3>

                   
                   <div className="container d-flex px-1">
               
                

                     <img src={DataPipelineFigure} alt="datapipelinefigure" width="85%" height="85%"/>


                   </div>

              </div>      

         </div>
         
      </div>
    );
  }
}

export default DataPipelineMethod;
