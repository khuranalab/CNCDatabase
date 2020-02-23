import React, { Component } from 'react';
import { connect } from 'react-redux';
import { Link } from 'react-router-dom';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
import Showcase from './Showcase';


class Landing extends Component {
  render() {
    return (
      <div>


        <Showcase/>



      </div>




      <div className="container col-md-10">

       <div>  
        <div className="p-5">
        </div>  

        <div className="text-center py-3"> 
        <h1>CNCDatabase</h1>
        Cornell Non-coding database
        </div>

        <div className="container col-md-8">
          <InputGroup>
          <Input placeholder="Search by gene name" />
            <Button color="primary">
              <FontAwesomeIcon icon="search"/>
            </Button>{' '}
          </InputGroup>

          <div className="pt-2 pb-5">
          Example: TERT or WDR74
          </div>
        </div>
      
        
      </div>

      <div className="container pt-5">
            <h2 className="text-center">Features</h2>
            <div className="d-flex justify-content-between pt-3">
              <div className="col-md-4">
                <h4>Feature One</h4>
                <p>feature explnation 1 </p>
              </div>
              <div className="col-md-4">
                <h4>Feature Two</h4>
                <p>feature explnation 1 </p>
              </div>
              <div className="col-md-4">
                <h4>Feature Three</h4>
                <p>feature explnation 1 </p>
              </div>
            </div>

      </div>  

        <div className="p-5">
        </div>  

      
      </div>
    );
  }
}

//export default Landing;
