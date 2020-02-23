import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
//import { UndrawChoose, UndrawChecklist, UndrawWorkChat } from 'react-undraw-illustrations';


class Citation extends Component {
  render() {
    return (
        <div className="container-fluid pt-3 pb-3 bg-light pl-5">
          <div className="row">
            
                <div className="col-lg-2 py-4">
                    <h1 className="text-left"><strong> Citation </strong></h1>
                </div>    
                <div className="col-lg-10">    
                      <span>   Liu et al, Database (2019) 
                          
                      </span>
                    
                   
                 </div>  

                

          </div>
        </div>
    );
  }
}

export default Citation;
