import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
import { UndrawBrowserStats } from 'react-undraw-illustrations';

class ShowcaseCallToAction extends Component {
  render() {
    return (
      <div className="py-5 bg-light">
         <div className="container col-lg-10">
            <div className="row">
               
               <div className="col-lg-6">
                  <UndrawBrowserStats
                    primaryColor='#4285f4'
                    height= '150px'
                  />
               </div>  


               <div className="col-lg-6 text-left">
                   <h5> Cornell Non-Coding cancer driver Database </h5>
                   <hr/>
                   <p>  Try query by gene name or by cancer type  .   
                        
                   </p>
                   <button type="submit" className="btn btn-primary btn-lg">
                        Get Started
                   </button>

               </div>
            </div>
         </div>
      </div>
    );
  }
}

export default ShowcaseCallToAction;
