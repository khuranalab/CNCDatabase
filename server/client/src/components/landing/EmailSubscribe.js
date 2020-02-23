import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';

class EmailSubscribe extends Component {
  






  render() {
    return (
      <div className="py-3">
         <div className="container col-lg-10">
            <div className="row">
               <div className="col-lg-3">

               </div>

               <div className="col-lg-3 text-right align-middle">
                   
                   <h5> Get Updates on Email </h5>
                   
                   
                   

               </div>

               


               <div className="col-lg-3 align-middle">
                   
                  <form onSubmit={null} className="input-group">
                      <input
                        placeholder="Enter your email"
                        className="form-control"
                        value={null}
                        onChange={null}
                      />
                      <span className="input-group-btn">
                        <button type="submit" className="btn btn-primary">
                          Subscribe
                        </button>
                      </span>
                  </form>
                   
               </div>

               <div className="col-lg-3 align-middle">
               </div>



            </div>
         </div>
      </div>
    );
  }
}

export default EmailSubscribe;
