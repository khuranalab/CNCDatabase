import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';

class DataSummary extends Component {
  render() {
    return (
      <div className="pt-3 pb-3">
         <div className="container col-lg-10">
            <h4 className="text-left">Database summary</h4>
            <br/>
            <div className="container d-flex justify-content-center px-5">
               <div className="col-lg-3 box-auto text-center">
                   <h1 className="text-dark"> 928 </h1>
                   <hr/>
                   <p>  non-coding drivers from computational predictions   
                        
                   </p>
                   
               </div>

               <div className="col-lg-3 box-auto text-center">
                   <h1> 9 </h1>
                   <hr/>
                   <p>  non-coding drivers with gene expression supports  
                        
                   </p>
                   
               </div>

               <div className="col-lg-3 box-auto text-center">
                   <h1> 18 </h1>
                   <hr/>
                   <p>  non-coding drivers with experimental validations   
                        
                   </p>
                   
               </div>

               <div className="col-lg-3 box-auto text-center">
                   <h1> 27 </h1>
                   <hr/>
                   <p>  cancer types   
                        
                   </p>
                   
               </div>

               <div className="col-lg-3 box-auto text-center">
                   <h1> 18 </h1>
                   <hr/>
                   <p>  publications   
                        
                   </p>
                   
               </div>





            </div>
         </div>
      </div>
    );
  }
}

export default DataSummary;
