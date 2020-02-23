import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
import { UndrawData } from 'react-undraw-illustrations';


class Showcase extends Component {
  render() {
    return (
      <div className="py-5 bg-light">
         <div className="container col-lg-10">
            <div className="row">
               <div className="col-lg-8 text-left">
                   <h1> CNCDatabase: </h1>
                   <h2> Cornell Non-Coding cancer driver Database </h2>
                   <hr/>
                   
                   <button type="submit" className="px-2 btn btn-primary btn-lg">
                        Get Started
                   </button>
               
               

               </div>

               <div className="col-lg-4">
                    <UndrawData
                        primaryColor='#6c68fb'
                        height= '250px'
                    />
               </div>   

            </div>
         </div>
      </div>
    );
  }
}

export default Showcase;
