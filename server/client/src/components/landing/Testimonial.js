import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';

class Testimonial extends Component {
  render() {
    return (
      <div className="py-5">
         <div className="container col-lg-10">
            <div className="row">
               <div className="col-lg-6 text-left">
                    <blockquote class="blockquote text-center">
                    <p class="mb-0">Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer posuere erat a ante.</p>
                    <footer class="blockquote-footer">Someone famous in <cite title="Source Title">Source Title</cite></footer>
                    </blockquote>
                   
                   

               </div>
            </div>
         </div>
      </div>
    );
  }
}

export default Testimonial;
