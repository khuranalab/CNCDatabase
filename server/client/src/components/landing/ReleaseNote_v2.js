import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { TwitterTimelineEmbed } from 'react-twitter-embed';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';

class ReleaseNote extends Component {
  render() {
    return (
      
      <div className="container-fluid pt-3 pb-3 pl-5">
            <div className="row">
               <div className="col-lg-6 text-left">
                   
                   <h1><strong> News </strong></h1>
                   <hr/>
                   <span className="text-info">October 20, 2019 <br/> </span>
                   <ul className="timeline"> 
                     <li><span className="text-info">October 20, 2019</span>
                         <p>Update search interface</p>
                     </li>       
                   </ul>
                   <span className="text-info">March 10, 2019 <br/> </span>
                   <ul> 
                     <li>Update non-coding cancer driver data. (version 0.2) </li>       
                   </ul>
                   <span className="text-info">Feburary 01, 2019 <br/> </span>
                   <ul> 
                     <li>Add auto-suggest feature in the search. </li>       
                   </ul>
                   <span className="text-info">Jan 01, 2019 <br/> </span>
                   <ul> 
                     <li>Update non-coding cancer driver data. (version 0.1) </li>       
                   </ul>
                   

               </div>

            </div>
      </div>
      
    );
  }
}

export default ReleaseNote;
