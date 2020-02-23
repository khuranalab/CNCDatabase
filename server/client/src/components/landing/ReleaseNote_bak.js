import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { TwitterTimelineEmbed } from 'react-twitter-embed';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';

class ReleaseNote extends Component {
  render() {
    return (
      <div className="py-3 bg-light">
         <div className="container col-lg-10">
            <div className="row">
               <div className="col-lg-6 text-left">
                   
                   <h5> Release notes </h5>
                   <hr/>
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

               <div className="col-lg-6 text-left">
                   <h5> News @cncdatabase</h5>
                   <hr/>

                   <TwitterTimelineEmbed
                      sourceType="profile"
                      screenName="ekta_khurana"
                      noHeader
                      options={{height: 400}}
                   />
                   
                   
               </div>

            </div>
         </div>
      </div>
    );
  }
}

export default ReleaseNote;
