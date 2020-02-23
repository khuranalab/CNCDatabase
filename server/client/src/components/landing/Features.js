import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
//import { UndrawChoose, UndrawChecklist, UndrawWorkChat } from 'react-undraw-illustrations';


class Features extends Component {
  render() {
    return (
      
        <div className="container-fluid bg-light pt-3 pb-3 pl-5">
                
                

             <div className="row">
                <div className="col-lg-8 item pt-4 pb-4">
                    <h1 ><strong> Features </strong></h1>
                    <hr className="style1"/>

                     <ul>
                        <li>
                           <p> Query non-coding cancer drivers by gene name, cancer types
                              , element types or evidence types that include computational
                               prediction method, gene expression associations from
                                RNA-seq and other experimental methods.</p>
                        </li>

                        <li>
                           <p> Illustrate the query results in summary plots.</p>
                        </li>

                        <li>
                           <p> Prioritize non-coding cancer driver candidates for follow-up experiments.</p>
                        </li>

                        <li>
                           <p> Users can submit curated non-coding cancer driver list to be included
                               in the new database release.</p>
                        </li>

                     </ul>
                   

               </div>
            </div>



                


        </div>

    );
  }
}

export default Features;
