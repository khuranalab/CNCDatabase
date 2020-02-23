import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
//import { UndrawChoose, UndrawChecklist, UndrawWorkChat } from 'react-undraw-illustrations';
import introduction_figure from './img/introduction_v6.png';

class Intorduction extends Component {
  render() {
    return (
            <div className="container-fluid">
                
              <div className="row">  


                 <div className="col-lg-10 pt-3 pb-3 pl-5">
                    <h1 className="text-left"><strong> Introduction </strong></h1>
                    
                    <hr className="style1"/>

                    <p>   This database provides the list of non-coding
                        cancer drivers in gene promoters, enhancers, ncRNAs
                         and CTCF-cohesin insulators published in over 25 studies.
                          Most mutations from whole genome sequencing data occur
                           in the non-coding regions with unknown impact to the
                            tumor development. The advances in the study of
                             regulatory regions have been uncovering a large
                              number of putative non-coding cancer drivers.
                               We applied text mining and manual curation to
                                the PubMed articles and then developed
                                 CNCDatabase, Cornell Non-Coding Cancer driver
                                  Database. CNCDatabase provides a helpful
                                   resource for researchers to explore the
                                    pathological role of non-coding alterations
                                     and their associations with gene expression
                                      in human cancers.
                    </p>
                    
                   
                   

                 </div>

                 <div className="container d-flex px-1 pt-1 pb-5">
               
                
                    <center>
                     <img src={introduction_figure} alt="introduction" width="70%" height="100%"/>
                    </center>

                 </div>


              </div>  


            </div>
        
    );
  }
}

export default Intorduction;
