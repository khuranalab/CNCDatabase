import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';

import WCM_Logo from './img/WCM_MB_LOGO_HZSS2L_CLR_RGB.png';

class Footer extends Component {
  render() {
    return (
        <div className="container-fluid py-3 pl-5">
          <div className="row">
            <div className="col-md-4 item">
              <h5><strong>Cornell Non-Coding Cancer Driver Database</strong></h5>
              
              <img src={WCM_Logo} alt="wcm_logo" width="40%" height="35%"/>
              
            </div>

            {/*<div className="col-md-2 item">
              <h5><strong>Navigation</strong></h5>
              <p>Home</p>
              <p>Search</p>
              <p>Download</p>
              <p>Submission</p>
              <p>Documentation</p>
              
             </div>*/}

            <div className="col-md-4 item align-bottom">
              <h5><strong>Contact</strong></h5>
              <p>
                
                <div className="row">
                  <a href="mailto:alm2069@med.cornell.edu;emliu.research@gmail.com?cc=ekk2003@med.cornell.edu&subject=cncdatabase suggestions">
                    <button className="btn btn-secondary ml-3">
                      Email suggestions
                    </button>
                  </a>
                  
                  <a target="_blank" href="https://khuranalab.med.cornell.edu/">
                    <button className="btn btn-secondary ml-3">
                      Khurana Lab
                    </button>
                  </a>
                </div>


              </p>  

              <p>Meyer Cancer Center, Weill Cornell Medicine, New York, New York 10065, USA</p>
              

            </div>

            <div className="col-md-2 item">
            </div>

          </div>

         <div className="row">
           <div className="col-lg-12">
              <p className="copyright">
                © 2018 - 2020, Weill Cornell Medicine. All Rights Reserved
              </p>
           </div> 
         </div> 


          

        </div>

    );
  }
}

export default Footer;
