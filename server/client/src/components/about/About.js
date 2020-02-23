import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';

class About extends Component {
  render() {
    return (
      <div className="container col-lg-10 py-3">
        <div className="pt-3 pb-5">
        <h1>About</h1>
        <hr/>
        Cornell Non-coding database (CNCDatabase) is developed by <a href="http://khuranalab.med.cornell.edu/">Khurana lab</a> at <a href="https://weill.cornell.edu/">Weill Cornell Medicine</a>. 
        </div>

        <div className="pt-2 pb-5">
        <h3>Citations</h3>
        <hr/>
        Liu et al, CNCDatabase: Cornall non-coding database. Database (2019)
        </div>

        <div className="pt-2 pb-5">
        <h3>Project members</h3>
        <hr/>
        <h5><b>Weill Cornell Medicine</b></h5>
        <li>Eric Minwei Liu</li>
        <li>Alexander Martinez-Fundichely</li>
        <li>Rajesh Bollapragada</li>
        <li>Maurice Spiewack</li>
        <li>Ekta Khurana</li>
        </div>

        <div className="py-3">
        <h3>Funding</h3>
        <hr/>
        This project is supported by the NIH through AAAAA (AAAAA)
        </div>
        
        <div className="py-5">
        </div>
      </div>
    );
  }
}

export default About;