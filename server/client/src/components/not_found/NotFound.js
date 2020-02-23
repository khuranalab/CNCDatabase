import React from 'react';
import { Link } from 'react-router-dom';
//import PageNotFound from '../assets/images/PageNotFound';



const NotFound = () => (

  
    <div className="py-5">

    <center>

    <span className="h1">Error 404</span>
    <br/>
    <span className="h5">- page not found</span>
    <br/>
      <div className="py-3">
        <Link to="/">
        <button className="btn btn-warning">
        Back to Home
        </button>
        </Link>
      </div>  
    </center>
    

    </div>
  
);

export default NotFound;