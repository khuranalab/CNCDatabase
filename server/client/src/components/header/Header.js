import React, { Component } from 'react';
import { connect } from 'react-redux';
import { Link } from 'react-router-dom';
import { Button } from 'reactstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import ModalExample2 from './LoginModal_v2';
//# Home, Search, Download, Help, Contact

class Header extends Component {


   





  renderContent(){
     switch (this.props.auth){
         case null:
           // user not login yet
           return;
         case false:
           // user already logged out
           return <li><a href="/auth/google">Login</a></li>;
           //return 'already logedout';
         default:
           // user already login
           return <li><a href="/api/logout">Logout</a></li>;
           //return 
     }
  }



  render() {
    
    //console.log(this.props);

    return (
      <nav className="navbar navbar-expand-lg navbar-dark bg-dark py-md-3">
        <FontAwesomeIcon icon="laptop-code" color="white" size="2x"/>
        <a className="navbar-brand pl-2" href="/">
          CNCDatabase
        </a>

        <button
          className="navbar-toggler"
          type="button"
          data-toggle="collapse"
          data-target="#navbarNavAltMarkup"
          aria-controls="navbarNavAltMarkup"
          aria-expanded="false"
          aria-label="Toggle navigation"
        >
          <span className="navbar-toggler-icon" />
        </button>
        <div className="collapse navbar-collapse" id="navbarNav">
          <ul className="navbar-nav">
            <li className="nav-item"><a className="nav-link active" href="/">Home <span className="sr-only">(current)</span></a></li>
            <li className="nav-item"><a className="nav-link" href="/search">Search</a></li>
            <li className="nav-item"><a className="nav-link" href="/dataset">Download</a></li>
            <li className="nav-item"><a className="nav-link" href="/submission">Submission</a></li>
            <li className="nav-item"><a className="nav-link" href="/documentation">Documentation</a></li>
            
          </ul>
          
         
        </div>

         {/* {this.renderContent()} */}

        <div>
         {/* <ModalExample2/> */}
        </div> 

      </nav>
    );
  }
}

function mapStateToProps( state){
   return { auth: state.auth };
}

//export default Header;
export default connect(mapStateToProps)(Header);

//<form className="form-inline">
//         <button className="btn btn-sm btn-outline-primary" type="button">Login</button>
//         </form>