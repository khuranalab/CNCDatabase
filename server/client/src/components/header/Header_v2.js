import React, { Component } from 'react';
import { connect } from 'react-redux';
import { Link } from 'react-router-dom';
//import { Button } from 'reactstrap';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import ModalExample2 from './LoginModal_v2';
//# Home, Search, Download, Help, Contact
import './css/header_setup.css';
import NavLink from './NavLink';
import logo from './img/cncdatabase_logo_v2.png';






class Header extends Component {


    constructor(props) {
        super(props);
        this.state = {
          links: [
            {path: "/", text: "Home", isActive: true},
            {path: "/search", text: "Search", isActive: false},
            {path: "/dataset", text: "Download", isActive: false},
            {path: "/submission", text: "Submission", isActive: false},
            {path: "/documentation", text: "Documentation", isActive: false},
          ]
        }
    }
    
            
    handleClick(i) {
        const links = this.state.links.slice(); 
        for (const j in links) {
          links[j].isActive = i == j ;
        }
        this.setState({links: links});
    }





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
      // py-md-1  
      <nav className=" d-flex flex-row navbar navbar-expand-lg navbar-custom">
        
        <Link to="/" onClick={() => {this.handleClick(0)}}>
           <img src={logo} alt="cncdatabase" width="40" height="30"/>
        </Link>

        <Link className="navbar-brand" to="/" onClick={() => {this.handleClick(0)}}>
           <strong>CNCDatabase</strong>
        </Link>


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
             {this.state.links.map((link, i) => 
               <NavLink
                 path={link.path} 
                 text={link.text} 
                 isActive={link.isActive}
                 key={link.path} 
                 onClick={() => {this.handleClick(i)}}
              /> 
              )}
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