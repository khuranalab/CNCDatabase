import React, { Component } from 'react';
import { connect } from 'react-redux';
import { Link } from 'react-router-dom';

class Header extends Component {
  render() {
    return (
      <nav className="grey lighten-5 z-depth-0">
        <div className="nav-wrapper">
          <Link to={'/'} className="left brand-logo black-text">
            <i class="material-icons">settings_system_daydream</i>
            CNCDatabase
          </Link>
          <ul id="nav-mobile" className="right hide-on-med-and-down">
            <li>
              <a href="#" className="black-text">
                Home
              </a>
            </li>

            <li>
              <a href="/search" className="black-text">
                Search
              </a>
            </li>
            <li>
              <a href="#" className="black-text">
                Download
              </a>
            </li>
            <li>
              <a href="#" className="black-text">
                Help
              </a>
            </li>
            <li>
              <a href="#" className="black-text">
                Contact
              </a>
            </li>
          </ul>
        </div>
      </nav>
    );
  }
}

export default Header;
