import React, { Component } from 'react';
import { reduxForm } from 'redux-form';
import { Link } from 'react-router-dom';

class SearchNew extends Component {
  render() {
    return (
      <div>
        <div className="fixed-action-btn">
          <Link to="/searchs/new" className="btn-floating btn-large red">
            <i className="material-icons">add</i>
          </Link>
        </div>
      </div>
    );
  }
}

export default SearchNew;
