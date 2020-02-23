import React, { Component } from 'react';
import { connect } from 'react-redux';
import { bindActionCreators } from 'redux';
import { fetchBygeneSymbol, fetchGeneSummary, fetchWeather } from '../../actions/index';
import EvidenceTypeInput from './evidenceTypeSelect';
import EvidenceTypeInput2 from './evidenceTypeSelect2';

class QueryForm2 extends Component {
  constructor(props) {
    super(props);

    this.state = { term: '' };

    this.onInputChange = this.onInputChange.bind(this);
    this.onFormSubmit = this.onFormSubmit.bind(this);
  }

  onInputChange(event) {
    this.setState({ term: event.target.value });
    console.log(event.target.value);
  }

  onFormSubmit(event) {
    event.preventDefault();
    
    console.log("watch",event.target);
    // We need to go and fetch data
    this.props.fetchBygeneSymbol(this.state.term);
    this.props.fetchGeneSummary(this.state.term);
    //this.props.fetchWeather(this.state.term);
    //this.setState({ term: '' });
  }

  render() {
    return (
      <form onSubmit={this.onFormSubmit} className="form-group">
        <input
          placeholder="cancer type:GBM"
          className="form-control"
          value={this.state.term}
          onChange={this.onInputChange}
        />
        <EvidenceTypeInput/>
        <button type="submit" className="btn btn-primary">
           Submit
        </button>
        
      </form>
    );
  }
}

function mapDispatchToProps(dispatch) {
  return bindActionCreators({ fetchBygeneSymbol, fetchGeneSummary, fetchWeather }, dispatch);
}

export default connect(null, mapDispatchToProps)(QueryForm2);

//export default SearchBar;
