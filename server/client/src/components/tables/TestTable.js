import React, { Component } from 'react';
import { connect } from 'react-redux';
import { fetchBygeneSymbol } from '../../actions';

class TestTable extends Component {
  componentDidMount() {
    this.props.fetchBygeneSymbol();
  }

  renderDataBygeneSmbol(dataBygeneSymbol) {
    return (
      <div className="card card-block">
        <h4 className="card-title">{dataBygeneSymbol.geneName}</h4>
        <p className="card-text">Chees Factory</p>
        <a className="btn btn-primary">email</a>
      </div>
    );
  }

  render() {
    console.log(this.props.dataBygeneSymbol);

    return (
      //<div>"test"</div>
      <div>{this.props.dataBygeneSymbol.map(this.renderDataBygeneSmbol)}</div>
    );
  }
}

// get state from redux store of reducers ( combineReducers ) and return object
// key used in the combineReducers 
function mapStateToProps({ dataBygeneSymbol }) {
  return { dataBygeneSymbol };
}

// pass action_creator (fetchBygenSymbol object) to Component as props
// connect calls dispatch function to put action into redux store
// short notion { fectchBygeneSymbol } = { fetchBygeneSymbol: fectchBygeneSymbol } 
export default connect(mapStateToProps, { fetchBygeneSymbol })(TestTable);

//export default TestTable;
