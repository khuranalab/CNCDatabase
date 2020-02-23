import React, { Component } from 'react';
//import { Field, reduxForm } from 'redux-form';
import { connect } from 'react-redux';
//import queryString from 'query-string';
import { bindActionCreators } from 'redux';
import { fetchBygeneSymbol, fetchGeneSummary, fetchWeather } from '../../actions/index';
import { fetchCancerDriverList } from '../../actions/index'
//import Autosuggest from 'react-autosuggest';
//import AutosuggestHighlightMatch from 'autosuggest-highlight/match';
//import AutosuggestHighlightParse from 'autosuggest-highlight/parse';
//import EvidenceTypeInput from './evidenceTypeSelect';
import QueryBarByGene from './Query_bar_by_gene';
import EvidenceTypeInput from './evidenceTypeSelect';


class QueryForm2 extends Component {
   
    


   renderInput({ input, label}){
      return(
          <div className="form-group">
            <label>{label}</label>
            <input className="form-control" {...input}/>
          </div>
      );
   }

   onSubmit = (formValues) => {
       console.log(formValues);
       this.props.fetchGeneSummary(formValues.geneSymbol);
       //console.log("test",stringify)
   }

   //onSubmit={this.props.handleSubmit(this.onSubmit)}
   render(){

        return(
             <form>
             {/*<form onSubmit={this.props.handleSubmit(this.onSubmit)}>*/}
                 <Field name="geneSymbol" component={this.renderInput} label="Enter gene name" />
                 <Field name="canerType" component={this.renderInput} label="Select cancer type" />
                 <Field name="evidenceType" component={this.renderInput} label="Select evidence type" />
                 <EvidenceTypeInput/>
                 <button type="submit" className="btn btn-primary">
                        Submit
                 </button>
             </form>     
        );
   }

}

QueryForm2 = connect(null, {fetchGeneSummary} )(QueryForm2);



