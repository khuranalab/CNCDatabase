import React, { Component } from 'react';
import { Field, reduxForm } from 'redux-form';
import { connect } from 'react-redux';
//import queryString from 'query-string';
//import { bindActionCreators } from 'redux';
import { fetchBygeneSymbol, fetchGeneSummary } from '../../actions/index';
//import { fetchCancerDriverList } from '../../actions/index'
//import Autosuggest from 'react-autosuggest';
//import AutosuggestHighlightMatch from 'autosuggest-highlight/match';
//import AutosuggestHighlightParse from 'autosuggest-highlight/parse';
//import EvidenceTypeInput from './evidenceTypeSelect';
//import QueryBarByGene from './Query_bar_by_gene';
//import EvidenceTypeInput from './evidenceTypeSelect';
//import EvidenceTypeInput2 from './evidenceTypeSelect2';
import RFReactSelect from './RFReactSelect';

class QueryForm extends Component {
   
    


   renderInput({ input, label}){
      return(
          <div className="form-group">
            <label>{label}</label>
            <input className="form-control" {...input}/>
          </div>
      );
   }

   onSubmit = (formValues) => {
       console.log("formValues",formValues);
       this.props.fetchGeneSummary(formValues.geneSymbol);
       this.props.fetchBygeneSymbol(formValues.geneSymbol);
       //console.log("test",stringify)
   }

   //onSubmit={this.props.handleSubmit(this.onSubmit)}
   render(){

       const myOptions = [
        { value: 'computational prediction', label: 'Computational prediction' },
        { value: 'gene expression association', label: 'Gene expression change' },
        { value: 'experimental validation', label: 'Expreriment' }
      ]



        return(
             
             <form onSubmit={this.props.handleSubmit(this.onSubmit)} className="form-group">
                 
                 <Field name="geneSymbol" component={this.renderInput} label="Enter Gene Name" />

                 <Field name="canerType" component={this.renderInput} label="Select Cancer Type" />

                 <label>Select Evidence Types</label>
                 <Field name="evidenceType" options={myOptions} component={RFReactSelect} label="Select evidence type" />
                 <button type="submit" className="btn btn-primary">
                        Submit
                 </button>
             </form>     
        );
   }

}

QueryForm = connect(null, {fetchGeneSummary, fetchBygeneSymbol} )(QueryForm);

export default reduxForm({
    form: 'queryCreate'
})(QueryForm);

