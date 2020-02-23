import React, { Component } from 'react';
import { Formik } from 'formik';
//import * as Yup from 'yup';
import { connect } from 'react-redux';
//import queryString from 'query-string';
import { bindActionCreators } from 'redux';
import { fetchBygeneSymbol, fetchGeneSummary } from '../../actions/index';
//import QueryBarByGene from './Query_bar_by_gene';
import { fetchCancerDriverList } from '../../actions/index'
//import Autosuggest from 'react-autosuggest';
//import AutosuggestHighlightMatch from 'autosuggest-highlight/match';
//import AutosuggestHighlightParse from 'autosuggest-highlight/parse';
//import EvidenceTypeInput from './evidenceTypeSelect';
//import QueryBarByGene from './Query_bar_by_gene';
//import EvidenceTypeInput from './evidenceTypeSelect';

import MyComponent from './CustomSelect';
//import EvidenceTypeInput2 from './evidenceTypeSelect2';
//import RFReactSelect from './RFReactSelect';

class QueryFormV7 extends Component {


   render(){



       return(


        <Formik
            initialValues={{ geneSymbol: '', elementType:'', cancerType:'', evidenceType:''}}
            onSubmit={(values, { setSubmitting }) => {
                setTimeout(() => {
                //alert(JSON.stringify(values, null, 2));
                console.log("form submit values", values)
                //this.props.fetchBygeneSymbol(values.geneSymbol);

                this.props.fetchGeneSummary(values.geneSymbol);

                this.props.fetchCancerDriverList(values);

                setSubmitting(false);
                }, 5);
            }}
            
        >

      {props => {
        const {
          values,
          touched,
          errors,
          dirty,
          isSubmitting,
          handleChange,
          handleBlur,
          handleSubmit,
          handleReset,
        } = props;

        return (
          <form onSubmit={handleSubmit} >
            <div className="form-group">
                <label for="geneSymbol">Gene name</label>
                <input
                    type="text"
                    className="form-control"
                    id="geneSymbol"
                    placeholder="Enter a gene name"
                    value={values.geneSymbol}
                    onChange={handleChange}
                    onBlur={handleBlur}
                
                />

                {errors.geneSymbol && 
                touched.geneSymbol &&
                    (<div style={{ color: 'red', marginTop: '.5rem' }}>{errors.email}</div>)
                }
            </div>

            <div className="form-group">
                <label for="elementType">Element Type</label>
                <input
                    type="text"
                    className="form-control"
                    id="elementType"
                    placeholder="Enter a element type"
                    value={values.elementType}
                    onChange={handleChange}
                    onBlur={handleBlur}
                
                />

            </div>



            <div className="form-group">
                <label for="cancerType">Cancer Type</label>
                <input
                    type="text"
                    className="form-control"
                    id="cancerType"
                    placeholder="Enter a cancer type"
                    value={values.cancerType}
                    onChange={handleChange}
                    onBlur={handleBlur}
                
                />

                {errors.geneSymbol && 
                touched.geneSymbol &&
                    (<div style={{ color: 'red', marginTop: '.5rem' }}>{errors.email}</div>)
                }
            </div>

            <div className="form-group">
                <label for="evidenceType">Evidence Type</label>
                <input
                    type="text"
                    className="form-control"
                    id="evidenceType"
                    placeholder="Enter a evidence type"
                    value={values.evidenceType}
                    onChange={handleChange}
                    onBlur={handleBlur}
                
                />

            </div>


            <div>
             <MyComponent/> 
            </div>  


            <button type="submit" disabled={isSubmitting} className="btn btn-primary mr-2">
              Submit
            </button>

            <button
              type="reset"
              className="btn btn-secondary"
              onClick={handleReset}
              disabled={!dirty || isSubmitting}
            >
              Reset
            </button>

            

            
          </form>
        );
      }}
    </Formik>
       );
   }
}

//QueryForm = connect(null, {fetchGeneSummary, fetchBygeneSymbol} )(QueryForm);

//export default QueryFormV5;

function mapDispatchToProps(dispatch) {
    return bindActionCreators({ fetchBygeneSymbol, fetchGeneSummary, fetchCancerDriverList }, dispatch);
}
  
export default connect(null, mapDispatchToProps)(QueryFormV7);
