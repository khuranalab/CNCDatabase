import React, { useState, useEffect } from 'react';
import useForm from 'react-hook-form';
import Select from 'react-select';

import { connect } from 'react-redux';
import { bindActionCreators } from 'redux';
import { fetchBygeneSymbol, fetchGeneSummary } from '../../actions/index';
import { fetchCancerDriverList } from '../../actions/index'


const options = [
    { value: 'promoter', label: 'promoter' },
    { value: 'enhancer', label: 'enhancer' },
    { value: 'lncRNA', label: 'lncRNA' },
    { value: 'CTCF-cohesin insulator', label: 'CTCF-cohesin insulator' },
];

const options2 = [
    { value: 'BRCA', label: 'BRCA' },
    { value: 'LUAD', label: 'LUAD' },
    { value: 'GBM', label: 'GBM' },
    { value: 'AML', label: 'AML' },
];

function QueryHookForm() {
  const { register, handleSubmit, setValue } = useForm(); // initialise the hook

  const onSubmit = data => {
    console.log(data);
    alert(JSON.stringify(data, null));
    
    const elementType = data.reactSelect.map(({ value }) => value)
    console.log( elementType)


    //this.props.fetchGeneSummary(vdata.geneSymbol);
    this.props.fetchCancerDriverList(elementType);

  }; // callback when validation pass

  const [values, setReactSelectValue] = useState({ 
      selectedOption: []
  });

  const handleMultiChange = selectedOption => {
    setValue("reactSelect", selectedOption)
    setReactSelectValue({ selectedOption })
  }
  
  useEffect(() => {
    register({ name: "reactSelect" })
  }, [register])


  return (
    <form onSubmit={handleSubmit(onSubmit)}>

      <div className="form-group">
       <label for ="geneSymbol">Gene name</label>   
       <Select
         value={values.selectedOption}
         options={options}
         onChange={handleMultiChange}
         isMulti
       />
      </div> 

      <div className="form-group">
       <label for ="elementType">Element Type</label>   
       <Select
         value={values.selectedOption}
         options={options}
         onChange={handleMultiChange}
         isMulti
       />
      </div> 

      <div className="form-group">
       <label for ="cancerType">Cancer Type</label>   
       <Select
         value={values.selectedOption}
         options={options}
         onChange={handleMultiChange}
         isMulti
       />
      </div> 

      <div className="form-group">
       <label for ="evidenceType">Evidence Type</label>   
       <Select
         value={values.selectedOption}
         options={options}
         onChange={handleMultiChange}
         isMulti
       />
      </div> 

      <button type="submit" className="btn btn-primary mr-2">
              Submit
      </button>

      <button
              type="reset"
              className="btn btn-secondary"
              
      >
              Reset
      </button>

    </form>
  );
}

//export default QueryHookForm;

function mapDispatchToProps(dispatch) {
    return bindActionCreators({ fetchBygeneSymbol, fetchGeneSummary, fetchCancerDriverList }, dispatch);
}
  
export default connect(null, mapDispatchToProps)(QueryHookForm);