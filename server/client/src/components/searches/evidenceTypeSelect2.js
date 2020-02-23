import React, { Component } from 'react'
import Select from 'react-select'

const options = [
  { value: 'chocolate', label: 'Computational prediction' },
  { value: 'strawberry', label: 'Gene expression change' },
  { value: 'vanilla', label: 'Expreriment' }
]

class EvidenceTypeInput2 extends Component{
    state = {
      selectedOption: null,
    }
    
    handleChange = (selectedOption) => {
      this.setState({ selectedOption });
      this.props.input.onChange(selectedOption.value)
      console.log(`Option selected:`, selectedOption);
    }
      
   handleBlur = () =>{
      this.props.input.onBlur(this.props.input.value)
   }


    render() {
        const { selectedOption } = this.state;
        const { input, name, id, onBlur, onFocus, value, ...custom } = this.props;
    
        return (
         <Select
          id={id}
          name={name}
          defaultValue={options[1]}
          value={this.props.input.value}
          onChange={this.handleChange}
          onFocus={onFocus}
          options={options}
         />
        );
    }

}

export default EvidenceTypeInput2;