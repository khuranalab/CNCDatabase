import React, { Component } from 'react'
import Select from 'react-select'

const options = [
  { value: 'chocolate', label: 'Chocolate' },
  { value: 'strawberry', label: 'Strawberry' },
  { value: 'vanilla', label: 'Vanilla' }
]

class EvidenceTypeInput extends Component{
    state = {
      selectedOption: null,
    }
    
    handleChange = (selectedOption) => {
      this.setState({ selectedOption });
      console.log(`Option selected:`, selectedOption);
    }
      
    render() {
        const { selectedOption } = this.state;
        //const { input, name, id, ...custom } = this.props;
    
        return (
         <Select
          defaultValue={options[0]}
          label="original"
          value={selectedOption}
          onChange={this.handleChange}
          options={options}
         />
        );
    }

}

export default EvidenceTypeInput;