import React, { Component } from 'react'
import Select from 'react-select'

const options = [
  { value: 'promoter', label: 'promoter' },
  { value: 'enhancer', label: 'enhancer' },
  { value: 'lncRNA', label: 'lncRNA' },
  { value: 'CTCF-cohesin insulator', label: 'CTCF-cohesin insulator' },
]

const MyComponent = () => (
  <Select 
      options={options} 
  />
)

export default MyComponent