import React, { useState, useEffect } from "react";
import ReactDom from "react-dom";
import useForm from "react-hook-form"
import Select from "react-select";

const options = [
    { value: "chocolate", label: "Chocolate" },
    { value: "strawberry", label: "Strawberry" },
    { value: "vanilla", label: "Vanilla" }
  ];

function hookForm() {

    const { register, handleSubmit, setValue } = useForm();
    
    const onSubmit = data => {
      alert(JSON.stringify(data, null));
    };
    
    return (
      <div className="App">
        <form onSubmit={handleSubmit(onSubmit)}>
          
          <div>
            <lable className="reactSelectLabel">React select</lable>
            <Select
              className="reactSelect"
              name="filters"
              placeholder="Filters"
              value={values.selectedOption}
              options={options}
              onChange={handleMultiChange}
              isMulti
            />
          </div>
  
          <button type="submit">Submit</button>
        </form>
      </div>
    );
  }

export default hookForm;