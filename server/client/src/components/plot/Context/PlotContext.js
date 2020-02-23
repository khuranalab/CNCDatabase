import React, { createContext, useState } from "react";


// react context API example
// https://www.codementor.io/@sambhavgore/an-example-use-context-and-hooks-to-share-state-between-different-components-sgop6lnrd

export const Context = createContext({});

export const Provider = props => {
  
  
  const {
     groupOptions: initialGroupOptions, 
     selectedFilters: initialSelectedFilters,
     children
  } = props;

  const [ groupOptions, setGroupOptions ] = useState(initialGroupOptions);
  const [ selectedFilters, setSelectedFilters ] = useState(initialSelectedFilters);

  const plotContext = {
     groupOptions,
     setGroupOptions,
     selectedFilters,
     setSelectedFilters
  };

  return (
     <Context.Provider value={plotContext}>
         {children}
     </Context.Provider>    
  );

};

export const { Consumer } = Context;