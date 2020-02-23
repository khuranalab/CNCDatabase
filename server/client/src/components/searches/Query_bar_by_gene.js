import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { bindActionCreators } from 'redux';
//import { fetchBygeneSymbol, fetchGeneSummary, fetchWeather } from '../../actions/index';
import Autosuggest from 'react-autosuggest';
import AutosuggestHighlightMatch from 'autosuggest-highlight/match';
import AutosuggestHighlightParse from 'autosuggest-highlight/parse';

// Imagine you have a list of languages that you'd like to autosuggest.
const people = [
    {
      first: '(GBM) - ',
      last: 'Glioblastoma Multiforme',
      twitter: 'dancounsell'
    },
    {
      first: '(BRCA) - ',
      last: 'Breast Invasive Carcinoma',
      twitter: 'mtnmissy'
    },
    {
      first: '(BLCA) - ',
      last: 'Bladder Urothelial Carcinoma',
      twitter: 'ladylexy'
    },
    {
      first: '(ESAD) - ',
      last: 'Esophageal Carcinoma',
      twitter: 'steveodom'
    }
];




// Teach Autosuggest how to calculate suggestions for any given input value.
//https://developer.mozilla.org/en/docs/Web/JavaScript/Guide/Regular_Expressions#Using_Special_Characters
function escapeRegexCharacters(str) {
  return str.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

function getSuggestions(value) {
  const escapedValue = escapeRegexCharacters(value.trim());
  
  if (escapedValue === '') {
    return [];
  }

  const regex = new RegExp('\\b' + escapedValue, 'i');
  
  return people.filter(person => regex.test(getSuggestionValue(person)));
}
  
  // When suggestion is clicked, Autosuggest needs to populate the input
  // based on the clicked suggestion. Teach Autosuggest how to calculate the
  // input value for every given suggestion.
  function getSuggestionValue(suggestion) {
    return `${suggestion.first} ${suggestion.last}`;
  }
  
  // Use your imagination to render suggestions.
  function renderSuggestion(suggestion, { query }) {
    const suggestionText = `${suggestion.first} ${suggestion.last}`;
    const matches = AutosuggestHighlightMatch(suggestionText, query);
    const parts = AutosuggestHighlightParse(suggestionText, matches);
  
    return (
      <span className={'suggestion-content ' + suggestion.twitter}>
        <span className="name">
          {
            parts.map((part, index) => {
              const className = part.highlight ? 'text-danger' : null;
  
              return (
                <span className={className} key={index}>{part.text}</span>
              );
            })
          }
        </span>
      </span>
    );
  }
  

  

  class QueryBarByGene extends Component {
    constructor() {
      super();
  
      // Autosuggest is a controlled component.
      // This means that you need to provide an input value
      // and an onChange handler that updates this value (see below).
      // Suggestions also need to be provided to the Autosuggest,
      // and they are initially empty because the Autosuggest is closed.
      this.state = {
        value: '',
        suggestions: []
      };
    }
  
    onChange = (event, { newValue, method }) => {
      this.setState({
        value: newValue
      });
    };
  
    // Autosuggest will call this function every time you need to update suggestions.
    // You already implemented this logic above, so just use it.
    onSuggestionsFetchRequested = ({ value }) => {
      this.setState({
        suggestions: getSuggestions(value)
      });
    };
  
    // Autosuggest will call this function every time you need to clear suggestions.
    onSuggestionsClearRequested = () => {
      this.setState({
        suggestions: []
      });
    };
  
    render() {
      const { value, suggestions } = this.state;
  
      // Autosuggest will pass through all these props to the input.
      const inputProps = {
        placeholder: 'cancer type: GBM',
        value,
        onChange: this.onChange
      };
  
      const theme = {
       container: 'autosuggest',
       input: 'form-control',
       suggestionsContainer: 'dropdown',
       suggestionsList: `dropdown-menu ${suggestions.length ? 'show' : ''}`,
       suggestion: 'dropdown-item',
       suggestionHighlighted: 'active'
     };



      // Finally, render it!
      return (
        <Autosuggest
          theme={theme}
          suggestions={suggestions}
          onSuggestionsFetchRequested={this.onSuggestionsFetchRequested}
          onSuggestionsClearRequested={this.onSuggestionsClearRequested}
          getSuggestionValue={getSuggestionValue}
          renderSuggestion={renderSuggestion}
          inputProps={inputProps}
        />
      );
    }
  }

export default QueryBarByGene;

//export default SearchBar;
