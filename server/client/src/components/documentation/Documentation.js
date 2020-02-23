import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
import ReactMarkdown from 'react-markdown';
import faqDocument from './faq.md';

//example to use react-markdown package
//https://stackoverflow.com/questions/42928530/how-do-i-load-a-markdown-file-into-a-react-component

class Documentation extends Component {
  constructor(props) {
    super(props)

    this.state = { terms: null }
  }

  componentWillMount() {
    fetch(faqDocument).then((response) => response.text()).then((text) => {
      this.setState({ terms: text })
    })
  }
  
  

  render() {
    
    //https://github.com/rexxars/react-markdown/issues/69
    function flatten(text, child) {
      return typeof child === 'string'
        ? text + child
        : React.Children.toArray(child.props.children).reduce(flatten, text)
    }
    
    function HeadingRenderer(props) {
      var children = React.Children.toArray(props.children)
      var text = children.reduce(flatten, '')
      //console.log("text",text)

      // remove metacharacter such as "?" at the end of a string
      var text2 = text.replace(/\W$/,'')
      //console.log("text2",text2)
      var slug = text2.toLowerCase().replace(/\W/g, '-')
      //console.log("slug",slug)
      return React.createElement('h' + props.level, {id: slug}, props.children)
    }





    return (
      <div className="container col-lg-10 py-3">
        <h3>Documentation</h3>
        <hr/>
        <ReactMarkdown 
           skipHtml={false}
           escapeHtml={true}
           source={this.state.terms}
           renderers={{heading: HeadingRenderer}}
        />
      </div>
    );
  }
}

export default Documentation;