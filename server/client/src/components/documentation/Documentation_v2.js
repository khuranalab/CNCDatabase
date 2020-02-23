import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
import { importMDX } from 'mdx.macro';
//import faqDocument from './faq.md';

//example to use react-markdown package
//https://stackoverflow.com/questions/42928530/how-do-i-load-a-markdown-file-into-a-react-component

const Content = importMDX.sync('./faq.md')


class Documentation2 extends Component {
  

  render() {
    
    

    return (
      <div className="container col-lg-10 py-3">
        <h1>Documentation</h1>
        <hr/>
        <Content
        />
      </div>
    );
  }
}

export default Documentation2;