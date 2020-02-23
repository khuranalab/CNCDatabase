import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
import { importMDX } from 'mdx.macro';

const Content = importMDX.sync('./dataset.md')


class DataSetDownload extends Component {
  render() {
    return (
      <div className="container-fluid py-5 pl-5">
        <div className="row">
           <div className="col-lg-10">
              <h1><strong>Download</strong></h1>
              <hr/>
              <Content
              />
           </div> 
        </div>    
      </div>
    );
  }
}

export default DataSetDownload;