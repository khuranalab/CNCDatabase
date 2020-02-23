import React, { Component } from 'react';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { importMDX } from 'mdx.macro';
//import faqDocument from './faq.md';

//example to use react-markdown package
//https://stackoverflow.com/questions/42928530/how-do-i-load-a-markdown-file-into-a-react-component

//const Content = importMDX.sync('./faq.md')

//import ReactGoogleSlides from 'react-google-slides';


class Documentation extends Component {
  

  render() {
    
    

    return (
      <div className="container-fluid py-5 pl-5">
           <div className="row">
                
               <div className="col-lg-10">  
                    <h1><strong>Documentation</strong></h1>
                    <hr className="style1"/>
                    {/*<h4><strong>Table of contents</strong></h4>
                    <ul>
                        <li>What is CNCDatabase?</li>
                        <li>How do I cite the CNCDatabase?</li>
                        <li>How do I get update on new non-coding cancer driver list?</li>
                    </ul>*/}
                    <h4><strong>What is CNCDatabase?</strong></h4>
                    <hr/>
                    <p>Although there are many existing databases that provide the list of cancer
                       drivers at protein coding regions in the genome, there is no central
                        database dedicated to collection of information of cancer drivers
                         at non-coding regions. CNCDatabase aims to collect non-coding cancer
                          driver lists from publications and to complement the understanding
                           of cancer beyond protein coding regions.</p>
                    <br/>
                    <h4><strong>How do I cite the CNCDatabase?</strong></h4>
                    <hr/>
                    <p>Please cite Liu et al, CNCDatabase: Cornell Non-coding Cancer driver database (2020)</p>

                    <br/>
                    <h4><strong>Acknowledgements</strong></h4>
                    <hr/>
                    <p>
                      We thank Douglas Duckworth from Weill Cornell Supercomputing Unit for helps in the
                       website deployment. We also thank Ann Palladino, Sandra Cohen, Tawny Cuykendall,
                        Erica Duo Xu, Andre Forbes at the Khurana lab for valuable suggestions.     
                    </p>
                    {/*<div>
                        <h4><strong>Tutorial 1: example queries</strong></h4>
                        <iframe src="https://docs.google.com/presentation/d/e/2PACX-1vTHblNfskQLj6NzNSE-9BT0XckLc_T7NqEC31OpBjvnRkxW_sE68FHnFgicnGDd7BYm7rq7zdAkGrwu/embed?start=false&loop=false&delayms=60000" frameborder="0" width="960" height="569" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>
                    </div>*/}
               </div>

               <div className="col-lg-2">
               </div>         
          </div>
      </div>
    );
  }
}

export default Documentation;