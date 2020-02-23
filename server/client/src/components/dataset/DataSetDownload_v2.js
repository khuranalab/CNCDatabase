import React, { Component } from 'react';


class DataSetDownload extends Component {
  render() {
    return (
      <div className="container-fluid py-5 pl-5">
        <div className="row">
           <div className="col-lg-10">
              <h1><strong>Download</strong></h1>
              <hr className="style1"/>

              <div className="card">
                <div className="card-header">
                   <h5><strong>PostgreSQL database content</strong></h5>
                </div>
                <div className="card-body">
                    
                    <span className="mr-2">PostgreSQL database file: </span>
                    <a download="cncdatabase_content_Jan_26_2020" target="_blank" href="https://github.com/khuranalab/CNCDatabase_data/raw/master/Jan_26_2020_release/cncdatabase_content_Jan_26_2020">file</a>
                    <hr/>
                    <span className="mr-2">Database schema design file (Entity-Relationship diagram): </span>
                    <a download="noncoding_cancer_driver_Jan_26_2020_w_COSMIC_cancer_gene_annotation.txt" target="_blank" href="https://github.com/khuranalab/CNCDatabase_data/raw/master/Jan_26_2020_release/noncoding_cancer_driver_Jan_26_2020_w_COSMIC_cancer_gene_annotation.txt">file</a>
                    <hr/>
                    <span className="mr-2">List of 25 publications: </span>
                    <a download="study.xlsx" target="_blank" href="https://github.com/khuranalab/CNCDatabase_data/raw/master/Jan_26_2020_release/study.xlsx">file</a>

                </div>
               </div>







              

           </div> 
        </div>    
      </div>
    );
  }
}

export default DataSetDownload;