import React from "react";
//import { render } from "react-dom";
import { withFormik } from "formik";
import * as yup from "yup";
//import Select from "react-select";
//import Papa from "papaparse";

import { connect } from 'react-redux';
import { bindActionCreators } from 'redux';
import { createSubmission } from '../../actions/index';
//import { fetchBygeneSymbol, fetchGeneSummary } from '../../actions/index';

import { ToastContainer, toast } from 'react-toastify';
import 'react-toastify/dist/ReactToastify.css';

//import geneNameData from './data/functional_element_with_snv_by_cancer_merged_May_15_2019_unique_geneName.csv';
//import geneNameData from './data/test.csv';


import emailjs from "emailjs-com";
import ReCAPTCHA from "reaptcha";

const FILE_SIZE = 10000 * 1024;
const SUPPORTED_FORMATS = [
   "text/csv"
];


const formikEnhancer = withFormik({

    
    validationSchema: yup.object().shape({
      name: yup.string()
         .required("Name is required."), 

      email: yup.string()
         .email("Invalid email address format.")
         .required("Email is required."),

      affiliation: yup.string()
         .required("Affiliation is required."), 

      file: yup.mixed()
         .required("A file is required.")
         .test(
          "fileSize",
          "File size cannot be larger than 10MB. ",
          value => {
            //console.log("yup test",value);
            return value && value.size <= FILE_SIZE; 
          })
         .test(
            "fileFormat",
            "Please upload .csv format text file",
            value => value && SUPPORTED_FORMATS.includes(value.type)
          ),
         
          
      recaptcha: yup.string()
         .required("reCAPTCHA is required.")    
    }),

    mapPropsToValues: props => ({
      name: [],
      email: [],
      affiliation: [],
      referenceUrl: [],
      comments: [],
      file: "",
      recaptcha: ""
    }),

    handleSubmit: (values, { props, setSubmitting }) => {

      const payload = {
        ...values,
        fileName: values.file.name,
        fileType: values.file.type,
        fileSize: `${values.file.size} bytes`
        //name: values.name.map(t => t.value),
        //elementType: values.elementType.map(t => t.value),
        //cancerType: values.cancerType.map(t => t.value),
        //evidenceType: values.evidenceType.map(t => t.value)
      };

      let formData = new FormData();
      
      formData.append("name", values.name);
      formData.append("email", values.email);
      formData.append("affiliation", values.affiliation);
      formData.append("referenceUrl", values.referenceUrl);
      formData.append("comments", values.comments);

      //for (let i =0; i<= values.files.length; i++){
      //     formData.append(`files[${i}]`, values.files[i]);
      //}
    
      formData.append("file",values.file)
      formData.append("fileName",values.file.name)
      formData.append("fileType",values.file.type)

      //
      let templateParams = {
        from_name: values.name,
        from_email: values.email,
        affiliation: values.affiliation,
        referenceUrl: values.referenceUrl,
        to_name: "cncdatabase administrator",
        fileName: values.file.name,
        message_html: values.comments  
      }

      emailjs.send("gmail","template_LnZJVI6L",templateParams,"user_Ff5atFvFpn9DsQwoNFehl")
        .then((response) => {
           //console.log("email sent success",response.status,response.text);
           toast.success('email send to administrator success')

        }, (error) => {
           //console.log("email sent with error",error.text);
           toast.error('email send with error')
        });


      //

      setTimeout(() => {
        //alert(JSON.stringify(payload, null, 2));
        //console.log("query payload", payload);
        //props.fetchGeneSummary(payload.geneSymbol);
        props.createSubmission(formData);

        setSubmitting(false);
      }, 1000);

    },
    displayName: "MyForm"
  });
  
const FileUploadForm = props => {
    
    const {
      values,
      touched,
      dirty,
      errors,
      handleChange,
      handleBlur,
      handleSubmit,
      handleReset,
      setFieldValue,
      setFieldTouched,
      isSubmitting
    } = props;


    return (

      <div className="container-fluid py-5 pl-5">
       
       <div className="row"> 
        

        <div className="col-lg-9">  



          <p className="h1"><strong>Submission</strong></p>
          <hr className="style1"/> 

            <div className="card">
            <h5 className="card-header"><strong>Submit non-coding cancer driver list</strong></h5> 
              <div className="card-body">
                  <p className="pt-3">Submit your curated non-coding cancer driver list to be included in the next data release.</p>
                  <div className="card">
                     <div className="card-body">
                       Submitted "non-coding cancer driver list" will require to be a csv text file separared by "," containing these seven columns.  
                       <br/>
                       <strong>PMID, Cancer Type, Gene Name, Element, Cohort Size, Mutated Samples, Evidence Type, Evidence Method</strong>
                     </div> 
                  </div>  
                  <form onSubmit={handleSubmit}>
                    <div className="form-row pt-3">
                      <div className="form-group col-md-6">
                        <label for="name">
                          <strong>Name</strong>
                          <span className="text-danger"> * </span>
                          <span className="text-secondary"><i> (required) </i></span>
                        </label>
                        <input
                          id="name"
                          className="form-control"
                          placeholder="Enter submitter full name"
                          type="text"
                          value={values.name}
                          onChange={handleChange}
                          onBlur={handleBlur}
                        />

                        {errors.name && touched.name && (
                          <div style={{ color: "red", marginTop: ".5rem" }}>{errors.name}</div>
                        )}
                      </div>

                      <div className="form-group col-md-6">
                        <label for="email">
                          <strong>Email</strong>
                          <span className="text-danger"> * </span>
                          <span className="text-secondary"><i> (required) </i></span>
                          
                        </label>
                        <input
                          id="email"
                          className="form-control"
                          placeholder="Enter submitter email"
                          type="email"
                          value={values.email}
                          onChange={handleChange}
                          onBlur={handleBlur}
                        />

                        {errors.email && touched.email && (
                          <div style={{ color: "red", marginTop: ".5rem" }}>{errors.email}</div>
                        )}
                      </div>

                      <div className="form-group col-md-6">
                        <label for="affiliation">
                          <strong>Affiliation</strong>
                          <span className="text-danger"> * </span>
                          <span className="text-secondary"><i> (required) </i></span>
                        </label>
                        <input
                          id="affiliation"
                          className="form-control"
                          placeholder="Enter submitter affiliation"
                          type="text"
                          value={values.affiliation}
                          onChange={handleChange}
                          onBlur={handleBlur}
                        />

                        {errors.affiliation && touched.affiliation && (
                          <div style={{ color: "red", marginTop: ".5rem" }}>{errors.affiliation}</div>
                        )}
                      </div>

                      <div className="form-group col-md-6">
                        <label for="referenceUrl">
                          <strong>Reference URL</strong>
                          <span className="text-secondary"><i> (optional) </i></span>
                        </label>
                        <input
                          id="referenceUrl"
                          className="form-control"
                          placeholder="Enter reference url"
                          type="text"
                          value={values.referenceUrl}
                          onChange={handleChange}
                          onBlur={handleBlur}
                        />

                        {errors.referenceUrl && touched.referenceUrl && (
                          <div style={{ color: "red", marginTop: ".5rem" }}>{errors.referenceUrl}</div>
                        )}
                      </div>

                      <div className="form-group col-md-12">
                        <label for="comments">
                          <strong>Comments</strong>
                          <span className="text-secondary"><i> (optional) </i></span>
                        </label>
                        <textarea
                          id="comments"
                          className="form-control"
                          placeholder="Enter comments"
                          type="text"
                          rows="5"
                          value={values.comments}
                          onChange={handleChange}
                          onBlur={handleBlur}
                        />

                        {errors.comments && touched.comments && (
                          <div style={{ color: "red", marginTop: ".5rem" }}>{errors.comments}</div>
                        )}
                      </div>

                      <div className="form-group col-md-12">
                        <label for="file">
                          <strong>Upload Your File</strong>
                          <span className="text-danger"> * </span>
                          <span className="text-secondary"><i> (required) </i></span>
                        </label>
                        <input
                            id="file" 
                            type="file" 
                            className="form-control" 
                            onChange={
                              (event) => {setFieldValue("file", event.currentTarget.files[0]);}
                            }
                            onBlur={handleBlur}
                        />

                        {errors.file && touched.file && (
                          <div style={{ color: "red", marginTop: ".5rem" }}>{errors.file}</div>
                        )}

                      </div>
                    
                    </div> 

                    <ToastContainer
                    /> 

                    <ReCAPTCHA
                       sitekey="6Lf2_9cUAAAAAEbCpWuWp13cC5vJpwNM59PB9Qwp"
                       //onVerify={handleChange}
                       reset={handleChange}
                       onVerify={(response) => { setFieldValue("recaptcha", response); console.log("recaptcha", response)}}
                    />

                    {errors.recaptcha && touched.recaptcha && (
                          <div style={{ color: "red", marginTop: ".5rem" }}>{errors.recaptcha}</div>
                    )}

                    <div className="py-3">
                        <button type="submit" disabled={isSubmitting} className="btn btn-primary mr-2">
                          Submit
                        </button>
                  
                        <button
                          type="reset"
                          className="btn btn-secondary"
                          onClick={handleReset}
                          disabled={!dirty || isSubmitting}
                        >
                          Reset
                        </button>
                    </div>

                  </form>
              </div>
            </div>     
        </div>
        
        <div className="col-lg-3">
        </div> 
       </div> 
      </div> 

    );
};
  
  



const FileUploadEnhancedForm = formikEnhancer(FileUploadForm);
  

//const QueryEnhancedForm = formikEnhancer(QueryForm);

//export default QueryEnhancedForm;


function mapDispatchToProps(dispatch) {
  return bindActionCreators({ createSubmission }, dispatch);
}

export default connect(null, mapDispatchToProps)(FileUploadEnhancedForm);

//export default FileUploadEnhancedForm
