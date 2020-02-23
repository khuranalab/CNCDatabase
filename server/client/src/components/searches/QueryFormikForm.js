import React from "react";
//import { render } from "react-dom";
import { withFormik } from "formik";
import * as Yup from "yup";
import Select from "react-select";
import Papa from "papaparse";

import { connect } from 'react-redux';
import { bindActionCreators } from 'redux';
import { fetchCancerDriverList } from '../../actions/index';
import { fetchBygeneSymbol, fetchGeneSummary } from '../../actions/index';

import geneNameData from './data/functional_element_with_snv_by_cancer_merged_May_15_2019_unique_geneName.csv';
//import geneNameData from './data/test.csv';









const formikEnhancer = withFormik({


    mapPropsToValues: props => ({
      geneSymbol: [],
      elementType: [],
      cancerType: [],
      evidenceType: [],
      pmid: []
    }),

    handleSubmit: (values, { props, setSubmitting }) => {
      const payload = {
        ...values,
        geneSymbol: values.geneSymbol.map(t => t.value),
        elementType: values.elementType.map(t => t.value),
        cancerType: values.cancerType.map(t => t.value),
        evidenceType: values.evidenceType.map(t => t.value),
        pmid: values.pmid.map(t => t.value),
      };
      setTimeout(() => {
        //alert(JSON.stringify(payload, null, 2));
        console.log("query payload", payload);
        props.fetchGeneSummary(payload.geneSymbol);
        props.fetchCancerDriverList(payload);

        setSubmitting(false);
      }, 1000);
    },
    displayName: "MyForm"
  });
  
const QueryForm = props => {
    
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
      <form onSubmit={handleSubmit}>
         <div className="form-row">
            <div className="form-group col-md-3">
              <GeneSymbolSelect
                className="form-control"
                value={values.geneSymbol}
                onChange={setFieldValue}
                onBlur={setFieldTouched}
                error={errors.geneSymbol}
                touched={touched.geneSymbol}
              />
            </div>

            <div className="form-group col-md-3">
              <ElementTypeSelect
                className="form-control"
                value={values.elementType}
                onChange={setFieldValue}
                onBlur={setFieldTouched}
                error={errors.elementType}
                touched={touched.elementType}
              />
            </div>
           

            <div className="form-group col-md-3">
              <CancerTypeSelect
                className="form-control"
                value={values.cancerType}
                onChange={setFieldValue}
                onBlur={setFieldTouched}
                error={errors.cancerType}
                touched={touched.cancerType}
              />
            </div>

            <div className="form-group col-md-3">
              <EvidenceTypeSelect
                className="form-control"
                value={values.evidenceType}
                onChange={setFieldValue}
                onBlur={setFieldTouched}
                error={errors.evidenceType}
                touched={touched.evidenceType}
              />
            </div>
        </div>

        <div className="form-row">
            <div className="form-group col-md-3">
              <PublicationSelect
                className="form-control"
                value={values.pmid}
                onChange={setFieldValue}
                onBlur={setFieldTouched}
                error={errors.pmid}
                touched={touched.pmid}
              />
            </div>
        </div>     

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


      </form>
    );
};
  
  const testData="aa,bb";
  console.log("testData string",testData)

  const testString = Papa.parse(testData)
   
  console.log("papa parsed string",testString)                   

  let geneSymbolOptions3;
  Papa.parse(geneNameData, {
                                header: true,
                                delimiter: ",",
                                download: true,
                                dynamicTyping: true,
                                complete: function(results) {
                                console.log("Parsing complete", results);
                                geneSymbolOptions3 = results.data
                                }
                            });
  console.log("papa parsed file", geneSymbolOptions3)
   

  
  const geneSymbolOptions = [
    { value: 'FOXA1', label: 'FOXA1' },
    { value: 'TERT', label: 'TERT' },
    { value: 'lncRNA', label: 'lncRNA' },
    { value: 'CTCF-cohesin insulator', label: 'CTCF-cohesin insulator' },
  ];



  const elementTypeOptions = [
    { value: 'promoter', label: 'promoter' },
    { value: '5UTR', label: '5\'UTR' },
    { value: '3UTR', label: '3\'UTR' },
    { value: 'enhancer', label: 'enhancer' },
    { value: 'lncRNA', label: 'lncRNA' },
    { value: 'lncRNA promoter', label: 'lncRNA promoter' },
    { value: 'CTCF-cohesin insulator', label: 'CTCF-cohesin insulator' },
    { value: 'splice site', label: 'splice site' },
    { value: 'small RNA', label: 'small RNA' },
    { value: 'rRNA', label: 'rRNA' },
    { value: 'misc_RNA', label: 'misc_RNA' },
    { value: 'snRNA', label: 'snRNA' },
    { value: 'miRNA', label: 'miRNA' },
  ];
  
  const evidenceTypeOptions = [
    { value: 'computational prediction', label: 'computational prediction' },
    { value: 'gene expression association', label: 'RNA-seq gene expression association' },
    { value: 'experimental validation', label: 'other experimental validation' },
  ];


  const cancerTypeOptions = [
    { value: 'PanCancer', label: 'PanCancer' },
    { value: 'Adenocarcinoma', label: 'Adenocarcinoma' },
    { value: 'Breast', label: 'Breast' },
    { value: 'Carcinoma', label: 'Adenocarcinoma' },
    { value: 'CNS', label: 'Central nerve system' },
    { value: 'Digestive tract', label: 'Digestive tract' },
    { value: 'Female reproductive system', label: 'Female reproductive system' },
    { value: 'Giloma', label: 'Giloma' },
    { value: 'Hematopoietic system', label: 'Hematopoietic system' },
    { value: 'Kidney', label: 'Kidney' },
    { value: 'Lymphatic system', label: 'Lymphatic system' },
    { value: 'Sarcoma', label: 'Sarcoma' },
    { value: 'Squamous', label: 'Squamous' },
    { value: 'B-ALL', label: 'B-ALL, B-Cell Acute Lymphoblastic Leukemia' },
    { value: 'BLCA', label: 'BLCA, Bladder Urothelial Carcinoma' },
    { value: 'BRCA', label: 'BRCA, Breast Invasive Carcinoma' },
    { value: 'CLL', label: 'CLL, Chronic Lymphoctytic Leukemia' },
    { value: 'COAD', label: 'COAD, Colon Adenocarcinoma' },
    { value: 'COADREAD', label: 'COADREAD, Colorectal Cancer' },
    { value: 'DLBC', label: 'DLBC, Lymphoid Neoplasm Diffuse Large B-cell Lymphoma' },
    { value: 'ESCA', label: 'ESCA, Esophageal Carcinoma' },
    { value: 'ETMR', label: 'ETMR, Embryonal Tumor with Multilayered Rosettes' },
    { value: 'GBM', label: 'GBM, Glioblastoma Multiforme' },
    { value: 'HGG', label: 'HGG, High-grade Glioma' },
    { value: 'HNSC', label: 'HNSC, Head and Neck Squamous Cell Carcinoma' },
    { value: 'KICH', label: 'KICH, Kidney Chromophobe' },
    { value: 'KIRC', label: 'KIRC, Kidney Renal Clear Cell Carcinoma' },
    { value: 'LAML', label: 'LAML, Acute Myeloid Leukemia' },
    { value: 'LGG', label: 'LGG, Lower Grade Glioma' },
    { value: 'LIHC', label: 'LIHC, Liver Hepatocellular Carcinoma' },
    { value: 'LUAD', label: 'LUAD, Lung Adenocarcinoma' },
    { value: 'LUSC', label: 'LUSC, Lung Squamous Cell Carcinoma' },
    { value: 'Lymph-BNHL', label: 'Lymph-BNHL, B-cell non-Hodgkin lymphoma' },
    { value: 'MALY', label: 'MALY, Malignant Lymphoma' },
    { value: 'MB', label: 'MB, Medullablastoma' },
    { value: 'OV', label: 'OV, Ovarian Serous Cystadenocarcinoma' },
    { value: 'PA', label: 'PA, Pilocytic Astrocytoma' },
    { value: 'PACA', label: 'PACA, Pancreatic Ductal Adenocarcinoma' },
    { value: 'Pan-Endocrine', label: 'Pan-Endocrine, Pancreatic Endocrine Neoplasmas' },
    { value: 'PRAD', label: 'PRAD, Prostate Adenocarcinoma' },
    { value: 'RECA', label: 'RECA, Kidney Renal Papillary Cell Carcinoma' },
    { value: 'RHAB', label: 'RHAB, Rhabdoid Tumor' },
    { value: 'SKCM', label: 'SKCM, Skin Cutaneous Melanoma' },
    { value: 'STAD', label: 'STAD, Stomach Adenocarcinoma' },
    { value: 'THCA', label: 'THCA, Papillary Thyroid Carcinoma' },
    { value: 'UCEC', label: 'UCEC, Uterine Corpus Endometrial Carcinoma' },
    

  ];

  const publicationOptions = [
    { value: '23348506', label: 'Huang et al., Science, 2013 (PMID: 23348506)' },
    { value: '24092746', label: 'Khurana et al., Science, 2013 (PMID: 24092746)' },
    { value: '25383969', label: 'Fredriksson et al., Nature Genetics, 2014 (PMID: 25383969)' },
    { value: '25261935', label: 'Weinhold et al., Nature Genetics, 2014 (PMID: 25261935)' },
    { value: '26304545', label: 'Lochovsky et al., Nucleic Acids Research, 2015 (PMID: 26304545)' },
    { value: '26053494', label: 'Melton et al., Nature Genetics, 2015 (PMID: 26053494)' },
    { value: '26200345', label: 'Puente et al., Nature, 2015 (PMID: 26200345)' },
    { value: '26691984', label: 'Araya et al., Nature Genetics, 2016 (PMID: 26691984)' },
    { value: '27571262', label: 'Bailey et al., Nature Genetics, 2016 (PMID: 27571262)' },
    { value: '27064257', label: 'Fujimoto et al., Nature Genetics, 2016 (PMID: 27064257)' },
    { value: '27311963', label: 'Mularonl et al., Genome Biology, 2016 (PMID: 27311963)' },
    { value: '27135926', label: 'Nik-Zainal et al., Nature, 2016 (PMID: 27135926)' },
    { value: '28467829', label: 'Hayward et al., Nature, 2017 (PMID: 28467829)' },
    { value: '28362259', label: 'Juul et al., eLIFE, 2017 (PMID: 28362259)' },
    { value: '28128360', label: 'Lanzos et al., Scientific Reports, 2017 (PMID: 28128360)' },
    { value: '28658208', label: 'Rheinbay et al., Nature, 2017 (PMID: 28658208)' },
    { value: '29489754', label: 'Grobner et al., Nature, 2018 (PMID: 29489754)' },
    { value: '29662167', label: 'Wedge et al., Nature Genetics, 2018 (PMID: 29662167)' },
    { value: '30718927', label: 'Frankell et al., Nature Genetics, 2019 (PMID: 30718927)' },
    { value: '31078526', label: 'Liu et al., Cell Systems, 2019 (PMID: 31078526)' },
    { value: '31597163', label: 'Shuai et al., Nature, 2019 (PMID: 31597163)' },
    { value: '31597162', label: 'Suzuki et al., Nature, 2019 (PMID: 31597162)' },
    { value: '32025015', label: 'Rheinbay et al., Nature, 2020 (PMID: 32025015)' },
    { value: '31974375', label: 'Zhou et al., Nature Communications, 2020 (PMID: 31974375)' },
    { value: '31954095', label: 'Zhu et al., Molecular Cell, 2020 (PMID: 3195409)' },
  ];






  class GeneSymbolSelect extends React.Component {

    handleChange = value => {
      // this is going to call setFieldValue and manually update values.topcis
      this.props.onChange("geneSymbol", value);
    };
  
    handleBlur = () => {
      // this is going to call setFieldTouched and manually update touched.topcis
      this.props.onBlur("geneSymbol", true);
    };
  
    render() {
      return (
        <div style={{ margin: "0.25rem 0" }}>
          <label for="geneSymbol">Gene name </label>
          <Select
            id="geneSymbol"
            options={geneSymbolOptions3}
            isMulti
            onChange={this.handleChange}
            onBlur={this.handleBlur}
            value={this.props.value}
          />
          {!!this.props.error && this.props.touched && (
            <div style={{ color: "red", marginTop: ".5rem" }}>
              {this.props.error}
            </div>
          )}
        </div>
      );
    }
}



  class ElementTypeSelect extends React.Component {

    handleChange = value => {
      // this is going to call setFieldValue and manually update values.topcis
      this.props.onChange("elementType", value);
    };
  
    handleBlur = () => {
      // this is going to call setFieldTouched and manually update touched.topcis
      this.props.onBlur("elementType", true);
    };
  
    render() {
      return (
        <div style={{ margin: "0.25rem 0" }}>
          <label for="elementType">Element type </label>
          <Select
            id="elementType"
            options={elementTypeOptions}
            isMulti
            onChange={this.handleChange}
            onBlur={this.handleBlur}
            value={this.props.value}
          />
          {!!this.props.error && this.props.touched && (
            <div style={{ color: "red", marginTop: ".5rem" }}>
              {this.props.error}
            </div>
          )}
        </div>
      );
    }
}

class CancerTypeSelect extends React.Component {

    handleChange = value => {
      // this is going to call setFieldValue and manually update values.topcis
      this.props.onChange("cancerType", value);
    };
  
    handleBlur = () => {
      // this is going to call setFieldTouched and manually update touched.topcis
      this.props.onBlur("cancerType", true);
    };
  
    render() {
      return (
        <div style={{ margin: "0.25rem 0" }}>
          <label for="cancerType">Cancer type </label>
          <Select
            id="cancerType"
            options={cancerTypeOptions}
            isMulti
            onChange={this.handleChange}
            onBlur={this.handleBlur}
            value={this.props.value}
          />
          {!!this.props.error && this.props.touched && (
            <div style={{ color: "red", marginTop: ".5rem" }}>
              {this.props.error}
            </div>
          )}
        </div>
      );
    }
}

class EvidenceTypeSelect extends React.Component {

  handleChange = value => {
    // this is going to call setFieldValue and manually update values.topcis
    this.props.onChange("evidenceType", value);
  };

  handleBlur = () => {
    // this is going to call setFieldTouched and manually update touched.topcis
    this.props.onBlur("evidenceType", true);
  };

  render() {
    return (
      <div style={{ margin: "0.25rem 0" }}>
        <label for="evidenceType">Evidence type </label>
        <Select
          id="evidenceType"
          options={evidenceTypeOptions}
          isMulti
          onChange={this.handleChange}
          onBlur={this.handleBlur}
          value={this.props.value}
        />
        {!!this.props.error && this.props.touched && (
          <div style={{ color: "red", marginTop: ".5rem" }}>
            {this.props.error}
          </div>
        )}
      </div>
    );
  }
}

class PublicationSelect extends React.Component {

  handleChange = value => {
    // this is going to call setFieldValue and manually update values.topcis
    this.props.onChange("pmid", value);
  };

  handleBlur = () => {
    // this is going to call setFieldTouched and manually update touched.topcis
    this.props.onBlur("pmid", true);
  };

  render() {
    return (
      <div style={{ margin: "0.25rem 0" }}>
        <label for="pmid">Publication PMID</label>
        <Select
          id="pmid"
          options={publicationOptions}
          isMulti
          onChange={this.handleChange}
          onBlur={this.handleBlur}
          value={this.props.value}
        />
        {!!this.props.error && this.props.touched && (
          <div style={{ color: "red", marginTop: ".5rem" }}>
            {this.props.error}
          </div>
        )}
      </div>
    );
  }
}

const QueryEnhancedForm = formikEnhancer(QueryForm);

//export default QueryEnhancedForm;


function mapDispatchToProps(dispatch) {
  return bindActionCreators({ fetchBygeneSymbol, fetchGeneSummary, fetchCancerDriverList }, dispatch);
}

export default connect(null, mapDispatchToProps)(QueryEnhancedForm);