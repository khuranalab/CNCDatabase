import React, { useState, useEffect, useRef } from "react"
import * as d3 from "d3"
//import { getTimelineData, getScatterData } from "./utils/dummyData"
import { connect } from "react-redux";
import { fetchCancerDriverList } from "../../actions";

//import Timeline from "./Timeline"
//import ScatterPlot from "./ScatterPlot"
import PiePlot from "./PiePlot_v8"
//import Histogram from "./Histogram"
// import Timeline from "./completed/Timeline"
// import ScatterPlot from "./completed/ScatterPlot"
// import Histogram from "./completed/Histogram"

//import "./css/styles.css"

//const parseDate = d3.timeParse("%m/%d/%Y")
//const dateAccessor = d => parseDate(d.date)
//const temperatureAccessor = d => d.temperature
//const humidityAccessor = d => d.humidity

const cancerTypeAccessor = d => d.cancertype
const elementTypeAccessor = d => d.element
const evidenceTypeAccessor = d => d.evidencetype
const evidenceMethodAccessor = d => d.evidencemethod

//const getData = () => ({
//  timeline: getTimelineData(),
//  scatter: getScatterData(),
//})





const TestPlot = ({dataCancerDriverList, fetchCancerDriverList}) => {
  //const generateData = (value, length = 8) =>
  // d3.range(length).map((item, index) => ({
  //   date: index,
  //    value: Math.random() * 100
  //  }));

  const data = dataCancerDriverList

  console.log("inputData",data)

  const inputData2 = [
    {date: 1, value: 20, cancertype:"GBM", element: "lncRNA", evidencetype: "computational prediction", evidencemethod: "ExinAtor"},
    {date: 2, value: 10, cancertype:"GBM", element: "promoter", evidencetype: "gene expression association", evidencemethod: "Wilcoxon rank-sum test"},
    {date: 3, value: 40, cancertype:"GBM", element: "promoter", evidencetype: "computational prediction", evidencemethod: "oncodriveFML"},
    {date: 1, value: 20, cancertype:"GBM", element: "promoter", evidencetype: "computational prediction", evidencemethod: "CNCDriver"},
    {date: 2, value: 10, cancertype:"GBM", element: "promoter", evidencetype: "computational prediction", evidencemethod: "CNCDriver"},
    {date: 3, value: 40, cancertype:"GBM", element: "splice site", evidencetype: "computational prediction", evidencemethod: "oncodriveFML"},
    {date: 1, value: 20, cancertype:"GBM", element: "splice site", evidencetype: "computational prediction", evidencemethod: "SMR"},
    {date: 2, value: 10, cancertype:"GBM", element: "splice site", evidencetype: "computational prediction", evidencemethod: "SMR"},
  ];



  //const [data, setData] = useState(generateData(1));
  //const [data, setData] = useState(inputData2);


  

  console.log("plot data",data)

  return (
    
      <div className="row">
        <div className="col-lg-3">
            <h1>pie chart </h1> 

            <PiePlot
                data={data}
                width={500}
                height={400}
                innerRadius={50}
                outerRadius={80}
                arcAccessor={cancerTypeAccessor}
                title="Element Type"
            />

        </div>

        <div className="col-lg-3">
            <h1>pie chart </h1> 

            <PiePlot
                data={data}
                width={500}
                height={400}
                innerRadius={50}
                outerRadius={80}
                arcAccessor={elementTypeAccessor}
                title="Element Type"
            />

        </div>

        <div className="col-lg-3">
            <h1>pie chart </h1> 

            <PiePlot
                data={data}
                width={500}
                height={400}
                innerRadius={50}
                outerRadius={80}
                arcAccessor={evidenceTypeAccessor}
                title="Element Type"
            />

        </div>

        <div className="col-lg-3">
            <h1>pie chart </h1> 

            <PiePlot
                data={data}
                width={500}
                height={400}
                innerRadius={50}
                outerRadius={80}
                arcAccessor={evidenceMethodAccessor}
                title="Element Type"
            />

        </div>


      </div>

    
  )
}

//export default TestPlot


const mapStateToProps = ( state ) => {
    return { dataCancerDriverList: state.dataCancerDriverList };
}
export default connect( mapStateToProps, { fetchCancerDriverList })(TestPlot);