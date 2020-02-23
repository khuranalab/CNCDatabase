import React, { useState, useEffect, useRef } from "react"
import * as d3 from "d3"
import { getTimelineData, getScatterData } from "./utils/dummyData"

//import Timeline from "./Timeline"
import ScatterPlot from "./ScatterPlot"
import PiePlot from "./PiePlot_v6"
//import Histogram from "./Histogram"
// import Timeline from "./completed/Timeline"
// import ScatterPlot from "./completed/ScatterPlot"
// import Histogram from "./completed/Histogram"

import "./css/styles.css"

const parseDate = d3.timeParse("%m/%d/%Y")
const dateAccessor = d => parseDate(d.date)
const temperatureAccessor = d => d.temperature
const humidityAccessor = d => d.humidity

const cancerTypeAccessor = d => d.cancertype
const elementTypeAccessor = d => d.elementtype
const evidenceTypeAccessor = d => d.evidencetype
const evidenceMethodAccessor = d => d.evidencemethod

const getData = () => ({
  timeline: getTimelineData(),
  scatter: getScatterData(),
})





const TestPlot = () => {
  const generateData = (value, length = 8) =>
    d3.range(length).map((item, index) => ({
      date: index,
      value: Math.random() * 100
    }));


  const inputData = [
    {date: 1, value: 20, cancertype:"GBM", elementtype: "lncRNA", evidencetype: "computational prediction", evidencemethod: "ExinAtor"},
    {date: 2, value: 10, cancertype:"GBM", elementtype: "promoter", evidencetype: "gene expression association", evidencemethod: "Wilcoxon rank-sum test"},
    {date: 3, value: 40, cancertype:"GBM", elementtype: "promoter", evidencetype: "computational prediction", evidencemethod: "oncodriveFML"},
    {date: 1, value: 20, cancertype:"GBM", elementtype: "promoter", evidencetype: "computational prediction", evidencemethod: "CNCDriver"},
    {date: 2, value: 10, cancertype:"GBM", elementtype: "promoter", evidencetype: "computational prediction", evidencemethod: "CNCDriver"},
    {date: 3, value: 40, cancertype:"GBM", elementtype: "splice site", evidencetype: "computational prediction", evidencemethod: "oncodriveFML"},
    {date: 1, value: 20, cancertype:"GBM", elementtype: "splice site", evidencetype: "computational prediction", evidencemethod: "SMR"},
    {date: 2, value: 10, cancertype:"GBM", elementtype: "splice site", evidencetype: "computational prediction", evidencemethod: "SMR"},
  ];



  //const [data, setData] = useState(generateData(1));
  const [data, setData] = useState(inputData);


  console.log("plot data",data)

  return (
    
      <div className="row">
        <div className="col-lg-3">
            <h1>pie chart </h1> 

            <PiePlot
                data={data}
                width={400}
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
                width={400}
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
                width={400}
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
                width={400}
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

export default TestPlot


