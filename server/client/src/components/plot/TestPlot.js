import React, { useState, useEffect, useRef } from "react"
import * as d3 from "d3"
import { getTimelineData, getScatterData } from "./utils/dummyData"

//import Timeline from "./Timeline"
import ScatterPlot from "./ScatterPlot"
//import Histogram from "./Histogram"
// import Timeline from "./completed/Timeline"
// import ScatterPlot from "./completed/ScatterPlot"
// import Histogram from "./completed/Histogram"

import "./css/styles.css"

const parseDate = d3.timeParse("%m/%d/%Y")
const dateAccessor = d => parseDate(d.date)
const temperatureAccessor = d => d.temperature
const humidityAccessor = d => d.humidity

const getData = () => ({
  timeline: getTimelineData(),
  scatter: getScatterData(),
})

const TestPlot = () => {
  const [data, setData] = React.useState(getData())

  return (
    <div className="App">
      <h1>
        Weather Dashboard
      </h1>
      <div className="App__charts">
        
        <ScatterPlot
          data={data.scatter}
          xAccessor={humidityAccessor}
          yAccessor={temperatureAccessor}
          xLabel="Humidity"
          yLabel="Temperature"
        />
        
      </div>
    </div>
  )
}

export default TestPlot


