import React, { useEffect, useRef } from "react";
import * as d3 from "d3";
//import { connect } from "react-redux";
//import { fetchCancerDriverList } from "../../actions";

const Pie = props => {
  
  const ref = useRef(null);

  //const out = d3.extent(props.data,props.arcAccessor)
  //console.log("out", out)

  useEffect(
    () => {  

  const dimensionCount = d3.nest()
     .key(props.arcAccessor)
     .rollup( v => v.length )
     .entries(props.data)
     .sort(function(a, b){ return d3.descending(a.value, b.value); });

  console.log("dimensionCount",dimensionCount);



  // create a new pie
  const createPie = d3
    .pie()
    .value(d => d.value);
    //.sort(function(a, b){ return d3.descending(a.value, b.value); })
    //.sort(null);

  const createArc = d3
    .arc()
    .innerRadius(props.innerRadius * 0.5)
    .outerRadius(props.outerRadius * 0.8)
    .cornerRadius(3)
    .padAngle(0.015)
    

  const createOuterArc = d3
    .arc()
    .innerRadius(props.innerRadius * 1.0)
    .outerRadius(props.outerRadius * 1.1)
    .cornerRadius(3)
    .padAngle(0.015)
    


  const createOuter2Arc = d3
    .arc()
    .innerRadius(props.innerRadius * 1.0)
    .outerRadius(props.outerRadius * 1.3)
    .cornerRadius(3)
    .padAngle(0.015)
    


  const createArcOver = d3
    .arc()
    .innerRadius(props.innerRadius * 0.5)
    .outerRadius(props.outerRadius * 0.8 + 10) 
    .cornerRadius(3)
    .padAngle(0.015)
    


  // color scale 
  //const colors = d3.scaleOrdinal(d3.schemeCategory10);
  const colors = d3.scaleOrdinal()
                   .domain(dimensionCount.map(d => d.value))
                   .range(d3.schemeCategory10);

  const format = d3.format("d");

  //useEffect(
  //  () => {


      //const pie = createPie(props.data);
      const pie = createPie(dimensionCount);

      const svg = d3.select(ref.current)
                    .attr("width",props.width)
                    .attr("height",props.height)
                    .append("g")
                    .attr("transform","translate(" + props.width / 2 + "," + props.height / 2 + ")");

      //svg.append('g').attr('class', 'slices');

      const groupWithData = svg
         .selectAll("g.slices")
         .data(pie) 

                    //group.append('g').attr('class', 'slices');
                    //group.append('g').attr('class', 'innerLabelName');
                    //group.append('g').attr('class', 'outerLabelName');
                    //group.append('g').attr('class', 'lines');
      
      //const groupWithData = group.selectAll("g.slices").data(pie);
      
      groupWithData.exit().remove();


      console.log("group with data",groupWithData);              
      
      const groupWithUpdate = groupWithData
         .enter()
         .append("g")
         .attr("class","slices");

      console.log("group with update",groupWithUpdate);  

      //const path = group.select('.slices')
      //  .datum(props.data).selectAll('path')
      //  .data(pie)

      const path = groupWithUpdate
        .append("path")
        .merge(groupWithData.select(".slices path")); 
 
      const debug = groupWithData.select(".slices path")   

      console.log("path after merge",debug);  

      path
        .attr("class","slices")  
        .attr("d", createArc)
        .attr("data-index", function(d, i){ return i;})
        .attr("fill", (d, i) => colors(i))
        .attr("stroke", "white")
        .style("stroke-width", "1px")
        

      const text = svg.select('.innerLabelName')
        .selectAll('text')
        .data(pie)
        .enter()
        .append("text")
        

      text
        .attr("text-anchor", "middle")
        .attr("alignment-baseline", "middle")
        .attr("transform", d => `translate(${createArc.centroid(d)})`)
        .style("fill", "white")
        .style("font-size", 12)
        .text(d => format(d.value));

      // create outer-text label
      const outerTextLabel = svg.select(".outerLabelName")
        .selectAll("text")
        .data(pie)
        .enter()
        .append("text")
        

      outerTextLabel
        .attr("text-anchor", function(d){
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2
              return (midangle < Math.PI ? 'start' : 'end')
        })
        .attr("transform", function(d, i) {
              const pos = createOuter2Arc.centroid(d);
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2
              pos[0] = props.outerRadius * 0.95 * (midangle < Math.PI ? 1 : -1);

              const percent = (d.endAngle - d.startAngle)/(2*Math.PI)*100
              if( percent<5) {
                  console.log("text label percent",percent)
                  console.log("text label index",i)
                  pos[1] -= (i-3)*12
              }

              console.log("text label pos:", pos);
              return 'translate(' + pos + ')';
        })
        .attr("data-index", function(d, i){ return i;})
        .style("font-size", 12)
        .text(function(d,i){
            //console.log("text log", d) 
            return d.data.key + ' ('+ format(d.value) + ')';
        });
     
      // create outer-text line
      const outerTextLine = svg.select(".lines")
      .selectAll("polyline")
      .data(pie)
      .enter()
      .append("polyline")
        

      outerTextLine
        .attr("stroke", "black")
        .style("fill", "none")
        .attr("stroke-width", 1)
        .attr('points', function(d,i) {
              const posA = createOuterArc.centroid(d) // line insertion in the slice
              const posB = createOuter2Arc.centroid(d) // line break: we use the other arc generator that has been built only for that
              const posC = createOuter2Arc.centroid(d); // Label position = almost the same as posB
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2 // we need the angle to see if the X position will be at the extreme right or extreme left
              posC[0] = props.outerRadius * 0.9 * (midangle < Math.PI ? 1 : -1); // multiply by 1 or -1 to put it on the right or on the left
              
              const percent = (d.endAngle - d.startAngle)/(2*Math.PI)*100
              if( percent<5) {
                  //console.log("text label percent",percent)
                  //console.log("text label index",i)
                  posC[1] -= (i-3)*12
              }
              
              posB[1] = posC[1]
              
              return [posA, posB, posC]
        })


      d3.selectAll('.outerLabelName text, .slices path').call(arcExtend)
        


      function arcExtend(selection) {
          
           selection.on('mouseover', function(d) {
                    console.log("move",this);

                    const currentEl = d3.select(this);
                    let segment, index;

                    if (currentEl.attr("class") === "slices") {
                      segment = currentEl.select("path");
                    } else {
                      index = currentEl.attr("data-index");
                      //const segment = d3.select("#" + pie.cssPrefix + "segment" + index);
                    }


                    d3.select(this).transition()
                    .duration(100)
                    .style("opacity",0.7)
                    .attr("d",createArcOver);
                    console.log("move class",currentEl.attr("class"));    
           });

           selection.on('mouseout', function(d) {
                    d3.select(this).transition()
                    .duration(100)
                    .style("opacity",1)
                    .attr("d",createArc);
           });
      }



    },
    [props.data]
  );

  return (
    <svg 
        width={props.width} 
        height={props.height}
        ref={ref}
    />
  );
};

export default Pie;
