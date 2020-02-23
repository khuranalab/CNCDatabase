import React, { useEffect, useRef, useContext } from "react";
import * as d3 from "d3";
//import { connect } from "react-redux";
//import { fetchCancerDriverList } from "../../actions";
import { PlotContext } from "./Context";
import { elementType } from "prop-types";


const Pie = props => {
  
  const ref = useRef(null);

  const plotContext = useContext(PlotContext); 
  const { groupOptions, selectedFilters, setSelectedFilters } = plotContext;
  //const out = d3.extent(props.data,props.arcAccessor)
  //console.log("out", out)
  
  //https://stackoverflow.com/questions/105034/create-guid-uuid-in-javascript
  //function uuidv4() {
  //  return ([1e7]+-1e3+-4e3+-8e3+-1e11).replace(/[018]/g, c =>
  //    (c ^ crypto.getRandomValues(new Uint8Array(1))[0] & 15 >> c / 4).toString(16)
  //  );
  //}

  function uuidv4() {
    return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
      var r = Math.random() * 16 | 0, v = c == 'x' ? r : (r & 0x3 | 0x8);
      return v.toString(16);
    });
  }

  //useEffect(
  //  () => {  

  // group by entries      
  const dimensionCount = d3.nest()
      .key(props.arcAccessor)
      .rollup( v => v.length )
      .entries(props.data)
      .sort(function(a, b){ return d3.descending(a.value, b.value); });

  // group by unique gene counts    
  const dimensionCount2 = d3.nest()
     .key(props.arcAccessor)
     .rollup( function(v){
                return d3.nest().key(d => d.gene).rollup(k => k.length).entries(v).length
            })
     .entries(props.data)
     .sort(function(a, b){ return d3.descending(a.value, b.value); });

  let i = 0, ln = dimensionCount2.length;
  for (i;i<ln;i++){
         dimensionCount2[i].uniqueId = uuidv4();
  }


  //console.log("dimensionCount2",dimensionCount2);

  i = 0; ln = dimensionCount.length;
  for (i;i<ln;i++){
      dimensionCount[i].uniqueId = uuidv4();
  }

  //console.log("dimensionCount",dimensionCount);

  const filteredData = ( groupOptions === "entries" ) ? dimensionCount : dimensionCount2;
  //console.log("plot group options", groupOptions);

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
    .startAngle(function(d) { return d.startAngle + Math.PI*(1); })
    .endAngle(function(d) { return d.endAngle + Math.PI*(1); });
    

  const createOuterArc = d3
    .arc()
    .innerRadius(props.innerRadius * 1.0)
    .outerRadius(props.outerRadius * 1.1)
    .cornerRadius(3)
    .padAngle(0.015)
    .startAngle(function(d) { return d.startAngle + Math.PI*(1); })
    .endAngle(function(d) { return d.endAngle + Math.PI*(1); });


  const createOuter2Arc = d3
    .arc()
    .innerRadius(props.innerRadius * 1.0)
    .outerRadius(props.outerRadius * 1.3)
    .cornerRadius(3)
    .padAngle(0.015)
    .startAngle(function(d) { return d.startAngle + Math.PI*(1); })
    .endAngle(function(d) { return d.endAngle + Math.PI*(1); });


  const createArcOver = d3
    .arc()
    .innerRadius(props.innerRadius * 0.5)
    .outerRadius(props.outerRadius * 0.8 + 10) 
    .cornerRadius(3)
    .padAngle(0.015)
    .startAngle(function(d) { return d.startAngle + Math.PI*(1); })
    .endAngle(function(d) { return d.endAngle + Math.PI*(1); });


  // color scale 
  //const colors = d3.scaleOrdinal(d3.schemeCategory10);
  const colors = d3.scaleOrdinal()
                   //.domain(dimensionCount.map(d => d.value))
                   .range(d3.schemeCategory10);

  const format = d3.format("d");

  useEffect(
    () => {


      //const pie = createPie(props.data);
      //const pie = createPie(dimensionCount2);
      const pie = createPie(filteredData);

      const svg = d3.select(ref.current)
                    //.attr("width",props.width)
                    //.attr("height",props.height)
                    //.append("g")
                    //.attr("transform","translate(" + props.width / 2 + "," + props.height / 2 + ")");

      const group = svg
                    //.append("g")
                    //.attr("transform","translate(" + props.width * (0.5) + "," + props.height * (0.15) + ")");


      //console.log("group info",ref.current)
            //group.append("g")
            //   .attr("class","slices");
            //group.append("g")
            //   .attr("class","innerLabelName");
            //svg.append("g")
            //   .attr("class","outerLabelName");
            //svg.append("g")
            //   .attr("class","lines");

            //svg.attr("transform","translate(" + props.width / 2 + "," + props.height / 2 + ")");

      //svg.append('g').attr('class', 'slices');

      const slicesWithData = group
         .selectAll("g.slices")
         .data(pie) 

                    //group.append('g').attr('class', 'slices');
                    //group.append('g').attr('class', 'innerLabelName');
                    //group.append('g').attr('class', 'outerLabelName');
                    //group.append('g').attr('class', 'lines');
      
      //const slicesWithData = svg.selectAll("g.slice").data(pie);
      //console.log("group with data before",slicesWithData); 
      slicesWithData.exit().remove();


      //console.log("group with data after",slicesWithData);              
      
      const slicesWithUpdate = slicesWithData
         .enter()
         .append("g")
         .attr("transform","translate(" + props.width * (0.5) + "," + props.height * (0.15) + ")")
         .attr("class","slices");

      //console.log("group with update",groupWithUpdate);  

      //const path = group.select('.slices')
      //  .datum(props.data).selectAll('path')
      //  .data(pie)

      const path = slicesWithUpdate
        .append("path")
        .merge(slicesWithData.select("path.slices")); 
 
      //const debug = groupWithData.select(".slices path")   

      //console.log("path after merge",debug);  

      path
        .attr("class","slices")  
        .attr("d", createArc)
        .attr("id", function(d,i){ return "segment" + d.data.uniqueId + "index" + i;})
        .attr("data-index", function(d, i){ return i;})
        .attr("fill", (d, i) => colors(i))
        .attr("stroke", "white")
        .style("stroke-width", "1px")
      
      const innerLabelText = slicesWithUpdate
        .append("text")
        .merge(slicesWithData.select("text.innerLabelText"))

      innerLabelText
        .attr("class","innerLabelText")
        .attr("data-index", function(d, i){ return i;})
        .attr("text-anchor", "middle")
        .attr("alignment-baseline", "middle")
        .attr("transform", d => `translate(${createArc.centroid(d)})`)
        .style("fill", "white")
        .style("font-size", 12)
        .text(d => format(d.value));  





      const outerLabelText = slicesWithUpdate
        .append("text")
        .merge(slicesWithData.select("text.outerLabelText"))
        

      outerLabelText
        .attr("class","outerLabelText")
        .attr("data-index", function(d, i){ return i;})
        .attr("text-anchor", function(d){
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2
              return (midangle < Math.PI ? 'end' : 'start')
        })
        .attr("transform", function(d,i) {
              const pos = createOuter2Arc.centroid(d);
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2
              pos[0] = props.outerRadius * 1 * (midangle < Math.PI ? -1 : 1);

              const percent = (d.endAngle - d.startAngle)/(2*Math.PI)*100
              if( percent<5) {
                  //console.log("text label percent",percent)
                  //console.log("text label index",i)
                  pos[1] += (i-1)*12
              }

              //console.log("text label pos:", pos);

              return 'translate(' + pos + ')';
        })
        .style("font-size", "0.7rem")
        .text(function(d,i){
            //console.log("text log", d) 
            return d.data.key + ' ('+ format(d.value) + ')';
        });

      // create outer-text line
      const outerLabelLine = slicesWithUpdate
      .append("polyline")
      .merge(slicesWithData.select("polyline"))

      outerLabelLine
        .attr("class","outerLabelLine")
        .attr("stroke", "black")
        .style("fill", "none")
        .attr("stroke-width", 1)
        .attr('points', function(d,i) {
              const posA = createOuterArc.centroid(d) // line insertion in the slice
              const posB = createOuter2Arc.centroid(d) // line break: we use the other arc generator that has been built only for that
              const posC = createOuter2Arc.centroid(d); // Label position = almost the same as posB
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2 // we need the angle to see if the X position will be at the extreme right or extreme left
              posC[0] = props.outerRadius * 0.9 * (midangle < Math.PI ? -1 : 1); // multiply by 1 or -1 to put it on the right or on the left
              
              const percent = (d.endAngle - d.startAngle)/(2*Math.PI)*100
              if( percent<5) {
                  //console.log("text label percent",percent)
                  //console.log("text label index",i)
                  posC[1] += (i-1)*12
              }
              
              posB[1] = posC[1]
              
              return [posA, posB, posC]
        })



        d3.selectAll('text.outerLabelText, text.innerLabelText, path.slices').call(arcExtend)
        


        function arcExtend(selection) {
            
             selection.on('mouseover', function(d) {
                      //console.log("move",this);
  
                      const currentEl = d3.select(this);
                      let segment, index;
  
                      console.log("selection",d)

                      if (currentEl.attr("class") === "slices") {
                        segment = currentEl;
                      } else {
                        index = currentEl.attr("data-index");
                        segment = d3.select("#" + "segment" + d.data.uniqueId + "index" + index);
                      }
  
  
                      //d3.select(this).transition()
                      segment.transition()
                        .duration(100)
                        .style("opacity",0.7)
                        .attr("d",createArcOver);
                        //console.log("move class",currentEl.attr("class"));
                        //console.log("move segment",segment);    
             });
  
             selection.on('mouseout', function(d) {
                      
                      const currentEl = d3.select(this);
                      let segment, index;

                      if (currentEl.attr("class") === "slices") {
                        segment = currentEl;
                      } else {
                        index = currentEl.attr("data-index");
                        segment = d3.select("#" + "segment" + d.data.uniqueId + "index" + index);
                      }
                      //d3.select(this).transition()
                      segment.transition()
                        .duration(100)
                        .style("opacity",1)
                        .attr("d",createArc);
             });
        }


        const labelBBoxes = d3.selectAll('text.outerLabelText')
                                  

        //console.log("text node",labelBBoxes)

        let text_y_coord = []; 

        labelBBoxes.each(function(d,i){
            
           let string = d3.select(this).attr("transform");
           let translate = string.substring(string.indexOf("(")+1, string.indexOf(")")).split(",");
           let x_coord = translate[0];
           let y_coord = translate[1];
           text_y_coord.push(y_coord);
           //console.log("bbox coord",y_coord);
        });

        let max_text_y_coord = Math.max(...text_y_coord);

        //console.log("array",max_text_y_coord);

        let new_y_boundary = max_text_y_coord>0 ? max_text_y_coord : 0; 

        let new_svg_height = new_y_boundary + props.height*0.15 + 30

        //console.log("check box coord",labelBBoxes)        


        svg
         .attr("width","500")
         .attr("height",new_svg_height) 

        
      
    },
    [props.data, groupOptions]
  );

  return (
    /*<svg width = {props.width} height = {props.height}>*/
    <svg 
        width = {props.width} 
        height = {props.height}
        ref={ref}
    >
    </svg>
  );
};

export default Pie;
