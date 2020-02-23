import React, { useEffect, useRef } from "react";
import * as d3 from "d3";

const Pie = props => {
  
  const ref = useRef(null);
  
  // create a new pie
  const createPie = d3
    .pie()
    .value(d => d.value)
    .sort(null);

  const createArc = d3
    .arc()
    .innerRadius(props.innerRadius * 0.5)
    .outerRadius(props.outerRadius * 0.8)

  const createOuterArc = d3
    .arc()
    .innerRadius(props.innerRadius * 1.0)
    .outerRadius(props.outerRadius * 1.1);

  const createOuter2Arc = d3
    .arc()
    .innerRadius(props.innerRadius * 1.1)
    .outerRadius(props.outerRadius * 1.2);

  const createArcOver = d3
    .arc()
    .innerRadius(props.innerRadius * 0.5)
    .outerRadius(props.outerRadius * 0.8 + 10);  
  
  const colors = d3.scaleOrdinal(d3.schemeCategory10);
  
  const format = d3.format(".2f");

  useEffect(
    () => {
      const data = createPie(props.data);
      const group = d3.select(ref.current).attr("transform","translate(" + props.width / 2 + "," + props.height / 2 + ")");
      const groupWithData = group.selectAll("g.arc").data(data);

      groupWithData.exit().remove();

      const groupWithUpdate = groupWithData
        .enter()
        .append("g")
        .attr("class", "arc");

      const path = groupWithUpdate
        .append("path")
        .merge(groupWithData.select("path.arc"));

      path
        .attr("class", "arc")
        .attr("d", createArc)
        .attr("fill", (d, i) => colors(i))
        .attr("stroke", "white")
        .style("stroke-width", "1px")
        .on("mouseover", function(d){
            d3.select(this).transition()
            .duration(100)
            .style("opacity",0.7)
            .attr("d",createArcOver);
        })
        .on("mouseout", function(d){
            d3.select(this).transition()
            .duration(100)
            .style("opacity",1)
            .attr("d",createArc);
        })

      const text = groupWithUpdate
        .append("text")
        .merge(groupWithData.select("text"));

      text
        .attr("text-anchor", "middle")
        .attr("alignment-baseline", "middle")
        .attr("transform", d => `translate(${createArc.centroid(d)})`)
        .style("fill", "white")
        .style("font-size", 10)
        .text(d => format(d.value));

      // create outer-text label
      const outerTextLabel = groupWithUpdate
        .append("text")
        .merge(groupWithData.select("text"));

      outerTextLabel
        .attr("text-anchor", function(d){
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2
              return (midangle < Math.PI ? 'start' : 'end')
        })
        .attr("transform", function(d) {
              const pos = createOuter2Arc.centroid(d);
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2
              pos[0] = props.outerRadius * 0.95 * (midangle < Math.PI ? 1 : -1);
              return 'translate(' + pos + ')';
        })
        .style("font-size", 15)
        .text(d => format(d.value));
     
      // create outer-text line
      const outerTextLine = groupWithUpdate
        .append("polyline")
        .merge(groupWithData.select("polyline"));

      outerTextLine
        .attr("stroke", "black")
        .style("fill", "none")
        .attr("stroke-width", 1)
        .attr('points', function(d) {
              const posA = createOuterArc.centroid(d) // line insertion in the slice
              const posB = createOuter2Arc.centroid(d) // line break: we use the other arc generator that has been built only for that
              const posC = createOuter2Arc.centroid(d); // Label position = almost the same as posB
              const midangle = d.startAngle + (d.endAngle - d.startAngle) / 2 // we need the angle to see if the X position will be at the extreme right or extreme left
              posC[0] = props.outerRadius * 0.9 * (midangle < Math.PI ? 1 : -1); // multiply by 1 or -1 to put it on the right or on the left
              return [posA, posB, posC]
        })


    },
    [props.data]
  );

  return (
    <svg width={props.width} height={props.height}>
      <g
        ref={ref}
        transform={`translate(${props.outerRadius} ${props.outerRadius})`}
      />
    </svg>
  );
};

export default Pie;
