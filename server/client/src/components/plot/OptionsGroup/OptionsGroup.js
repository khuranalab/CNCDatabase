import React, {useState, useCallback, useContext} from 'react';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { fetchGeneSummary } from '../../../actions';
import { PlotContext } from "../Context";



const OptionsGroup = () => {
    
    const plotContext = useContext(PlotContext);

    const { groupOptions, setGroupOptions } = plotContext;

    const [showLegend, setShowLegend] = useState(false);
    const [option1, setOption1] = useState({name: "genes", isActive: true})
    const [option2, setOption2] = useState({name: "entries", isActive: false})

    const showLegendHandler = useCallback( () => {
        setShowLegend( prevShowLegend => !prevShowLegend );
    },[]);

    const option1Handler = useCallback( () => {
        setOption1( prevOption1 => {
            const newStatus = !prevOption1.isActive;
            //console.log("option1 new:",newStatus)
            return { ...prevOption1, isActive: newStatus }

        });
        setOption2( prevOption2 => {
            const newStatus = !prevOption2.isActive;
            //console.log("option1 new:",newStatus)
            return { ...prevOption2, isActive: newStatus }

        });

        setGroupOptions(option1.name);
        //console.log("option1:",option1.isActive)
    },[]);

    const option2Handler = useCallback( () => {
        setOption2( prevOption2 => {
            const newStatus = !prevOption2.isActive;
            //console.log("option2 new:",newStatus)
            return { ...prevOption2, isActive: newStatus }

        });
        setOption1( prevOption1 => {
            const newStatus = !prevOption1.isActive;
            //console.log("option1 new:",newStatus)
            return { ...prevOption1, isActive: newStatus }

        });

        setGroupOptions(option2.name);
        //console.log("option2:",option2.isActive)
        
    },[]);

    console.log("groupOptions status:",groupOptions);

    return(
         <div>
           <div className="d-flex flex-row"> 
            <span className="mr-2 ">Display numbers by:
              <button className="iconButton" onClick={showLegendHandler}>
               <FontAwesomeIcon icon="info-circle" color="grey" size="1x" onClick={()=>{}}/> 
              </button>
            </span> 
            <button type="button" 
                    className={"btn " + (option1.isActive ? "btn-dark": "btn-outline-dark") + " btn-sm mr-1"}
                    onClick={option1Handler}
            >
                {option1.name}
            </button>
            <button type="button" 
                    className={"btn " + (option2.isActive ? "btn-dark": "btn-outline-dark") + " btn-sm mr-1"}
                    onClick={option2Handler}
            >
                {option2.name}
            </button>
           </div>
           
           { showLegend && <div className="d-flex flex-column mt-3">  
            <p><strong>Display numbers by genes:</strong> Numbers are grouped by numbers of unique genes per category in pie chart<br/>
               <strong>Display numbers by entries:</strong> Numbers are grouped by numbers of entries per category in pie chart</p>       
           </div>}
         </div>
    );
   

};

export default OptionsGroup