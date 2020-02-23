import React, { Component} from 'react';
import TabList from './TabList';


class TabContent extends Component {
     render(){
         return(
             <div className="container up">
                <TabList>
                    <div label="Gene" className="tab-content">
                    </div>
                          search by gene
                    <div label="Cancer Type" className="tab-content">
                    </div>
                          search by cancer type
                    <div label="Element type" className="tab-content">
                    </div>        
                          search by element type
                </TabList>
             </div>
         )
     }
}

export default TabContent;