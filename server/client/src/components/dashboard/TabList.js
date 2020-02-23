import React, { Component } from 'react';
import Tab from './Tab';

class TabList extends Component{
    

    constructor(props){
       super(props);

       this.state = {
           activeTab: this.props.children[0].props.label, 
       };
    }

    onClickTabItem = (tab) => {
        this.setState({ activeTab: tab});
    }

    render(){
        const { onClickTabItem, props: {children}, state: {activeTab}} = this;

        return(
            <div>
                 <nav>
                    <div className="nav nav-pills">
                        { children.map ((child)=> {
                            const { label } = child.props;
                            
                            return(
                                <Tab
                                    activeTab={activeTab}
                                    key={label}
                                    label={label}
                                    onclick={onClickTabItem} 
                                />
                            );
                        })}
                    </div>         
                 </nav>

                 <div className="tab-content">
                   {children.map((child)=>{
                       if(child.props.label !== activeTab) return undefined;
                       return child.props.children
                   })}
                 </div>

            </div> 
        );

    }


}

export default TabList;