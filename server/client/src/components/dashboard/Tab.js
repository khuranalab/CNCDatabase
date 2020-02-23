import React, { Component } from 'react';

class Tab extends Component {
    onClick = () => {
        const { label, onClick } = this.props;
        onClick(label);
    }

    render(){
        const { onClick, props: { activeTab, label }} = this;

        let className = 'nav-item nav-link';

        if(activeTab === label) {
            className += ' active';
        } 

        return(
            <a className={className} onClick={onClick}>         
              {label}
            </a>  
        )

    }

}

export default Tab;