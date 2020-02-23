import React from 'react';
import './react_collapsible.css';

const PanelGroup = props => {
    return (
        <div className="panel-group" role="group">
            {props.children}
        </div>
    )
}

export default PanelGroup;