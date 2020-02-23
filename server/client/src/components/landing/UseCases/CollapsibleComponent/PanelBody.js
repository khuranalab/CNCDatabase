import React from 'react';
import './react_collapsible.css';

const PanelBody = props => {
    return (
        <div className="panel__body" aria-hidden={props.isExpanded}>
            {props.children}
        </div>
    )
}

export default PanelBody;