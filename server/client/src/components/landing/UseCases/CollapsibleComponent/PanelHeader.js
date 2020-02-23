import React from 'react';
import MinusIcon from './MinusIcon';
import PlusIcon from './PlusIcon';
import './react_collapsible.css';

const PanelHeader = props => {
    return (
        <button
            className={props.isExpanded ? "panel__header active" : "panel__header"}
            onClick={props.handleToggle}
            onKeyDown={props.handleKeyDown}
            aria-expanded={props.isExpanded}
        >
               
            {props.isExpanded ? <MinusIcon/> : <PlusIcon/>}
            {props.children}
        </button>
    )
}

export default PanelHeader;