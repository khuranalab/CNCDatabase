import React from 'react';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import './react_collapsible.css';

const MinusIcon = () => {
    return (
        <span className="panel__header-icon mr-2">
            <FontAwesomeIcon icon="chevron-down" color="grey" size="1x"/>
        </span>
    )    
}

export default MinusIcon;