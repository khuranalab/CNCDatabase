import React from 'react';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import './react_collapsible.css';

const PlusIcon = () => {
    return (
        <span className="panel__header-icon mr-2">
            <FontAwesomeIcon icon="chevron-right" color="grey" size="1x"/>
        </span>
    )
}

export default PlusIcon;