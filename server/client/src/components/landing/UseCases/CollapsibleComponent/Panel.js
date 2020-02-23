import React, { Component } from 'react';
import './react_collapsible.css';
import PanelHeader from './PanelHeader';
import PanelBody from './PanelBody';

// code example: https://codepen.io/levimcg/pen/ZvEPgm/

class Panel extends Component {
    constructor(props) {
        super(props)
        
        this.state = {
            isExpanded: props.openDefault
        }
        
        this.handleToggle = this.handleToggle.bind(this)
        this.handleKeyDown = this.handleKeyDown.bind(this)
    }
    
    handleToggle() {
        this.setState({
            isExpanded: !this.state.isExpanded 
        })
    }
    
    handleKeyDown(event) {
        if(event.keyCode == 40) {
            event.preventDefault();
            this.setState({
                isExpanded: true
            })
        }
        
        if(event.keyCode == 38) {
            event.preventDefault();
            this.setState({
                isExpanded: false
            })
        }
    }
    
    render() {
        return (
            <div className="panel">
                <PanelHeader
                    handleToggle={this.handleToggle}
                    handleKeyDown={this.handleKeyDown}
                    isExpanded={this.state.isExpanded}
                >
                    {this.props.title}
                </PanelHeader>
                <PanelBody isExpanded={!this.state.isExpanded}>
                    {this.props.children}
                </PanelBody>
            </div>
        )
    }
}

export default Panel;



