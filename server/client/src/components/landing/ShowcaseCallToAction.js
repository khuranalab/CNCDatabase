import React, { Component } from 'react';
import Collapsible from 'react-collapsible';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
//import { UndrawBrowserStats } from 'react-undraw-illustrations';

class ShowcaseCallToAction extends Component {
  render() {
    return (
      <div className="py-5">
         <div className="container col-lg-10">
            <div className="row">
              <div className="col-lg-6 text-left">
                   
                   <h3><strong> Examples </strong></h3>
                    <Collapsible trigger="TERT promoter">
                      <p>This is the collapsible content. It can be any element or React component you like.</p>
                      <p>It can even be another Collapsible component. Check out the next section!</p>
                    </Collapsible>

                    <Collapsible trigger="ESR1 enhancer">
                      <p>This is the collapsible content. It can be any element or React component you like.</p>
                      <p>It can even be another Collapsible component. Check out the next section!</p>
                    </Collapsible>

                    <Collapsible trigger="non-coding RNA">
                      <p>This is the collapsible content. It can be any element or React component you like.</p>
                      <p>It can even be another Collapsible component. Check out the next section!</p>
                    </Collapsible>

              </div>      

            </div>
         </div>
      </div>
    );
  }
}

export default ShowcaseCallToAction;
