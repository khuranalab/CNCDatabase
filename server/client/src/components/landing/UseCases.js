import React, { Component } from 'react';
//import Collapsible from 'react-collapsible';
//import { connect } from 'react-redux';
//import { Link } from 'react-router-dom';
//import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
//import { InputGroup, InputGroupAddon, InputGroupText, Input, Button } from 'reactstrap';
//import { UndrawBrowserStats } from 'react-undraw-illustrations';
import Panel from './UseCases/CollapsibleComponent/Panel';

class UseCases extends Component {
  render() {
    return (
      <div className="container-fluid pl-5 pb-5">
         <div className="row">
              <div className="col-lg-10">
                  <div className="pb-3"> 
                    <h3><strong> Non-coding cancer driver examples </strong></h3>
                  </div>
                   
                   <Panel title="TERT promoter">
                     <hr/>
                     <p>Optio asperiores sit et dolores est in qui. Et dolores aut. Sunt non sint optio vero esse reiciendis velit aliquam. Exercitationem minus id corporis. Facere qui nisi asperiores laboriosam ea nihil. Fugiat facilis suscipit voluptatem recusandae sit quia quisquam molestiae.</p>
                     <p>At quae sint temporibus recusandae velit <a href="#0">non officia</a>. A optio suscipit non aliquam aut exercitationem ipsam. Autem rerum aperiam eos. Dolore necessitatibus harum nisi reprehenderit ea quia molestiae sint et. Cupiditate nostrum quia iure hic nulla vel. Quisquam ducimus nulla explicabo labore cum.</p>
                   </Panel>
            
                   <Panel title="ESR1 enhancer">
                     <hr/>
                     <p>Optio asperiores sit et dolores est in qui. Et dolores aut. Sunt non sint optio vero esse reiciendis velit aliquam. Exercitationem minus id corporis. Facere qui nisi asperiores laboriosam ea nihil. Fugiat facilis suscipit voluptatem recusandae sit quia quisquam molestiae.</p>
                     <p>At quae sint temporibus recusandae velit <a href="#0">non officia</a>. A optio suscipit non aliquam aut exercitationem ipsam. Autem rerum aperiam eos. Dolore necessitatibus harum nisi reprehenderit ea quia molestiae sint et. Cupiditate nostrum quia iure hic nulla vel. Quisquam ducimus nulla explicabo labore cum.</p>
                   </Panel>

                   <Panel title="LncRNA NEAT1">
                     <hr/>
                     <p>Optio asperiores sit et dolores est in qui. Et dolores aut. Sunt non sint optio vero esse reiciendis velit aliquam. Exercitationem minus id corporis. Facere qui nisi asperiores laboriosam ea nihil. Fugiat facilis suscipit voluptatem recusandae sit quia quisquam molestiae.</p>
                     <p>At quae sint temporibus recusandae velit <a href="#0">non officia</a>. A optio suscipit non aliquam aut exercitationem ipsam. Autem rerum aperiam eos. Dolore necessitatibus harum nisi reprehenderit ea quia molestiae sint et. Cupiditate nostrum quia iure hic nulla vel. Quisquam ducimus nulla explicabo labore cum.</p>
                   </Panel>

              </div>      

         </div>
         
      </div>
    );
  }
}

export default UseCases;
