//import 'materialize-css/dist/css/materialize.css';
import 'bootstrap/dist/css/bootstrap.min.css';

//import $ from 'jquery';
//import Popper from 'popper.js';
//import 'bootstrap/dist/js/bootstrap.bundle.min';

//import 'react-bootstrap-table-next/dist/react-bootstrap-table2.min.css';
import 'react-tabulator/lib/styles.css';
//import 'react-tabulator/lib/css/tabulator.min.css'; // standard theme
import 'react-tabulator/lib/css/tabulator_simple.min.css'; // simple theme
//import 'react-tabulator/lib/css/tabulator_modern.min.css'; // modern theme
//import 'react-tabulator/lib/css/bootstrap/tabulator_bootstrap4.min.css'; // bootstrap4 theme

//import 'antd/dist/antd.css'
import './index.css'
import './custom_theme.css'

import React from 'react';
import ReactDOM from 'react-dom';
import { Provider } from 'react-redux';
import { createStore, applyMiddleware, compose } from 'redux';
import reduxThunk from 'redux-thunk';

//import InfoTable from './InfoTable';
//import App from './components/App';
//import registerServiceWorker from './registerServiceWorker';

//ReactDOM.render(<InfoTable />, document.getElementById('root'));
//registerServiceWorker();

import App from './components/App';
import reducers from './reducers';

//const store = createStore(reducers, {}, applyMiddleware(reduxThunk));
const composeEnhancers = window.__REDUX_DEVTOOLS_EXTENSION_COMPOSE__ || compose;
const store = createStore(reducers, composeEnhancers(applyMiddleware(reduxThunk)) );
//const store = createStore(reducers, applyMiddleware(reduxThunk));

//const createStoreWithMiddleware = applyMiddleware(reduxThunk)(createStore);

ReactDOM.render(
  //<Provider store={store}>
  //<Provider store={createStoreWithMiddleware(reducers)}>
  <Provider store={store}>
    <App />
  </Provider>,
  document.querySelector('#root')
);
