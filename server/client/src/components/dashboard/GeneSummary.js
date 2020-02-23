import React from 'react';
//import 'react-tabulator/lib/styles.css';
//import { ReactTabulator} from 'react-tabulator';
//import { Link } from 'react-router-dom';
import { connect } from 'react-redux';
import { fetchGeneSummary } from '../../actions';
import { fetchWeather } from '../../actions';
//import _ from 'lodash';

class GeneSummary extends React.Component {
     
    componentDidMount() {
        // action
        this.props.fetchGeneSummary();
        this.props.fetchWeather();
    }

    renderGeneSummary(cityData) {

        //const data = geneSummaryData

        //const c1 = Object.values(data)
        //const data2 = [{ name: 'Hot', name2: 'Hot2'}]
        //console.log("cityData", data2[0].name2)
        console.log("cityData", cityData.map(a => a.symbol))
        //console.log("gene summary datagene symbol", )

        return(
             

             <tr>
                 <li><b>Symbol:</b> { cityData.map(a => a.external_gene_name) } </li>
                 <li><b>Approved name:</b> { cityData.map(a => a.hgnc_symbol) } </li>
                 <li><b>HGNC ID:</b>{ cityData.map(a => a.hgncid) } </li>
                 <li><b>Ensembl ID:</b> { cityData.map(a => a.ensemblid) } </li>
                 <li><b>NCBI Gene ID:</b> { cityData.map(a => a.entrezgene) } </li>
                 <li><b>Locus type:</b> { cityData.map(a => a.gene_type) } </li>
                 <li><b>Description:</b> { cityData.map(a => a.description) } </li>
             </tr>
        );
    }    

    render (){
        console.log("test2",this.props.geneSummaryData)

        return(
            //<div>{ this.props.geneSummaryData.symbol } </div>
            <div>{ this.renderGeneSummary(this.props.geneSummaryData) }</div>
            //<div>{ this.renderGeneSummary(this.props.weather) }</div>
            //<div>{ this.renderGeneSummary(this.props.geneSummaryData) }</div>
        );

    }

}

//export default ReactBootstrapTable;

//
const mapStateToProps = ( state ) => {
    return {
      geneSummaryData: state.geneSummaryData,  
      weather: state.weather  
    };
}
  
export default connect(mapStateToProps, { fetchGeneSummary, fetchWeather })(GeneSummary);

//export default ReactTabulatorTable;