import React from 'react';
//import 'react-tabulator/lib/styles.css';
import { ReactTabulator} from 'react-tabulator';
import { Link } from 'react-router-dom';
import { connect } from 'react-redux';
import { fetchBygeneSymbol } from '../../actions';


class ReactTabulatorTable extends React.Component {
     
    //componentDidMount() {
        // action
        //this.props.fetchBygeneSymbol();
    //}

    renderDataBygeneSmbol2(dataBygeneSymbol) {

        const data = dataBygeneSymbol

        

        return (
            <div className="container-fluid" style={{ marginTop: 20 }}>
               
            </div>

        );

    }    

    render (){
      
      const data = [
        {id:1, name:"Oli Bob", age:"12", col:"red", dob:""},
        {id:2, name:"Mary May", age:"1", col:"blue", dob:"14/05/1982"},
        {id:3, name:"Christine Lobowski", age:"42", col:"green", dob:"22/05/1982"},
        {id:4, name:"Brendon Philips", age:"125", col:"orange", dob:"01/08/1980"},
        {id:5, name:"Margret Marmajuke", age:"16", col:"yellow", dob:"31/01/1999"},

      ];

      const pmidFormatter = function(cell, formatterParams){
        //cell - the cell component
        //formatterParams - parameters set for the column

        const data = cell.getData();
        console.log("check",data);
        const href='https://mapquest.com/latlng/' + data.col; 
        return "<a target='_blank' href='https://www.ncbi.nlm.nih.gov/pubmed/" + data.col + "'>" + cell.getValue() + "</a>";
        //return "<Link to='https://mapquest.com/latlng/" + data.col + "'>" + cell.getValue() + "</Link>";
}

      const columns = [
        { title: "Name", field: "name", width: 150 },
        { title: "Age", field: "age", align: "left", formatter: "progress" },
        { title: "Favourite Color", field: "col", formatter: pmidFormatter },
        { title: "Date Of Birth", field: "dob", align: "center" },
        { title: "Rating", field: "rating", align: "center", formatter: "star" },
        { title: "Passed?", field: "passed", align: "center", formatter: "buttonCross" }
      ];

      const options = {
        pagination: "local",
        paginationSize: 3,
      };

      return(
          <div className="container-fluid" style={{ marginTop: 20 }}>
                  <ReactTabulator
                      data={data}
                      columns={columns}
                      tooltips={true}
                      layout={"fitData"}
                      options={options}
                  />
          </div>
      )

        //return(
        //    <div>{ this.renderDataBygeneSmbol2(this.props.dataBygeneSymbol) }</div>
        //);


    }

}

//export default ReactBootstrapTable;

//
//const mapStateToProps = ( state ) => {
//    return { dataBygeneSymbol: state.dataBygeneSymbol };
//}
  
//export default connect(mapStateToProps, { fetchBygeneSymbol })(ReactBootstrapTable);

export default ReactTabulatorTable;