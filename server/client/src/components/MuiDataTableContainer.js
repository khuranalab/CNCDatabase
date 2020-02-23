import React, { Component } from 'react';
import MUIDataTable from 'mui-datatables';

class MuiDataTableContainer extends Component {

    //https://appdividend.com/2018/05/04/react-bootstrap-table2-example/

    state = {
       dataSources: []; 
       columns:[{

       }]

    }

    
    render(){

        const options = {
            filterType: 'dropdown',
            responsive: 'stacked'
          };

        return(
            <div className="container" style={{ marginTop: 50 }}>
                <MUIDataTable
                title={'Search results'}
                data={ this.state.data }
                columns={ this.state.columns }
                options={options}
                />
            </div>
        );

    }


}

export default MuiDataTableContainer;