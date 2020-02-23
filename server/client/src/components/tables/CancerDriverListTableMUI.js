import React from 'react';
import { connect } from 'react-redux';
import { fetchCancerDriverList } from '../../actions';
import MUIDataTable from 'mui-datatables';

class CancerDriverListTableMUI extends React.Component {

  componentDidMount() {
    //this.props.fetchBygeneSymbol();
  }

  renderDataBygeneSmbol2(dataBygeneSymbol) {
    //console.log(dataBygeneSymbol);

    const dataObj = [{
        dataField: 'id',
        text: 'Product ID',
        sort: true
      },
      {
        dataField: 'name',
        text: 'Product Name',
        sort: true,
      }, {
        dataField: 'price',
        text: 'Product Price',
        sort: true
      }]

      const peopleObj = {
        jim: {
          name: "jim",
          age: 20
        },
        tim: {
          name: "tim",
          age: 22
        }
      }
      //const peopleArray = Object.values(peopleObj)
      //const dataArray = dataObj.map(function(o) {return Object.values(o);})
      //const dataArray =JSON.stringify(dataObj)
      //const dataArray = dataObj.map(function(o) {return Object.values(o);})
      //const dataArrayKey = Object.keys(dataObj[0])
      const object1 = {
        geneName: "TERT",
        ensemblId: "ENSG00000164362",
        ncrType: "promoter",
        chr: "5",
        startbp: 1295227,
        endbp: 1295254,
        width: 28,
        hotspot: true,
        "trinucleotide.context": "unknown",
        mutationType: "unspecified single nucleotide substitution",
        mutationConsequence: "unknown",
        numMutations: 54,
        numSamples: 10,
        cancerType: "BLCA",
        "adj.pValue": 3.07e-153,
        detectionTechnique: "All mutations within 50 bp of one another were merged into hotspot clusters until no cluster was within 50 bp of another cluster. Clusters with fewer than two mutations were discarded. A p-value was calculated for each cluster using a negative binomial distribution that considered the length, number of mutations, and background mutation rate associated with the cluster. The background mutation rate was calculated as the mean background mutation probability for each sample that had a mutation represented in the cluster. The background mutation probability of each sample was calculated as the total number of mutations divided by the size of the genome. The p-values were adjusted using the Benjamini-Hochberg method.",
        validationTechnique: "none",
        totalNumSamples: 0,
        publicationName: "Genome-wide analysis of noncoding regulatory mutations in cancer",
        journal: "10.1038/ng.3101",
        PMID: 25261935,
        dataSources: "CGHub;Broad Genome Data Analysis Center (GDAC) Firehose;Alexandrov et al"
        };
        
        const object2 = dataBygeneSymbol[0]
        //console.log(Object.getOwnPropertyNames(object2));
        //console.log(object2)

      //const columns = ['C1','C2','C3']
      const dataArrayKey = [
        "pmid",
        "cancerType",
        "gene",
        "element",
        "evidencetype",
        "evidencemethod"
      ]
      const dataArray = dataBygeneSymbol.map(function(o) {return Object.values(o);})
      //const dataArrayKey = Object.keys(dataBygeneSymbol[0])

    //const JSONItems = JSON.parse(columns);
    //console.log(dataArrayKey);
    //console.log(object2) 

    const options = {
        filterType: 'dropdown',
        responsive: 'stacked'
      };

    return (
        //<div> test </div>
        <div>
        <MUIDataTable
                title={''}
                data={ dataArray }
                columns={ dataArrayKey }
                options={options}
                />
        </div>

    );
  }


  render() {
    //console.log(this.props.dataBygeneSymbol);

    

    return (
      //<div>"test"</div>
      //<div>{this.props.dataBygeneSymbol.map(this.renderDataBygeneSmbol)}</div>
      <div>{ this.renderDataBygeneSmbol2(this.props.dataBygeneSymbol) }</div>

      //      <div className="container" style={{ marginTop: 50 }}>
      //          <MUIDataTable
      //          title={'Search results'}
      //          data={ data }
      //          columns={ columns }
      //          options={options}
      //          />
      //      </div>


    );
  }
}

const mapStateToProps = ( state ) => {
    return { dataBygeneSymbol: state.dataCancerDriverList };
}
  
export default connect(mapStateToProps, { fetchCancerDriverList })(CancerDriverListTableMUI);

//export default TestTable;