import React, { Component } from 'react';
import { connect } from 'react-redux';
import { fetchBygeneSymbol } from '../../actions';
import MUIDataTable from 'mui-datatables';

class TestTable2 extends Component {

   





  componentDidMount() {
    this.props.fetchBygeneSymbol();
  }

  renderDataBygeneSmbol(dataBygeneSymbol) {

    return (
      <div className="card card-block">
        <h4 className="card-title">{dataBygeneSymbol.geneName}</h4>
        <p className="card-text">Chees Factory</p>
        <a className="btn btn-primary">email</a>
      </div>
    );
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
        "geneName",
        "ensemblId",
        "ncrType",
        "chr",
        "startbp",
        "endbp",
        "width",
        "hotspot",
        "trinucleotide.context",
        "mutationType",
        "mutationConsequence",
        "numMutations",
        "numSamples",
        "cancerType",
        "adj.pValue",
        "detectionTechnique",
        "validationTechnique",
        "totalNumSamples",
        "publicationName",
        "journal",
        "PMID",
        "dataSources"
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
                title={'Search results'}
                data={ dataArray }
                columns={ dataArrayKey }
                options={options}
                />
        </div>

    );
  }


  render() {
    //console.log(this.props.dataBygeneSymbol);

    const options = {
        filterType: 'dropdown',
        responsive: 'stacked'
      };

      const columns = ['Name', 'Title', 'Location', 'Age', 'Salary'];

      const data = [
        ['Gabby George', 'Business Analyst', 'Minneapolis', 30, '$100,000'],
        ['Aiden Lloyd', 'Business Consultant', 'Dallas', 55, '$200,000'],
        ['Jaden Collins', 'Attorney', 'Santa Ana', 27, '$500,000'],
        ['Franky Rees', 'Business Analyst', 'St. Petersburg', 22, '$50,000'],
        ['Aaren Rose', 'Business Consultant', 'Toledo', 28, '$75,000'],
        ['Blake Duncan','Business Management Analyst','San Diego',65,'$94,000'],
        ['Frankie Parry', 'Agency Legal Counsel', 'Jacksonville', 71, '$210,000'],
        ['Lane Wilson', 'Commercial Specialist', 'Omaha', 19, '$65,000'],
        ['Robin Duncan', 'Business Analyst', 'Los Angeles', 20, '$77,000'],
        ['Mel Brooks', 'Business Consultant', 'Oklahoma City', 37, '$135,000'],
        ['Harper White', 'Attorney', 'Pittsburgh', 52, '$420,000'],
        ['Kris Humphrey', 'Agency Legal Counsel', 'Laredo', 30, '$150,000'],
        ['Frankie Long', 'Industrial Analyst', 'Austin', 31, '$170,000'],
        ['Brynn Robbins', 'Business Analyst', 'Norfolk', 22, '$90,000'],
        ['Justice Mann', 'Business Consultant', 'Chicago', 24, '$133,000'],
        ['Addison Navarro','Business Management Analyst','New York',50,'$295,000'],
        ['Jesse Welch', 'Agency Legal Counsel', 'Seattle', 28, '$200,000'],
        ['Eli Mejia', 'Commercial Specialist', 'Long Beach', 65, '$400,000'],
        ['Gene Leblanc', 'Industrial Analyst', 'Hartford', 34, '$110,000'],
        ['Danny Leon', 'Computer Scientist', 'Newark', 60, '$220,000'],
        ['Lane Lee', 'Corporate Counselor', 'Cincinnati', 52, '$180,000'],
        ['Jesse Hall', 'Business Analyst', 'Baltimore', 44, '$99,000'],
        ['Danni Hudson', 'Agency Legal Counsel', 'Tampa', 37, '$90,000'],
        ['Terry Macdonald', 'Commercial Specialist', 'Miami', 39, '$140,000'],
        ['Justice Mccarthy', 'Attorney', 'Tucson', 26, '$330,000'],
        ['Silver Carey', 'Computer Scientist', 'Memphis', 47, '$250,000'],
        ['Franky Miles', 'Industrial Analyst', 'Buffalo', 49, '$190,000'],
        ['Glen Nixon', 'Corporate Counselor', 'Arlington', 44, '$80,000'],
        ['Gabby Strickland','Business Process Consultant','Scottsdale',26,'$45,000'],
        ['Mason Ray', 'Computer Scientist', 'San Francisco', 39, '$142,000']
      ];

      //console.log(data);

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

function mapStateToProps({ dataBygeneSymbol }) {
  return { dataBygeneSymbol };
}

export default connect(mapStateToProps, { fetchBygeneSymbol })(TestTable2);

//export default TestTable;