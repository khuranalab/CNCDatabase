import React from 'react';
import { connect } from 'react-redux';
import { fetchBygeneSymbol } from '../../actions';
import { Table, Divider, Tag } from 'antd';


//import BootstrapTable from 'react-bootstrap-table-next';
//import paginationFactory from 'react-bootstrap-table2-paginator';


class AntdTableContainer2 extends React.Component {
     
    componentDidMount() {
        // action
        this.props.fetchBygeneSymbol();
    }

    renderDataBygeneSmbol2(dataBygeneSymbol) {

        //const data = dataBygeneSymbol.map(function(o) {return Object.values(o);})
        const data = dataBygeneSymbol
        //console.log(data)

        const columns = [{
            title: 'PMID',
            dataIndex: 'pmid',
            key: 'pmid',
            width: 120,
            fixed: 'left',
          }, 
          {
            title: 'Cancer Type',
            dataIndex: 'cancertype',
            key: 'cancertype',
            width: 250,
            sorter: (a, b) => a.ensemblId - b.ensemblId,
          },
          {
            title: 'Type',
            dataIndex: 'cancertype2',
            key: 'cancertype2',
            width: 120,
          },
          {
            title: 'GeneName',
            dataIndex: 'gene',
            key: 'gene',
            width: 120,
          },
          {
            title: 'Ensembl ID',
            dataIndex: 'ensemblID',
            key: 'ensemblID',
            width: 120,
          },
          {
            title: 'Mutation Type',
            dataIndex: 'mutationtype',
            key: 'mutationtype',
            width: 120,
          },
          {
            title: 'Element',
            dataIndex: 'element',
            key: 'element',
            width: 120,
          },
          {
            title: 'Mutated Sample Size',
            dataIndex: 'mutatedsamplesize',
            key: 'mutatedsamplesize',
            width: 120,
          },
          {
            title: 'Cohort Size',
            dataIndex: 'cohortsize',
            key: 'cohortsize',
            width: 200,
          },
          {
            title: 'Mutation Rate',
            dataIndex: 'mutationrate',
            key: 'mutationrate',
            width: 120,
          },
          {
            title: 'Evidence Type',
            dataIndex: 'evidencetype',
            key: 'evidencetype',
            width: 200,
          },
          {
            title: 'Evidence Method',
            dataIndex: 'evidencemethod',
            key: 'evidencemethod',
            width: 120,
          }
        ];



        return (
            <Table
            dataSource={ data }
            columns={ columns }
            bordered
            size = "middle"
            scroll={{ x: 1500, y: 500 }}
            />

        );

    }

    render (){
      
      const data = [{
        "geneName": "TERT",
        "ensemblId": "ENSG00000164362",
        "ncrType": "promoter",
        "chr": "5",
        "startbp": 1295227,
        "endbp": 1295254,
        "width": 28,
        "hotspot": true,
        "trinucleotide.context": "unknown",
        "mutationType": "unspecified single nucleotide substitution",
        "mutationConsequence": "unknown",
        "numMutations": 54,
        "numSamples": 10,
        "cancerType": "BLCA",
        "adj.pValue": 3.07e-153,
        "detectionTechnique": "All mutations within 50 bp of one another were merged ",
        "validationTechnique": "none",
        "totalNumSamples": 0,
        "publicationName": "Genome-wide analysis of noncoding regulatory mutations in cancer",
        "journal": "10.1038/ng.3101",
        "PMID": 25261935,
        "dataSources": "CGHub;Broad Genome Data Analysis Center (GDAC) Firehose;Alexandrov et al"
        },
        {
            "geneName": "TERT",
            "ensemblId": "ENSG00000164362",
            "ncrType": "promoter",
            "chr": "5",
            "startbp": 1295227,
            "endbp": 1295254,
            "width": 28,
            "hotspot": true,
            "trinucleotide.context": "unknown",
            "mutationType": "unspecified single nucleotide substitution",
            "mutationConsequence": "unknown",
            "numMutations": 54,
            "numSamples": 10,
            "cancerType": "BLCA",
            "adj.pValue": 3.07e-153,
            "detectionTechnique": "All mutations within 50 bp of one another were merged",
            "validationTechnique": "none",
            "totalNumSamples": 0,
            "publicationName": "Genome-wide analysis of noncoding regulatory mutations in cancer",
            "journal": "10.1038/ng.3101",
            "PMID": 25261935,
            "dataSources": "CGHub;Broad Genome Data Analysis Center (GDAC) Firehose;Alexandrov et al"
            }
        ];


         return (
            <div>{ this.renderDataBygeneSmbol2(this.props.dataBygeneSymbol) }</div>
         );

    }

}

const mapStateToProps = ( state ) => {
    return { dataBygeneSymbol: state.dataBygeneSymbol };
}
  
export default connect(mapStateToProps, { fetchBygeneSymbol })(AntdTableContainer2);
  

//export default AntdTableContainer2;