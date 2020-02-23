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
            title: 'GeneName',
            dataIndex: 'geneName',
            key: 'geneName',
            width: 120,
            fixed: 'left',
          }, 
          {
            title: 'EnsemblID',
            dataIndex: 'ensemblId',
            key: 'ensemblId',
            width: 250,
            sorter: (a, b) => a.ensemblId - b.ensemblId,
          },
          {
            title: 'ElementType',
            dataIndex: 'ncrType',
            key: 'ncrType',
            width: 120,
          },
          {
            title: 'Chr',
            dataIndex: 'chr',
            key: 'chr',
            width: 120,
          },
          {
            title: 'StartBp',
            dataIndex: 'startbp',
            key: 'startbp',
            width: 120,
          },
          {
            title: 'EndBp',
            dataIndex: 'endbp',
            key: 'endbp',
            width: 120,
          },
          {
            title: 'Width',
            dataIndex: 'width',
            key: 'width',
            width: 120,
          },
          {
            title: 'hotspot',
            dataIndex: 'hotspot',
            key: 'hotspot',
            width: 120,
          },
          {
            title: 'trinucleotide.context',
            dataIndex: 'trinucleotide.context',
            key: 'trinucleotide.context',
            width: 200,
          },
          {
            title: 'mutationType',
            dataIndex: 'mutationType',
            key: 'mutationType',
            width: 120,
          },
          {
            title: 'mutationConsequence',
            dataIndex: 'mutationConsequence',
            key: 'mutationConsequence',
            width: 200,
          },
          {
            title: 'numMutaitons',
            dataIndex: 'numMutations',
            key: 'numMutations',
            width: 120,
          },
          {
            title: 'numSamples',
            dataIndex: 'numSamples',
            key: 'numSamples',
            width: 120,
            sorter: (a, b) => a.numSamples - b.numSamples,
          },
          {
            title: 'cancerType',
            dataIndex: 'cancerType',
            key: 'cancerType',
            width: 120,
          },
          {
            title: 'Adj.pValue',
            dataIndex: 'adj.pValue',
            key: 'adj.pValue',
            width: 120,
          },
          {
            title: 'DetectionTechnique',
            dataIndex: 'detectionTechnique',
            key: 'detectionTechnique',
            width: 700,
          },
          {
            title: 'ValidationTechnique',
            dataIndex: 'validationTechnique',
            key: 'validationTechnique',
            width: 200,
          },
          {
            title: 'TotalNumSamples',
            dataIndex: 'totalNumSamples',
            key: 'totalNumSamples',
            width: 200,
          },
          {
            title: 'PublicationTitle',
            dataIndex: 'publicationName',
            key: 'publicationName',
            width: 300,
          },
          {
            title: 'Journal',
            dataIndex: 'journal',
            key: 'journal',
            width: 120,
          },
    
          {
            title: 'PMID',
            dataIndex: 'PMID',
            key: 'PMID',
            width: 120,
    
          },
          {
            title: 'dataSources',
            dataIndex: 'dataSources',
            key: 'dataSources',
            width: 300,
            
          }
        ];



        return (
            <Table
            dataSource={ data }
            columns={ columns }
            bordered
            size = "middle"
            scroll={{ x: 4500, y: 500 }}
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