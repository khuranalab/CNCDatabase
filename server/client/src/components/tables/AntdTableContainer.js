import React from 'react';
import { Table, Divider, Tag } from 'antd';


//import BootstrapTable from 'react-bootstrap-table-next';
//import paginationFactory from 'react-bootstrap-table2-paginator';


class AntdTableContainer extends React.Component {
     
    



    render (){
      
      const columns = [{
        title: 'GeneName',
        dataIndex: 'geneName',
        key: 'geneName',
        
      }, 
      {
        title: 'EnsemblID',
        dataIndex: 'ensemblId',
        key: 'ensemblId',
      },
      {
        title: 'ElementType',
        dataIndex: 'ncrType',
        key: 'ncrType',
      },
      {
        title: 'Chr',
        dataIndex: 'chr',
        key: 'chr',
      },
      {
        title: 'StartBp',
        dataIndex: 'startbp',
        key: 'startbp',
      },
      {
        title: 'EndBp',
        dataIndex: 'endbp',
        key: 'endbp',
      },
      {
        title: 'Width',
        dataIndex: 'width',
        key: 'width',
      },
      {
        title: 'hotspot',
        dataIndex: 'hotspot',
        key: 'hotspot',
      },
      {
        title: 'trinucleotide.context',
        dataIndex: 'trinucleotide.context',
        key: 'trinucleotide.context',
      },
      {
        title: 'mutationType',
        dataIndex: 'mutationType',
        key: 'mutationType',
      },
      {
        title: 'mutationConsequence',
        dataIndex: 'mutationConsequence',
        key: 'mutationConsequence',
      },
      {
        title: 'numMutaitons',
        dataIndex: 'numMutations',
        key: 'numMutations',
      },
      {
        title: 'numSamples',
        dataIndex: 'numSamples',
        key: 'numSamples',
      },
      {
        title: 'cancerType',
        dataIndex: 'cancerType',
        key: 'cancerType',
      },
      {
        title: 'Adj.pValue',
        dataIndex: 'adj.pValue',
        key: 'adj.pValue',
      },
      {
        title: 'DetectionTechnique',
        dataIndex: 'detectionTechnique',
        key: 'detectionTechnique',
      },
      {
        title: 'ValidationTechnique',
        dataIndex: 'validationTechnique',
        key: 'validationTechnique',
      },
      {
        title: 'TotalNumSamples',
        dataIndex: 'totalNumSamples',
        key: 'totalNumSamples',
      },
      {
        title: 'PublicationTitle',
        dataIndex: 'publicationName',
        key: 'publicationName',
      },
      {
        title: 'Journal',
        dataIndex: 'journal',
        key: 'journal',
      },
      {
        title: 'PMID',
        dataIndex: 'PMID',
        key: 'PMID',

      },
      {
        title: 'dataSources',
        dataIndex: 'dataSources',
        key: 'dataSources',
        
      }
      ];
      
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
        "detectionTechnique": "All mutations within 50 bp of one another were merged into hotspot clusters until no cluster was within 50 bp of another cluster. Clusters with fewer than two mutations were discarded. A p-value was calculated for each cluster using a negative binomial distribution that considered the length, number of mutations, and background mutation rate associated with the cluster. The background mutation rate was calculated as the mean background mutation probability for each sample that had a mutation represented in the cluster. The background mutation probability of each sample was calculated as the total number of mutations divided by the size of the genome. The p-values were adjusted using the Benjamini-Hochberg method.",
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
            "detectionTechnique": "All mutations within 50 bp of one another were merged into hotspot clusters until no cluster was within 50 bp of another cluster. Clusters with fewer than two mutations were discarded. A p-value was calculated for each cluster using a negative binomial distribution that considered the length, number of mutations, and background mutation rate associated with the cluster. The background mutation rate was calculated as the mean background mutation probability for each sample that had a mutation represented in the cluster. The background mutation probability of each sample was calculated as the total number of mutations divided by the size of the genome. The p-values were adjusted using the Benjamini-Hochberg method.",
            "validationTechnique": "none",
            "totalNumSamples": 0,
            "publicationName": "Genome-wide analysis of noncoding regulatory mutations in cancer",
            "journal": "10.1038/ng.3101",
            "PMID": 25261935,
            "dataSources": "CGHub;Broad Genome Data Analysis Center (GDAC) Firehose;Alexandrov et al"
            }
        ];


         return (
                <Table
                dataSource={ data }
                columns={ columns }
                />

         );

    }

}

export default AntdTableContainer;