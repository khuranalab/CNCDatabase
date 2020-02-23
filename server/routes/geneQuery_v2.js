const pool = require('../db');

module.exports = app => {
  app.get('/persons', (request, response, next) => {
    pool.query('SELECT * FROM "noncodingmutations"', (err, res) => {
      if (err) return next(err);

      response.json(res);
    });
  });

  app.get('/geneName/:genename', (request, response, next) => {
    const { genename } = request.params;
    pool.query(
      'SELECT * FROM "noncodingmutations" WHERE "geneName" = $1',
      [genename],
      (err, res) => {
        if (err) return next(err);

        response.json(res);
      }
    );
  });

  app.get('/api/geneSymbol', (request, response, next) => {
    //res.json({ geneSymbol: req.query.query });

    //const { genename } = request.query.hgncID;
    //response.json({ name: request.query.hgncID });
    //

    console.log("geneSymbol server side ",request.query)

    pool.query(
      'SELECT pmid, cancertype, gene, element, mutatedsamplesize, cohortsize, evidencetype, evidencemethod FROM "noncoding_cancer_driver" WHERE "gene" = $1',
      //["TERT"],
      //'SELECT * FROM "noncodingmutations" WHERE "geneName" = $1',
      [request.query.geneSymbol],
      (err, res) => {
        if (err) return next(err);

        response.json(res);
      }
    );
  });

//######

app.get('/api/cancer_driver_list', (request, response, next) => {
  //res.json({ geneSymbol: req.query.query });

  //const { genename } = request.query.hgncID;
  //response.json({ name: request.query.hgncID });
  //

  console.log("cancer driver server side ",request.query)

  const params = [];
  //let sql = 'SELECT pmid, cancertype, gene, element, mutatedsamplesize, cohortsize, evidencetype, evidencemethod FROM "noncoding_cancer_driver" WHERE "gene" = $1 AND "cancertype" = $2';
  let sql = 'SELECT pmid, cancertype, gene, cosmic, element, mutatedsamplesize, cohortsize, evidencetype, evidencemethod FROM "noncoding_cancer_driver_Feb_20_2020" where 1 = 1';

  //if( request.query.geneSymbol !== ""){
  //     sql += ' and gene = ?';
  //     params.push(request.query.geneSymbol);
  //}

  if( typeof request.query.geneSymbol !== "undefined"){
    console.log("check geneSymbol",request.query.geneSymbol)
    sql += ' and gene = ANY (?)';
    params.push(request.query.geneSymbol);
  } 
  
  console.log("check elementType array1",request.query.elementType)

  if( typeof request.query.elementType !== "undefined"){
    console.log("check elementType array2",request.query.elementType)
    sql += ' and element = ANY (?)';
    params.push(request.query.elementType);
  } 


  // check whether it is empty array 
  if( typeof request.query.cancerType !== "undefined" ){
    console.log("check array",request.query.cancerType)

    sql += ' and cancertype = ANY (?)';
    params.push(request.query.cancerType);
  }

  if( typeof request.query.evidenceType !== "undefined"){
    sql += ' and evidencetype = ANY (?)';
    params.push(request.query.evidenceType);
  }

  if( typeof request.query.pmid !== "undefined"){
    sql += ' and pmid = ANY (?)';
    params.push(request.query.pmid);
  }

  // replace ? to oridnal parameter by toOrdinal from pg-parameterzie
  // https://github.com/bergur/pg-parameterize 
  let index = 1;
  ordinalString = sql;
  while(ordinalString.indexOf('?') !== -1 ){
    ordinalString = ordinalString.replace('?','$' + index)
    index++
  }


  console.log("sql", sql)
  console.log("params", params)
  console.log("ordinalString",ordinalString)

  pool.query(
    //'SELECT pmid, cancertype, gene, element, mutatedsamplesize, cohortsize, evidencetype, evidencemethod FROM "noncoding_cancer_driver" WHERE "gene" = $1 AND "cancertype" = $2',
    //["TERT"],
    //'SELECT * FROM "noncodingmutations" WHERE "geneName" = $1',
    ordinalString,
    params,
    (err, res) => {
      if (err) return next(err);

      response.json(res);
    }
  );
});





//######

  app.get('/api/gene_summary', (request, response, next) => {
    //res.json({ geneSymbol: req.query.query });

    //const { genename } = request.query.hgncID;
    //response.json({ name: request.query.hgncID });
    //
    pool.query(
      //'SELECT symbol, name, locus_group, locus_type,status, ensembl_gene_id, hgnc_id, entrez_id, alias_symbol, gene_family,refseq_accession,omim_id FROM "gene_summary_test" WHERE "symbol" = $1',
      'SELECT external_gene_name, hgnc_symbol, hgncid, ensemblid, entrezgene, gene_type, description FROM "gene_summary" WHERE "external_gene_name" = $1',
      //["TERT"],
      //'SELECT * FROM "noncodingmutations" WHERE "geneName" = $1',
      [request.query.geneSymbol],
      (err, res) => {
        if (err) return next(err);

        response.json(res);
      }
    );
  });




//########

 app.get('/api/v1.0/tumorType', (request, response, next) => {
    //res.json({ geneSymbol: req.query.query });

    //const { genename } = request.query.hgncID;
    //response.json({ name: request.query.hgncID });
    //
    pool.query(
      //'SELECT * FROM "codingmutations" WHERE "tumorType" = $1',
      'SELECT * FROM "noncodingdrivers" WHERE "cancerType" = $1',
      [request.query.tumorType],
      (err, res) => {
        if (err) return next(err);

        response.json(res.rows);
      }
    );
  });



//#########


  app.get('/chr/:chr', (request, response, next) => {
    const { chr } = request.params;
    pool.query(
      'SELECT * FROM "noncodingmutations" WHERE chr = $1',
      [chr],
      (err, res) => {
        if (err) return next(err);

        response.json(res);
      }
    );
  });

  app.get('/', (req, res) => {
    res.send({ hi: 'there' });
  });
};
