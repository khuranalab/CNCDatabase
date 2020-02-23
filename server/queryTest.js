const express = require('express');
const promise = require('bluebird');

const initOptions = {
  promiseLib: promise
};

const databaseConfig = {
  host: 'localhost',
  port: 5432,
  database: 'dev_db',
  user: 'devUser',
  password: 'testDev'
};

const pgp = require('pg-promise')(initOptions);
const db = pgp(databaseConfig);
//const db = pgp('postgresql://devUser:testDev@localhost:5432/dev_db');



const app = express();
const port = 5100;

db.connect()
    .then(obj => {
        obj.done(); // success, release the connection;
    })
    .catch(error => {
        console.log('ERROR:', error.message || error);
    });


app.listen(port);
