const { Pool } = require('pg');
const keys = require('../config/keys');

console.log(keys.user);

const pool = new Pool({
  user: keys.user,
  host: keys.host,
  database: keys.database,
  password: keys.password,
  port: keys.port
});

module.exports = pool;
