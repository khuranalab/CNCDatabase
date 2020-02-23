const { Pool } = require('pg');

const pool = new Pool({
  user: 'devUser',
  host: 'localhost',
  database: 'dev_db',
  password: 'testDev',
  port: 5432
});

pool.query('SELECT * FROM persons', (err, res) => {
  if (err) return console.log(err);

  console.log(res);
});
