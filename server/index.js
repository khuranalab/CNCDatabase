const express = require('express');

//const mongoose = require('mongoose');
//const cookieSession = require('cookie-session');
//const passport = require('passport');
//const keys = require('./config/keys');


const multer = require('multer');
const cors = require('cors');


//const upload = require('./upload');



// cross-origin requests


//require('./models/User');
//require('./services/passport');

const bodyParser = require('body-parser');

//const pool = require('./db');
const routes = require('./routes/geneQuery_v2');

//const authRoutes = require('./routes/authRoutes');

//mongoose.Promose = global.Promise;

//if (process.env.NODE_ENV !== 'test') {
//  mongoose.connect(keys.mongoURI);
//}

const app = express();

// file upload section
//const corsOptions ={
//    origin: '*',
//    optionsSuccessStatus: 200,
//}

//app.use(cors(corsOptions))

app.use(cors())

const storage = multer.diskStorage({
    destination: function (req, file, cb) {
        cb(null, 'public')
    },
      filename: function (req, file, cb) {
        cb(null, Date.now() + '-' +file.originalname )
    }
})

const upload = multer({ storage: storage }).single('file')

app.post('/api/submissions',function(req, res){
   
  upload(req, res, function (err) {
    console.log("Request ---", req.body);
    console.log("Request file ---", req.file);//Here you get file.

    if (err instanceof multer.MulterError) {
        return res.status(500).json(err)
    } else if (err) {
        return res.status(500).json(err)
    }
    
    return res.status(200).send(req.file)

  })



});




// user login section

//app.use(
//  cookieSession({
//     maxAge: 30 * 24 * 60 * 60 * 1000,
//     keys: [keys.cookieKey]
//  })
//);

//app.use(passport.initialize());
//app.use(passport.session());

app.use(bodyParser.urlencoded({ extended: false }));
app.use(bodyParser.json());

routes(app);

//require('./routes/authRoutes')(app);
//require('./routes/studyQuery')(app);
//require('./routes/submissionRoutes')(app);



//app.use((err, req, res, next) => {
//  res.status(404).send({ error: err.message });
//});





//app.get('/', (req, res) => res.sendFile('auth.html', { root : __dirname}));
/*
app.get('/persons', (request, response, next) => {
  pool.query('SELECT * FROM "nonCodingMutations"', (err, res) => {
    //if (err) return next(err);

    response.json(res);
  });
});

app.get('/', (req, res) => {
  res.send({ hi: 'there' });
});
*/

app.get('*', (req, res) => {
  res.status(404).send({ message: '404 error' });
  //console.log("bug")
  //res.status(404).render('index');
  //res.sendFile(path.join(__dirname, './index.html'));
});

//app.use((req, res, next) => {
//   res.status(404).send("404 Not Found");
//  var err = new Error('404 Not Found');
//  err.status = 404;
//  next(err);
//});


//

if (process.env.NODE_ENV === 'production'){
    
    // Express will serve up production assets
    // like main.js and main.css
    app.use(express.static('client/build'));    

    // Express will serve up the index.html file
    // if it doesn't recognize the route    
    const path = require('path');
    app.get('*', (req,res) => {
         res.sendFile(path.resolve(__dirname, 'client', 'build', 'index.html'));
    }); 

}





const PORT = process.env.PORT || 5200;
//PORT = 5000;
//app.listen(PORT);

app.listen(PORT , () => console.log('App listening on port ' + PORT));
