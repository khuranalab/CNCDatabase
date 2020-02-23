const SubmissionController = require('../controllers/submissions_controller');




module.exports = app => {

    //app.get('/api/submission', SubmissionController.greeting);
    
    app.post('/api/submissions',SubmissionController.create)

};  