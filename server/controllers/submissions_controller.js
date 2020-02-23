const Submission = require('../models/submission');




module.exports = {
    greeting(req, res){
        res.send({ hi: 'submission'});
    },

    create(req, res, next) {
        const submissionProps = req.body;

        Submission.create(submissionProps)
          .then(submission => res.send(submission))
          .catch(next);

    }
}