const mongoose = require('mongoose');
const { Schema } = mongoose;

const SubmissionSchema = new Schema({
   
    name: {
        type: String,
        required: false
    },
    
    email: {
        type: String,
        required: false
    },

    affiliation: {
        type: String,
        required: false
    },

    referenceUrl: {
        type: String,
        required: false
    },

    comments: {
        type: String,
        required: false
    },

    file: {
        type: Buffer,
        contentType: String,
        required: false 
    }

});

const Submission = mongoose.model('submission', SubmissionSchema);

module.exports = Submission;