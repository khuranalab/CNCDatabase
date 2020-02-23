import axios from 'axios';
import { toast } from 'react-toastify';

import { FETCH_BYGENESYMBOL, 
         FETCH_ALL_STUDIES, FETCH_GENE_SUMMARY, 
         FETCH_WEATHER, FETCH_CANCER_DRIVER_LIST, CREATE_SUBMISSION } from './types';

const API_KEY = '2ad2c95c897fdc10fc483bd4107eec90';
const ROOT_URL = `http://api.openweathermap.org/data/2.5/forecast?APPID=${API_KEY}`;



export const fetchWeather = city => async dispatch => {
   const url = `${ROOT_URL}&q=${city},us`;
   const response = await axios.get(url);

   console.log("fetch weather",response)

   dispatch({ type: FETCH_WEATHER,
              payload: response  
   }); 

};




export const fetchBygeneSymbol = geneSymbol => async dispatch => {
  const response = await axios.get('/api/geneSymbol', {
    params: {
      geneSymbol: geneSymbol
    }
  });

  //const exampleData = [
  //  { geneName: 'Jane' },
  //  { geneName: 'Alex' },
  //  { geneName: 'Jim' }
  // ];

  //console.log('Received request', hgncID, exampleData);
  //dispatch({ type: FETCH_BYGENESYMBOL, payload: exampleData });
  console.log('Received request', geneSymbol, response);

    // only dispatch the content we have interest into the payload
    dispatch({ type: FETCH_BYGENESYMBOL, 
               payload: response.data.rows 
    });
};

//export const fetchUser = () => async dispatch => {
//     const res = await axios.get('/api/current_user');
//     
//     dispatch({
//             type: FETCH_USER,
//             payload: res.data
//          });
//};  


export const fetchAllStudies = () => async dispatch => {
  const res = await axios.get('/api/studies');
  
  dispatch({
          type: FETCH_ALL_STUDIES,
          payload: res.data.rows
       });
}; 

export const fetchGeneSummary = geneSymbol => async dispatch => {
  const response = await axios.get('/api/gene_summary', {
    params: {
      geneSymbol: geneSymbol
    }
  });
  
  console.log('Received gene summary action', geneSymbol, response);

  dispatch({
          type: FETCH_GENE_SUMMARY,
          payload: response.data.rows
       });
};  

export const fetchCancerDriverList = formValues => async dispatch => {
  const response = await axios.get('/api/cancer_driver_list', {
    params: {
      geneSymbol: formValues.geneSymbol,
      elementType: formValues.elementType,
      cancerType: formValues.cancerType,
      evidenceType: formValues.evidenceType,
      pmid: formValues.pmid
    }
  });
  
  console.log('Received fetchCancerDriverList action', formValues, response);

  dispatch({
          type: FETCH_CANCER_DRIVER_LIST,
          payload: response.data.rows
       });
};  


export const createSubmission = formData => async dispatch => {





  const response = await axios.post('/api/submissions', formData)
                       
  
  
  toast.success('upload success')
  

  console.log('Received createSubmission action', formData, response);

  dispatch({
          type: CREATE_SUBMISSION,
          payload: response
       });
};  