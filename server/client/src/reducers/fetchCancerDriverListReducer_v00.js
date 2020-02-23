import { FETCH_CANCER_DRIVER_LIST } from '../actions/types';
import axios from 'axios';

async function getInitialCancerDriverList() {
  try {
    // fetch data from a url endpoint
    const response = await axios.get('/api/cancer_driver_list', {
      params: {
        geneSymbol: [],
        elementType: [],
        cancerType: [],
        evidenceType: [],
        pmid: []
      }
    });

    const data = await response.data.rows;

    return data;

  } catch (error) {
    alert(error); // catches both errors
  }
}

const initialCancerDriverList = [{
  type: FETCH_CANCER_DRIVER_LIST,
  payload: getInitialCancerDriverList()
}]

//const initialCancerDriverList = dispatch(FETCH_CANCER_DRIVER_LIST)

export default function(state = [], action) {
  switch (action.type) {
    case FETCH_CANCER_DRIVER_LIST:
      console.log('Action received(driver list)', action);
      return action.payload;
    default:
      return state;
  }
}
