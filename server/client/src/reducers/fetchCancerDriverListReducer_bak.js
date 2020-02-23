import { FETCH_CANCER_DRIVER_LIST } from '../actions/types';

export default function(state = [], action) {
  switch (action.type) {
    case FETCH_CANCER_DRIVER_LIST:
      console.log('Action received(driver list)', action);
      return action.payload;
    default:
      return state;
  }
}
