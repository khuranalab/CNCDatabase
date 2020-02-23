import { FETCH_ALL_STUDIES } from '../actions/types';

export default function(state = [], action) {
  switch (action.type) {
    case FETCH_ALL_STUDIES:
      //console.log('Action received', action);
      return action.payload;
    default:
      return state;
  }
}