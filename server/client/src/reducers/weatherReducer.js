import { FETCH_WEATHER } from '../actions/types';

export default function(state = [], action) {
  switch (action.type) {
    case FETCH_WEATHER:
      //console.log('Action received', action);
      return action.payload.data;
    default:
      return state;
  }
}