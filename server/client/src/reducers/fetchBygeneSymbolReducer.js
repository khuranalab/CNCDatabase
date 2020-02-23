import { FETCH_BYGENESYMBOL } from '../actions/types';

export default function(state = [], action) {
  switch (action.type) {
    case FETCH_BYGENESYMBOL:
      //console.log('Action received', action);
      return action.payload;
    default:
      return state;
  }
}
