import { CREATE_SUBMISSION } from '../actions/types';

export default function(state = {}, action) {
  switch (action.type) {
    case CREATE_SUBMISSION:
      console.log('Action received (submission)', action);
      return { ...state, [action.payload.id]: action.payload };
    default:
      return state;
  }
}
