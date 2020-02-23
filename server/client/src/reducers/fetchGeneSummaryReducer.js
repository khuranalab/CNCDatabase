import { FETCH_GENE_SUMMARY } from '../actions/types';

export default function(state = [], action) {
  switch (action.type) {
    case FETCH_GENE_SUMMARY:
      //console.log('Action received', action);
      return action.payload;
    default:
      return state;
  }
}
