#include "bareItems.hpp"


/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from three Point ID's
  \param  bool is false if orientation  has been changed.
  \param i is a Point ID 
  \param j is a Point ID
  \param k is a Point ID

  To be used for triangular faces.
  \pre i, j and k >0. i!=j!=k

*/
pair<BareFace,bool>
makeBareFace(ID const i, ID  const j, ID const k) {
  if(i<j && i < k) {
    if(j<k){
      return pair<BareFace,bool>(BareFace(i,j,k),true);
    }
    else{
      return pair<BareFace,bool>(BareFace(i,k,j),false);
    }
  }
  else if(j<k && j < i) {
    if(k<i){
      return pair<BareFace,bool>(BareFace(j,k,i),true);
    }
    else{
      return pair<BareFace,bool>(BareFace(j,i,k),false);
    }
  }
  else 
    {
      if(i<j){
	return pair<BareFace,bool>(BareFace(k,i,j),true);
      }
      else{
	return pair<BareFace,bool>(BareFace(k,j,i),false);
      }
    }
}

/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from four Point ID's
  \param  bool is false if orientation  has been changed.
  \param i is a Point ID 
  \param j is a Point ID
  \param k is a Point ID
  \param l is a Point ID
  \pre i, j, k and l >0. i!=j!=k!=l
  
  To be used with Quad faces.

\remarks For quad faces the construction process is more complex. We start from
  the smallest vertex and we take the first three vertices in the
   sequence. We then procede as for the triangles.
*/

pair<BareFace,bool>
makeBareFace(ID const i, ID  const j, ID const k, ID const l) {
  vector<ID> helper(4);
  helper[0]=i;  helper[1]=j;  helper[2]=k;  helper[3]=l;
  vector<ID>::iterator vi=max_element(helper.begin(),helper.end());
  rotate(helper.begin(),vi,helper.end());
  return makeBareFace(helper[1],helper[2],helper[3]);
}
