#include <iostream>
#include <utility>

/*
  a stupid code that generates automatically  the pattern of a
  "standard" finite element (just cut and past the result in the file
  feDef.h, and don't forget to delete the last "," ).
                                               Jean-Fred Gerbeau 11/99.
*/
int main()
{
  int nbNode  = 10;
  int nbPattern = nbNode*nbNode;
  int nbDiag  = nbNode;
  int nbUpper = nbNode*(nbNode-1)/2;
  int ip,i,j;
  std::pair<int,int> _pattern[nbPattern];
  // first : diagonal terms
  for(ip=1;ip<=nbDiag;ip++)
    _pattern[ip-1].first = _pattern[ip-1].second = ip;
  // second : upper terms and third : lower terms
  ip = nbDiag+1;
  for(i=1;i<nbNode;i++){
    for(j=i+1;j<=nbNode;j++){
      _pattern[ip-1].first   = _pattern[ip-1+nbUpper].second = i;
      _pattern[ip-1].second  = _pattern[ip-1+nbUpper].first  = j;
      ip++;
    }
  }

  std::cout << "// First : rows " << std::endl;
  std::cout << "static int _pf[" << nbPattern << "] = {";
  for(ip=0;ip<nbDiag;ip++) std::cout << _pattern[ip].first << "," ;
  std::cout << "// diag  entries" << std::endl;
  for(ip=nbDiag;ip<nbDiag+nbUpper;ip++) std::cout << _pattern[ip].first << "," ;
  std::cout << "// upper  entries" << std::endl;
  for(ip=nbDiag+nbUpper;ip<nbPattern;ip++)std::cout << _pattern[ip].first << "," ;
  std::cout << "}; // lower  entries" << std::endl << std::endl;

  std::cout << "// Second : column " << std::endl;
  std::cout << "static int _ps[" << nbPattern << "] = {";
  for(ip=0;ip<nbDiag;ip++) std::cout << _pattern[ip].second << "," ;
  std::cout << "// diag  entries" << std::endl;
  for(ip=nbDiag;ip<nbDiag+nbUpper;ip++) std::cout << _pattern[ip].second << "," ;
  std::cout << "// upper  entries" << std::endl;
  for(ip=nbDiag+nbUpper;ip<nbPattern;ip++)std::cout << _pattern[ip].second << "," ;
  std::cout << "}; // lower  entries" << std::endl;

}

