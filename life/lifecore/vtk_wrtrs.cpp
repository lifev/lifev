#include "vtk_wrtrs.hpp"


void wr_vtk_ascii_scalar(string fname, string name, Real* U, int Usize,
			 string look_up_table)
{
 int i,j;
 ofstream ofile(fname.c_str(),ios::app);

 ASSERT(ofile,"Error: Output file cannot be open");
 
 ofile << "SCALARS " << name << " " << "float " << "1" << endl;
 ofile << "LOOKUP_TABLE " << look_up_table << endl;

 for (i=0;i<Usize;++i)
  ofile << U[i] << endl;

 if (look_up_table!="default")
   {
    ofile << "LOOKUP_TABLE " << look_up_table << " " << Usize << endl; 
    ifstream lutfile(look_up_table.c_str());
    ASSERT(lutfile,"Error: LookUp Table  file cannot be open"); 
    // assert
    Real r,g,b,a;
    for (j=0;j<Usize;++j){
     lutfile >> r >> g >> b >> a;
     ofile << r << " " << g << " " << b << " " << a << endl;   
    }
   } 
}
void wr_vtk_ascii_vector(string fname, string name, Real* U, int Usize)
{
 int i,j;
 ofstream ofile(fname.c_str(),ios::app);
 int nbcomp = 3;

 ASSERT(ofile,"Error: Output file cannot be open"); 


 ofile << "VECTORS " << name << " " << "float " << endl;
 for (i=0;i< Usize/nbcomp; ++i){
     for (j=0; j<nbcomp; j++)
       ofile << U[i+j*Usize/nbcomp] << " ";
     ofile << endl;
   }
}

//----------------------------------------------------------------------
// obsolete ?
void wr_vtk_ascii_scalar(string fname, string name, vector<Real> U, string look_up_table)
{
 unsigned int i,j;
 ofstream ofile(fname.c_str(),ios::app);

 ASSERT(ofile,"Error: Output file cannot be open"); 


 ofile << "POINT_DATA " << U.size() << endl;
 ofile << "SCALARS " << name << " " << "float " << "1" << endl;
 ofile << "LOOKUP_TABLE " << look_up_table << endl;

 for (i=0;i<U.size();++i)
  ofile << U[i] << endl;

 if (look_up_table!="default")
   {
    ofile << "LOOKUP_TABLE " << look_up_table << " " << U.size() << endl; 
    ifstream lutfile(look_up_table.c_str());
    ASSERT(lutfile,"Error: LookUp Table  file cannot be open"); 
    // assert
    Real r,g,b,a;
    for (j=0;j<U.size();++j){
     lutfile >> r >> g >> b >> a;
     ofile << r << " " << g << " " << b << " " << a << endl;   
    }
   } 
}

