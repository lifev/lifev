#define N_PHY_UNK 1
#define DEBUG

#define CSR_ORD // ensuring Ordered CSR

#include <vector>
#include <algorithm>
#include "lifeV.hpp"
#include "readMesh3D.hpp"
#include "chrono.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "refFE.hpp"
#include "elemOper_ext.hpp"
#include "dof.hpp"
#include "markers.hpp"

#include "values.hpp"
#include "assemb.hpp"

#include "dataAztec.hpp"

//===================================
//
// DIFFERENTIAL OPERATORS WE WANT TO SOLVE FOR, WRAPPED IN A SUITABLE CLASS
//
typedef EOExpr<Real,Stiff>  EOStiff;


// f
class SourceFct
{
public:
  inline Real operator()(Real x,Real y,Real z,int ic=0) const {
    return 0;
  }
};

/////////////////
// AZTEC's STUFF
////////////////
void init_options(int options[], double params[])
{

  /*
   * Choose among AZTEC options (see User's Guide).
   */

 AZ_defaults(options, params);

 options[AZ_solver]   = AZ_gmres;
 options[AZ_scaling]  = AZ_none;
 options[AZ_conv]     = AZ_r0;
 options[AZ_output]   = 1; //AZ_warnings;
 options[AZ_pre_calc] = AZ_calc;
 options[AZ_max_iter] = 1550;
 options[AZ_poly_ord] = 5;
 options[AZ_overlap]  = AZ_none;
 options[AZ_kspace]   = 40;
 options[AZ_aux_vec]  = AZ_resid;
 options[AZ_precond] =  AZ_dom_decomp;
 options[AZ_subdomain_solve] = AZ_ilut;
 options[AZ_keep_info] = 0;
 params[AZ_tol]       = 1.00e-10;
 params[AZ_drop]      = 1.00e-4;
 params[AZ_ilut_fill] = 5;
 params[AZ_omega]     = 1.;

} /* init_options */
//////////////////////////////////////////////////

