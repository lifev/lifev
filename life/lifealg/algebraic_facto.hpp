/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file algebraic_facto.h
  \author Alain Gauthier.
  \date 04/02

  \brief Matrices-vector product for AZTEC. Preconditioner for AZTEC

  \version 18/09/02, addition of: other pointers in DataFactorisation, two other Schur
  complement preconditionner my_precSchur_PC anc my_precSchur_CC, Alain.

  \version 18/01/03, fully Dirichlet bc conditions treatment (without define). New treatment of the
  Aztec options for the Schur complement system using DataAztec (without including main.h). Miguel

  #purposes: definition of various functions which are used as
            matrix-vector products with AZTEC. The aim is the implementation
            of exact/inexact algebraic factorisation methods for NS equation

             construction of preconditioners used within AZTEC according
            to the algebraic factorisation method.
*/

#ifndef _ALGEBRAIC_FACTO_HH
#define _ALGEBRAIC_FACTO_HH

#ifndef DATA_NAME_AZTEC
#define DATA_NAME_AZTEC 1
#endif

#ifndef DATA_NAME_AZTEC1
#define DATA_NAME_AZTEC1 10
#endif
#ifndef DATA_NAME_AZTEC2
#define DATA_NAME_AZTEC2 20
#endif
#ifndef DATA_NAME_AZTEC3
#define DATA_NAME_AZTEC3 30
#endif
#ifndef DATA_NAME_AZTEC_MATS
#define DATA_NAME_AZTEC_MATS 40
#endif

#include "sparseArray.hpp"
#include "dataAztec.hpp"

namespace LifeV
{
//////////////////////
//!class for passing matrix factorisation information to AZTEC
/*!
   contains pointers to matrices and other informations which are
   given to AZTEC and used in the matrix-vector product.
 */
template <typename MatrixTypeC, typename MatrixTypeD, typename MatrixTypeDtr,
	  typename MatrixTypeH, typename MatrixTypeMpLp, typename VectorType>
class DataFactorisation
{
public:

  //! Constructor which points on external data, version for Cahouet-Chabart
  /*!
    \param ex_C       Bloc C of the Stokes matrix
    \param ex_D       Bloc D of the Stokes matrix
    \param ex_trD     Bloc D^T of the Stokes Matrix
    \param ex_H       Approximation of C, can be stored in a Vector type if diagonal
    \param ex_HinvC   Matrix H^{-1}*C
    \param ex_HinvDtr Matrix H^{-1}*D^T
    \param ex_Mp      Matrix of pressure mass
    \param ex_vec     Vector which can be use to exchange informations with AZTEC matrix-vector product
    \param ex_idFacto  identifier equal to 1 if a former factorization of Lp was already performed, equal to 0 otherwise
    \param ex_mu      copy of the viscosity value (not always used)
    \param dataAztec_i  points to a DataAztec object holding options and parameters for the linear system
    \param dataAztec_s  points to a DataAztec object holding options and parameters for the linear system
    \param fullEssential true if full Dirichlet conditions are involved
   */
  DataFactorisation(MatrixTypeC const & ex_C, MatrixTypeD const & ex_D,
  		    MatrixTypeDtr const & ex_trD, MatrixTypeH const & ex_H,
		    MatrixTypeC const & ex_HinvC, MatrixTypeDtr const & ex_HinvDtr,
		    MatrixTypeMpLp const & ex_Mp, VectorType & ex_vec,
		    UInt & ex_idFacto, double const ex_mu, DataAztec& dataAztec_i,
		    DataAztec& dataAztec_s, bool fullEssential, const UInt recur);

  //! Constructor which points on external data, version for ACT pressure-corrected
  /*!
    \param ex_C       Bloc C of the Stokes matrix
    \param ex_D       Bloc D of the Stokes matrix
    \param ex_trD     Bloc D^T of the Stokes Matrix
    \param ex_H       Approximation of C, can be stored in a Vector type if diagonal
    \param ex_HinvC   Matrix H^{-1}*C
    \param ex_HinvDtr Matrix H^{-1}*D^T
    \param ex_vec     Vector which can be use to exchange informations with AZTEC matrix-vector product
    \param dataAztec_i  points to a DataAztec object holding options and parameters for the linear system
    \param dataAztec_s  points to a DataAztec object holding options and parameters for the linear system
    \param fullEssential true if full Dirichlet conditions are involved
   */
  DataFactorisation(MatrixTypeC const & ex_C, MatrixTypeD const & ex_D,
  		    MatrixTypeDtr const & ex_trD, MatrixTypeH const & ex_H,
		    MatrixTypeC const & ex_HinvC, MatrixTypeDtr const & ex_HinvDtr,
		    VectorType & ex_vec, DataAztec& dataAztec_i,
		    DataAztec& dataAztec_s, bool fullEssential, const UInt recur);


  //! Constructor which points on external data, version for pressure matrix method
  //------------NOT USED YET--------------
  /*!
    \param ex_C       Bloc C of the Stokes matrix
    \param ex_D       Bloc D of the Stokes matrix
    \param ex_trD     Bloc D^T of the Stokes Matrix
    \param ex_H       Approximation of C, can be stored in a Vector type if diagonal
    \param ex_Lp      Matrix of Laplacian of pressure (useful for preconditioning S=D*H^{-1}D^T)
    \param ex_vec     Vector which can be use to exchange informations with AZTEC matrix-vector product
    \param ex_idFacto  identifier equal to 1 if a former factorization of Lp was already performed, equal to 0 otherwise
   \param dataAztec_i  points to a DataAztec object holding options and parameters for the linear system
   \param dataAztec_s  points to a DataAztec object holding options and parameters for the linear system
   \param fullEssential true if full Dirichlet conditions are involved
   */
  DataFactorisation(MatrixTypeC const & ex_C, MatrixTypeD const & ex_D,
  		    MatrixTypeDtr const & ex_trD, MatrixTypeH const & ex_H,
		    MatrixTypeMpLp const & ex_Lp, VectorType & ex_vec,
		    UInt & ex_idFacto, DataAztec& dataAztec_i,
		    DataAztec& dataAztec_s, bool fullEssential, const UInt recur);

  //! Constructor which points on external data, version for pressure matrix method
  //------------ will be replaced by the previous one ------------
  /*!
    \param ex_C       Bloc C of the Stokes matrix
    \param ex_D       Bloc D of the Stokes matrix
    \param ex_trD     Bloc D^T of the Stokes Matrix
    \param ex_H       Approximation of C, can be stored in a Vector type if diagonal
    \param ex_vec     Vector which can be use to exchange informations with AZTEC matrix-vector product
    \param dataAztec_i  points to a DataAztec object holding options and parameters for the linear system
    \param dataAztec_s  points to a DataAztec object holding options and parameters for the linear system
    \param fullEssential true if full Dirichlet conditions are involved
   */
  DataFactorisation(MatrixTypeC const & ex_C, MatrixTypeD const & ex_D,
  		    MatrixTypeDtr const & ex_trD, MatrixTypeH const & ex_H,
		    VectorType & ex_vec, DataAztec& dataAztec_i,
		    DataAztec& dataAztec_s, bool fullEssential, const UInt recur);

  //! point to matrix C or an equivalent approximation
  const MatrixTypeC* _C;

  //! point to matrix D
  const MatrixTypeD* _D;

  //! point to matrix trD
  const MatrixTypeDtr* _trD;

  //! point to matrix H: a diagonal approximation of C easy to invert and
  //! useful for preconditioning
  const MatrixTypeH* _H;

  //! point to matrix H^{-1}*C
  const MatrixTypeC* _HinvC;

  //! point to matrix H^{-1}*D^T
  const MatrixTypeDtr* _HinvDtr;

  //! point to pressure mass matrix Mp
  const MatrixTypeMpLp* _Mp;

  //! point to laplacian of pressure mass matrix Lp
  const MatrixTypeMpLp* _Lp;

  //! point to a vector which will contain eventually in output
  //! the product C^{-1}*D^T*P
  VectorType* _vec;

  //! idFacto which allow to record that a factorisation of Lp
  //! has been was already done
  UInt* _idFacto;

  //! a copy of the viscosity, used mainly for Cahouet-Chabart preconditioner
  double _mu;

  //! points to a DataAztec object holding options and parameters for the linear system
  DataAztec* _dataAztec_i;

  //! points to a DataAztec object holding options and parameters for the linear system
  DataAztec* _dataAztec_s;

  //! true if full Dirichlet conditions are involved
  bool _fullEssential;

  //! recursion level
  UInt _recur;

};

template <typename MatrixTypeC, typename MatrixTypeD, typename MatrixTypeDtr,
	  typename MatrixTypeH, typename MatrixTypeMpLp, typename VectorType>
DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>
::DataFactorisation(MatrixTypeC const & ex_C, MatrixTypeD const & ex_D,
  		    MatrixTypeDtr const & ex_trD, MatrixTypeH const & ex_H,
		    MatrixTypeC const & ex_HinvC, MatrixTypeDtr const & ex_HinvDtr,
		    MatrixTypeMpLp const & ex_Mp, VectorType & ex_vec,
		    UInt & ex_idFacto, double const ex_mu, DataAztec& dataAztec_i,
		    DataAztec& dataAztec_s, bool fullEssential, const UInt recur=1) {
    _C=&ex_C;
    _D=&ex_D;
    _trD=&ex_trD;
    _H=&ex_H;
    _HinvC=&ex_HinvC;
    _HinvDtr=&ex_HinvDtr;
    _Mp=&ex_Mp;
    _vec=&ex_vec;
    _idFacto=&ex_idFacto;
    _mu=ex_mu;
    _dataAztec_i=&dataAztec_i;
    _dataAztec_s=&dataAztec_s;
    _fullEssential = fullEssential;
    _recur=recur;
}


template <typename MatrixTypeC, typename MatrixTypeD, typename MatrixTypeDtr,
	  typename MatrixTypeH, typename MatrixTypeMpLp, typename VectorType>
DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>
::DataFactorisation(MatrixTypeC const & ex_C, MatrixTypeD const & ex_D,
  		    MatrixTypeDtr const & ex_trD, MatrixTypeH const & ex_H,
		    MatrixTypeC const & ex_HinvC, MatrixTypeDtr const & ex_HinvDtr,
		    VectorType & ex_vec, DataAztec& dataAztec_i,
		    DataAztec& dataAztec_s, bool fullEssential, const UInt recur=1) {
    _C=&ex_C;
    _D=&ex_D;
    _trD=&ex_trD;
    _H=&ex_H;
    _HinvC=&ex_HinvC;
    _HinvDtr=&ex_HinvDtr;
    _vec=&ex_vec;
    _idFacto=0;
    _mu = 1.;
    _dataAztec_i=&dataAztec_i;
    _dataAztec_s=&dataAztec_s;
    _fullEssential = fullEssential;
    _recur=recur;
}


template <typename MatrixTypeC, typename MatrixTypeD, typename MatrixTypeDtr,
	  typename MatrixTypeH, typename MatrixTypeMpLp, typename VectorType>
DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>
::DataFactorisation(MatrixTypeC const & ex_C, MatrixTypeD const & ex_D,
  		    MatrixTypeDtr const & ex_trD, MatrixTypeH const & ex_H,
		    MatrixTypeMpLp const & ex_Lp, VectorType & ex_vec,
		    UInt & ex_idFacto, DataAztec& dataAztec_i,
		    DataAztec& dataAztec_s, bool fullEssential, const UInt recur=1) {
    _C=&ex_C;
    _D=&ex_D;
    _trD=&ex_trD;
    _H=&ex_H;
    _Lp=&ex_Lp;
    _vec=&ex_vec;
    _idFacto=&ex_idFacto;
    _mu=1.;
    _dataAztec_i=&dataAztec_i;
    _dataAztec_s=&dataAztec_s;
    _fullEssential=fullEssential;
    _recur=recur;
}


template <typename MatrixTypeC, typename MatrixTypeD, typename MatrixTypeDtr,
	  typename MatrixTypeH, typename MatrixTypeMpLp, typename VectorType>
DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>
::DataFactorisation(MatrixTypeC const & ex_C, MatrixTypeD const & ex_D,
  		    MatrixTypeDtr const & ex_trD, MatrixTypeH const & ex_H,
		    VectorType & ex_vec, DataAztec& dataAztec_i,
		    DataAztec& dataAztec_s, bool fullEssential, const UInt recur=1) {
    _C=&ex_C;
    _D=&ex_D;
    _trD=&ex_trD;
    _H=&ex_H;
    _vec=&ex_vec;
    _mu=1.;
    _dataAztec_i=&dataAztec_i;
    _dataAztec_s=&dataAztec_s;
    _fullEssential=fullEssential;
    _recur=recur;
}

////////////////////////
/*!
   \brief Matrix-vector product function passed through AZTEC.

   Computes the product C * p where
   the matrix C is stored in the class DataFactorisation.

   \param p
   \param ap
   \param Amat
   \param proc_config contains information on processor for AZTEC, initialized by the function AZ_set_proc_config.

   In output the result is set in the vector ap.
*/
template <typename MatrixTypeC,
  typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
void my_Cmatvec(double *p, double *ap, AZ_MATRIX * Amat,
		     int proc_config[])
{
  // Extraction of informations stored in the structure AZ_MATRIX
  DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>*
    my_data = static_cast< DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,
    MatrixTypeH,MatrixTypeMpLp,VectorType>* >(AZ_get_matvec_data(Amat));

  // product C*p:
  operMatVec(ap,*(my_data->_C),p);
}

////////////////////////
/*!
   Supply local matrix (without ghost node columns) for rows given by
   requested_rows[0 ... N_requested_rows-1].  Return this information in
   'row_lengths, columns, values'.  If there is not enough space to complete
   this operation, return 0. Otherwise, return 1.

   Purpose : necessary for using AZTEC preconditioner with USER-DEFINED
             matrix.

   \param Amat             On input, points to user's data containing
                           matrix values.
   \param N_requested_rows On input, number of rows for which nonzero are
                           to be returned.
   \param requested_rows   On input, requested_rows[0...N_requested_rows-1]
                           give the row indices of the rows for which nonzero
			   values are returned.
   \param row_lengths      On output, row_lengths[i] is the number of
                           nonzeros in the row 'requested_rows[i]'
   \param columns,
   \param values           On output, columns[k] and values[k] contains
                           the column number and value of a matrix nonzero
                           where all nonzeros for requested_rows[i] appear
			   before requested_rows[i+1]'s nonzeros.
			   NOTE: Arrays are of size 'allocated_space'.
   \param allocated_space  On input, indicates the space available in
                           'columns' and 'values' for storing nonzeros.
			   If more space is needed, return 0.
   \return An integer value, 0 if more memory space is needed, 1 otherwise.

   REMARK: VERSION VALID FOR MixedMatr MATRIX, it's not necessary for
           other matrix types (C would be an MSR matrix in this case).
 */
template <typename MatrixTypeC,
  typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
int  matC_getrow( int columns[], double values[], int row_lengths[],
		  struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
		  int requested_rows[], int allocated_space)
{
  // Extraction of informations stored in the structure AZ_MATRIX
  DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,
    MatrixTypeMpLp,VectorType>* my_data = static_cast< DataFactorisation<MatrixTypeC,MatrixTypeD,
	      MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>* >(AZ_get_matvec_data(Amat));

  UInt nz_ptr, row, jcol, nnz;
  UInt m=0,n=0;
  int k;
  nnz=0;
  Container coldata, position;

  coldata.resize(my_data->_C->Patt()->nCols());
  position.resize(my_data->_C->Patt()->nCols());

  // loop on requested rows
  nz_ptr= 0;
  for (k = 0; k < N_requested_rows; k++) {

    row = requested_rows[k];

    nnz=my_data->_C->Patt()->row(row, coldata.begin(), position.begin());

    if ( (int)(nz_ptr + nnz) > allocated_space) return(0);

    //TODO : (Alain) I should define a method of the matrix C which
    //               would return the values so we woulden't have to
    //               access explicitely to values or bvalues like here !
    for (jcol= 0; jcol < nnz; jcol++)
      {
	extract_pair(my_data->_C->Patt()->
		     locateElBlock(row, coldata[jcol]), m, n);
	values[nz_ptr]= my_data->_C->bValues(m,n)[position[jcol]];
	columns[nz_ptr++]= coldata[jcol];
      }

    row_lengths[k] = nnz;
   }
   return(1);
}

////////////////////////
/*!
   \brief Matrix-vector product function passed through AZTEC.

   The aim is to compute the product (D*C^(-1)*trD) * p where
   the matrices C, D and trD are stored in the class DataFactorisation.
   In output the result is set in the vector ap.

   REMARK: this function correspond to the EXACT algebraic factorisation
           method.
 */
template <typename MatrixTypeC,
  typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
void my_matvec(double *p, double *ap, AZ_MATRIX * Amat, int proc_config[])
{
  // Extraction of C, D and trD stored in the structure AZ_MATRIX
  DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,
    MatrixTypeMpLp,VectorType>* my_data = static_cast< DataFactorisation<MatrixTypeC,MatrixTypeD,
	      MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,
	      VectorType>* >(AZ_get_matvec_data(Amat));

  // introduction of vectors
  UInt dim  = my_data->_trD->Patt()->nCols();
  UInt dim2 = my_data->_trD->Patt()->nRows();

  // product C^(-1)*p1= p2 where p2= D^T * p:

  // AZTEC specifications
  int    proc_config_i[AZ_PROC_SIZE];// Processor information:
  //  proc_config[AZ_node] = node name
  //  proc_config[AZ_N_procs] = # of nodes
  int    options_i[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_i[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_i[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.


  AZ_set_proc_config(proc_config_i, AZ_NOT_MPI);

  //AZTEC matrix and preconditioner
  AZ_MATRIX *A_i;
  AZ_PRECOND *prec_i;

  int N_eq_i= dim2;
  A_i= AZ_matrix_create(N_eq_i);
  // data containing the matrices C, D, trD and H as pointers
  // are passed through A_i and prec_i:
  AZ_set_MATFREE(A_i, my_data,
		 my_Cmatvec<MatrixTypeC,
		 MatrixTypeD,
		 MatrixTypeDtr,
		 MatrixTypeH,
		 MatrixTypeMpLp,
		 VectorType>);
  // necessary for almost all preconditioner construted by AZTEC:

  AZ_set_MATFREE_getrow(A_i,my_data,
			matC_getrow<MatrixTypeC,
			MatrixTypeD,
			MatrixTypeDtr,
			MatrixTypeH,
			MatrixTypeMpLp,
			VectorType>,NULL,0,proc_config_i);

  prec_i= AZ_precond_create(A_i, AZ_precondition, NULL);

  my_data->_dataAztec_i->aztecOptionsFromDataFile(options_i,params_i);

  //initialisation of first recursion level for AZTEC memory manager
  options_i[AZ_recursion_level]=my_data->_recur;
  //without writing output
  options_i[AZ_output]=AZ_none;

  //recall options used for the first linear system solving (i)
  options_i[AZ_pre_calc]=AZ_reuse;
  A_i->data_org[AZ_name]= DATA_NAME_AZTEC;
  options_i[AZ_keep_info]=1; // keep information

  std::vector<double> p1(A_i->data_org[AZ_N_internal]+A_i->data_org[AZ_N_border]+A_i->data_org[AZ_N_external]);
  std::vector<double> p2(A_i->data_org[AZ_N_internal]+A_i->data_org[AZ_N_border]);
  p1.assign( p1.size(), 0.0 );
  p2.assign( p1.size(), 0.0 );

  // product trD*p:
  operMatVec(p1, *(my_data->_trD), p);

  // solve the system C*p2= p1 with AZTEC.
  AZ_iterate( &p2.front(), &p1.front(), options_i, params_i, status_i,
	     proc_config_i, A_i, prec_i, NULL);

  //return to zero recursion level for AZTEC memory manager
  //--options_i[AZ_recursion_level];

  // destroy A_i
  AZ_matrix_destroy(&A_i);

  // Result p2 is saved and will be used to correct the velocity
  // in the step (iii) of algebraic factorization
  for (UInt i=0; i< dim2; i++)
    (*(my_data->_vec))[i]= p2[i];

  // product p=D*p2:
  operMatVec(ap, *my_data.D(), &p2.front());

  if (my_data->_fullEssential)
    ap[dim-1]=p[dim-1]; // diagonalisation of the last row.
}

////////////////////////
/*!
   \brief Matrix-vector product function passed through AZTEC.

   SPECIALISATION IF C IS OF TYPE MSR

   The aim is to compute the product (D*C^(-1)*trD) * p where
   the matrices C, D and trD are stored in the class DataFactorisation.
   In output the result is set in the vector ap.

   REMARK: this function correspond to the EXACT algebraic factorisation
           method.
 */
template <typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
void my_matvec(double *p, double *ap, AZ_MATRIX * Amat, int proc_config[])
{
  // Extraction of C, D and trD stored in the structure AZ_MATRIX
  DataFactorisation<MSRMatr<double>,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,
    MatrixTypeMpLp,VectorType>* my_data=static_cast< DataFactorisation<MSRMatr<double>,MatrixTypeD,
	      MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,
	      VectorType>* >(AZ_get_matvec_data(Amat));

  // introduction of vectors
  UInt dim_u = my_data->_trD->Patt()->nRows();

  // product C^(-1)*p1= p2 where p2= D^T * p:
  // AZTEC specifications
  int    data_org_i[AZ_COMM_SIZE];   // data organisation
  int    proc_config_i[AZ_PROC_SIZE];// Processor information:
  int    options_i[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_i[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_i[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.


  AZ_set_proc_config(proc_config_i, AZ_NOT_MPI);

  // data_org assigned "by hands" while no parallel computation is performed
  data_org_i[AZ_N_internal]= dim_u;
  data_org_i[AZ_N_border]= 0;
  data_org_i[AZ_N_external]= 0;
  data_org_i[AZ_N_neigh]= 0;
  data_org_i[AZ_name]= DATA_NAME_AZTEC;

  //AZTEC matrix and preconditioner
  AZ_MATRIX *A_i;
  AZ_PRECOND *prec_i;

  int N_eq_i= dim_u;

  A_i= AZ_matrix_create(N_eq_i);
  AZ_set_MSR(A_i, (int*) my_data->_C->Patt()->giveRaw_bindx(),
	     (double*) my_data->_C->giveRaw_value(),
	     data_org_i, 0, NULL, AZ_LOCAL);

  prec_i= AZ_precond_create(A_i, AZ_precondition, NULL);

  my_data->_dataAztec_i->aztecOptionsFromDataFile(options_i,params_i);

  //initialisation of first recursion level for AZTEC memory manager
  options_i[AZ_recursion_level]=my_data->_recur;


  //without writing output
  options_i[AZ_output]=AZ_none;

  //recall options used for the first linear system solving (i)
  options_i[AZ_pre_calc]=AZ_reuse;
  //A_i->data_org[AZ_name]= DATA_NAME_AZTEC;
  options_i[AZ_keep_info]=1; // keep information


  std::vector<double> p1( dim_u );
  std::vector<double> p2( dim_u );
  p1.assign( dim_u, 0.0 );
  p2.assign( dim_u, 0.0 );

  // product trD*p:
  operMatVec( &p1.front(), *(my_data->_trD), p);

  // solve the system C*p2= p1 with AZTEC.
  AZ_iterate( &p2.front(), &p1.front(), options_i, params_i, status_i,proc_config_i, A_i, prec_i, NULL);



  // destroy A_i and prec_i
  AZ_matrix_destroy(&A_i);
  AZ_precond_destroy(&prec_i);

  // Result p2 is saved and will be used to correct the velocity
  // in the step (iii) of algebraic factorization
  for (UInt i=0; i< dim_u; i++)
      (*(my_data->_vec))[i]= p2[i];

  // product p=D*p2:
  operMatVec(ap, *(my_data->_D), &p2.front());

  if ( my_data->_fullEssential ) {
    const UInt dim  = my_data->_trD->Patt()->nCols();
    ap[dim-1]=p[dim-1]; // diagonalisation of the last row.
  }


}

////////////////////////
/*!
   \brief Matrix-vector product function passed through AZTEC.

   SOLVE THE SYSTEM FOR C FOR EACH BLOCK SEPARATLY (WORK WITH MIXEDMATR TYPE)

   The aim is to compute the product (D*C^(-1)*trD) * p where
   the matrices C, D and trD are stored in the class DataFactorisation.
   In output the result is set in the vector ap.

   REMARK: this function correspond to the EXACT algebraic factorisation
           method.
 */
template <typename MatrixTypeC,
  typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
void my_matvec_block(double *p, double *ap, AZ_MATRIX * Amat,
		     int proc_config[])
{
  // Extraction of C, D and trD stored in the structure AZ_MATRIX
  DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,
    MatrixTypeMpLp,VectorType>* my_data = static_cast< DataFactorisation<MatrixTypeC,MatrixTypeD,
	      MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,
	      VectorType>* >(AZ_get_matvec_data(Amat));

  // introduction of vectors
  const UInt dim2 = my_data->_trD->Patt()->nRows();

  std::vector<double> p1( dim2 );
  std::vector<double> p2( dim2 );
  p1.assign( dim2, 0.0 );
  p2.assign( dim2, 0.0 );

  // product trD*p:
  operMatVec( &p1.front(), *(my_data->_trD), p);

  // product C^(-1)*p1= p2:
  // solve the system C*p2= p1 with AZTEC.

  // AZTEC specifications for each system
  int    data_org_i1[AZ_COMM_SIZE];   // data organisation for C1
  int    data_org_i2[AZ_COMM_SIZE];   // data organisation for C2
  int    data_org_i3[AZ_COMM_SIZE];   // data organisation for C3
  // identical for each system
  int    proc_config_i[AZ_PROC_SIZE];// Processor information:
  int    options_i[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_i[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_i[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.

  //
  AZ_set_proc_config(proc_config_i, AZ_NOT_MPI);

  //AZTEC matrix and preconditioner
  AZ_MATRIX *C1, *C2, *C3;
  AZ_PRECOND *prec_C1, *prec_C2, *prec_C3;

  int N_eq_i= dim2/3; // number of DOF for each component

  // data_org assigned "by hands" while no parallel computation is performed
  data_org_i1[AZ_N_internal]= N_eq_i;
  data_org_i1[AZ_N_border]= 0;
  data_org_i1[AZ_N_external]= 0;
  data_org_i1[AZ_N_neigh]= 0;
  data_org_i1[AZ_name]= DATA_NAME_AZTEC1;

  data_org_i2[AZ_N_internal]= N_eq_i;
  data_org_i2[AZ_N_border]= 0;
  data_org_i2[AZ_N_external]= 0;
  data_org_i2[AZ_N_neigh]= 0;
  data_org_i2[AZ_name]= DATA_NAME_AZTEC2;

  data_org_i3[AZ_N_internal]= N_eq_i;
  data_org_i3[AZ_N_border]= 0;
  data_org_i3[AZ_N_external]= 0;
  data_org_i3[AZ_N_neigh]= 0;
  data_org_i3[AZ_name]= DATA_NAME_AZTEC3;

  // set each block
  C1= AZ_matrix_create(N_eq_i);
  C2= AZ_matrix_create(N_eq_i);
  C3= AZ_matrix_create(N_eq_i);

  AZ_set_MSR(C1, (int*) my_data->_C->Patt()->block_ptr(0,0)->giveRaw_bindx(),
	     (double*) my_data->_C->giveRaw_value(0,0),
	     data_org_i1, 0, NULL, AZ_LOCAL);
  AZ_set_MSR(C2, (int*) my_data->_C->Patt()->block_ptr(1,1)->giveRaw_bindx(),
	     (double*) my_data->_C->giveRaw_value(1,1),
	     data_org_i2, 0, NULL, AZ_LOCAL);
  AZ_set_MSR(C3, (int*) my_data->_C->Patt()->block_ptr(2,2)->giveRaw_bindx(),
	     (double*) my_data->_C->giveRaw_value(2,2),
	     data_org_i3, 0, NULL, AZ_LOCAL);

  // create preconditioner for each block
  prec_C1= AZ_precond_create(C1, AZ_precondition, NULL);
  prec_C2= AZ_precond_create(C2, AZ_precondition, NULL);
  prec_C3= AZ_precond_create(C3, AZ_precondition, NULL);

  my_data->_dataAztec_i->aztecOptionsFromDataFile(options_i,params_i);

  //initialisation of first recursion level for AZTEC memory manager
  options_i[AZ_recursion_level]=my_data->_recur;
  //without writing output
  options_i[AZ_output]=AZ_none;

  //recall options used for the first linear system solving (i)
  options_i[AZ_pre_calc]=AZ_reuse;
  options_i[AZ_keep_info]=1; // keep information

  //for each block
  AZ_iterate( &p2.front(), &p1.front(), options_i, params_i, status_i,
	     proc_config_i, C1, prec_C1, NULL);
  AZ_iterate( &p2.front()+N_eq_i, &p1.front()+N_eq_i, options_i, params_i,
	     status_i, proc_config_i, C2, prec_C2, NULL);
  AZ_iterate( &p2.front()+2*N_eq_i, &p1.front()+2*N_eq_i, options_i,
	     params_i, status_i, proc_config_i, C3, prec_C3, NULL);

  //return to zero recursion level for AZTEC memory manager
  //--options_i[AZ_recursion_level];

  // destroy Ci and prec_Ci
  AZ_matrix_destroy(&C1);
  AZ_precond_destroy(&prec_C1);
  AZ_matrix_destroy(&C2);
  AZ_precond_destroy(&prec_C2);
  AZ_matrix_destroy(&C3);
  AZ_precond_destroy(&prec_C3);

  // Result p2 is saved and will be used to correct the velocity
  // in the step (iii) of algebraic factorization
  for (UInt i=0; i< dim2; i++)
    (*(my_data->_vec))[i]= p2[i];

  // product p=D*p2:
  operMatVec( ap, *(my_data->_D), &p2.front());

  if ( my_data->_fullEssential ) {
    const UInt dim  = my_data->_trD->Patt()->nCols();
    ap[dim-1]=p[dim-1]; // diagonalisation of the last row.
  }
}

////////////////////////
/*!
   \brief Matrix-vector product function passed through AZTEC.

   The aim is to compute the product (D*H^{-1}*trD) * p where
   the matrices H, D and trD are stored in the class DataFactorisation.
   In output the result is set in the vector ap.
   H corresponds to an approximation of the matrix C.

   REMARK: H is diagonal and then stored in vector type.
 */
template <typename MatrixTypeC,
  typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
void my_approxmatvec(double *p, double *ap, AZ_MATRIX * Amat, int proc_config[])
{
  // Extraction of H, , HinvDtr, D and trD stored in the structure AZ_MATRIX
  DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,
    MatrixTypeMpLp,VectorType>* my_data = static_cast< DataFactorisation<MatrixTypeC,MatrixTypeD,
	      MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,
	      VectorType>* >(AZ_get_matvec_data(Amat));

  // introduction of vectors
  UInt dim2 = my_data->_trD->Patt()->nRows();

  std::vector<double> p1(dim2);

  // product p1= H^{-1}*trD * p :
  operMatVec( &p1.front(), *(my_data->_HinvDtr), p);

  // product D*p1 :
  operMatVec(ap, *(my_data->_D), &p1.front());


  if(my_data->_fullEssential) {
    UInt dim  = my_data->_trD->Patt()->nCols();
    ap[dim-1]=p[dim-1]; // diagonalisation of the last row.
  }
}

////////////////////////
/*!
   \brief Preconditioner of the Schur complement passed to AZTEC.

   The aim is to solve the system (D*H^{-1}*trD) * z = r where
   in input z=r.
   the matrices H, D and trD are stored in the class DataFactorisation.
   In output the result is set in the vector z.
   H corresponds to an approximation of the matrix C in the case of
   Chorin-Temam's method, Yosida one or inexact LU precondioning.

   REMARK: H is diagonal and then stored in vector type.
 */
template <typename MatrixTypeC,
  typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
void my_precSchur(double *z, int *options, int *proc_config, double *params,
		  AZ_MATRIX *Amat, AZ_PRECOND *prec)
{
  // Extraction of H, D and trD stored in the structure AZ_MATRIX
  DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,
    MatrixTypeMpLp,VectorType>* my_data =static_cast< DataFactorisation<MatrixTypeC,MatrixTypeD,
    MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>* >(AZ_get_precond_data(prec));

  UInt dim  = my_data->_trD->Patt()->nCols();
  double lastz=z[dim-1];

  // AZTEC specifications
  int    options_s[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_s[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_s[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                    // indicating success or failure.
  //

  //AZTEC matrices
  AZ_MATRIX *SchurPrec;

  int N_eq= dim;
  SchurPrec= AZ_matrix_create(N_eq);
  // data containing the matrices C, H and trD as pointers
  // are passed through S :
  AZ_set_MATFREE(SchurPrec, my_data, my_approxmatvec<
		 MatrixTypeC,
		 MatrixTypeD,
		 MatrixTypeDtr,
		 MatrixTypeH,
		 MatrixTypeMpLp,
		 VectorType>);

  my_data->_dataAztec_s->aztecOptionsFromDataFile(options_s,params_s);

  //initialisation of first recursion level for AZTEC memory manager
  options_s[AZ_recursion_level]=my_data->_recur;

  AZ_iterate(z, z, options_s, params_s, status_s,
	     proc_config, SchurPrec, NULL, NULL);

  //return to zero recursion level for AZTEC memory manager
  //--options_s[AZ_recursion_level];

  // destroy AZTEC matrices
  AZ_matrix_destroy(&SchurPrec);

  if (my_data->_fullEssential) {
    z[dim-1]=lastz; // diagonalisation of the last row.
  }
}

////////////////////////
/*!
   \brief Pressure corrected Preconditioner of the Schur complement passed to AZTEC.

   The aim is to solve the system (S*B^{-1}*S) * z = r where
   S=(D*H^{-1}*trD), B= D*H^{-1}*C*H^{-1}*D^T and
   in input z=r.
   the matrices H, D, HinvC, HinvDtr and trD are stored in the
   class DataFactorisation.
   In output the result is set in the vector z.
   H corresponds to an approximation of the matrix C in the case of
   Chorin-Temam's method, Yosida one or inexact LU precondioning.

   REMARK: H is diagonal and then stored in vector type.
 */
template <typename MatrixTypeC,
  typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
void my_precSchur_PC(double *z, int *options, int *proc_config, double *params,
		     AZ_MATRIX *Amat, AZ_PRECOND *prec)
{
  // Extraction of H, D, HinvC, HinvDtr and trD stored in
  // the structure AZ_MATRIX
  DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,
    MatrixTypeMpLp,VectorType>* my_data = static_cast< DataFactorisation<MatrixTypeC,MatrixTypeD,
    MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>* >(AZ_get_precond_data(prec));

  UInt dim  = my_data->_trD ->Patt()->nCols();
  double lastz=z[dim-1];

  // AZTEC specifications
  int    options_s[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_s[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_s[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                    // indicating success or failure.
  //

  //AZTEC matrix S=D*H^{-1}*D^T
  AZ_MATRIX *matS;

  int N_eq= dim;
  matS= AZ_matrix_create(N_eq);
  // data containing the matrices H and trD as pointers
  // are passed through S :
  AZ_set_MATFREE(matS, my_data, my_approxmatvec<
		 MatrixTypeC,
		 MatrixTypeD,
		 MatrixTypeDtr,
		 MatrixTypeH,
		 MatrixTypeMpLp,
		 VectorType>);

  my_data->_dataAztec_s->aztecOptionsFromDataFile(options_s,params_s);

  //----------------------------
  //(i) Solve the system S*z1= z
  //----------------------------

  // temporary vector z1
  Vector z1(dim);
  z1=0.;


  // keep factorisation for reusing it
  matS->data_org[AZ_name]= DATA_NAME_AZTEC_MATS;     // identification
  options_s[AZ_keep_info]= 1;                        // keep information

  //initialisation of first recursion level for AZTEC memory manager
  options_s[AZ_recursion_level]=my_data->_recur;


  AZ_iterate(&z1[0], z, options_s, params_s, status_s,
	     proc_config, matS, NULL, NULL);


  //-------------------------------
  //(ii) Solve the system S*z= B*z1
  //-------------------------------

  //compute rhs B*z1
  Vector Bz1(dim);
  Bz1= (*(my_data->_D))*( (*(my_data->_HinvC)) *
			( (*(my_data->_HinvDtr))*z1 ) );

  //recall options used for the first linear system solving (i)
  options_s[AZ_pre_calc]=AZ_reuse;
  matS->data_org[AZ_name]= DATA_NAME_AZTEC_MATS;
  options_s[AZ_keep_info]= 0;


  AZ_iterate(z, &Bz1[0], options_s, params_s, status_s,
	     proc_config, matS, NULL, NULL);

  //return to zero recursion level for AZTEC memory manager
  //--options_s[AZ_recursion_level];

  // destroy matS
  AZ_matrix_destroy(&matS);

  if (my_data->_fullEssential)
    z[dim-1]=lastz; // diagonalisation of the last row.

}

////////////////////////
/*!
   \brief Cahouet-Chabart Preconditioner of the Schur complement passed to AZTEC.

   The aim is to solve the system (nu*Mp^{-1}+S^{-1})^{-1} * z = r where
   S=(D*H^{-1}*trD), Mp = pressure mass matrix and
   in input z=r.
   the matrices Mp, H, D, HinvC, HinvDtr and trD are stored in the
   class DataFactorisation.
   In output the result is set in the vector z.
   H corresponds to an approximation of the matrix C.

   REMARK: H is diagonal and then stored in vector type.
 */
template <typename MatrixTypeC,
  typename MatrixTypeD,
  typename MatrixTypeDtr,
  typename MatrixTypeH,
  typename MatrixTypeMpLp,
  typename VectorType>
void my_precSchur_CC(double *z, int *options, int *proc_config, double *params,
		     AZ_MATRIX *Amat, AZ_PRECOND *prec)
{
  // Extraction of Mp, H, D, HinvC, HinvDtr and trD stored in
  // the structure AZ_MATRIX
  DataFactorisation<MatrixTypeC,MatrixTypeD,MatrixTypeDtr,MatrixTypeH,
  MatrixTypeMpLp,VectorType>* my_data = static_cast< DataFactorisation<MatrixTypeC,MatrixTypeD,
    MatrixTypeDtr,MatrixTypeH,MatrixTypeMpLp,VectorType>* >(AZ_get_precond_data(prec));

  UInt dim  = my_data->_trD->Patt()->nCols();
  double lastz=z[dim-1];

  // ----------------------------
  // (i) Solve the system S*z1= z
  // ----------------------------

  // AZTEC specifications
  int    options_s[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_s[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_s[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                    // indicating success or failure.
  //

  //AZTEC matrix S=D*H^{-1}*D^T
  AZ_MATRIX *matS;

  int N_eq= dim;
  matS= AZ_matrix_create(N_eq);
  // data containing the matrices H and trD as pointers
  // are passed through S :
  AZ_set_MATFREE(matS, my_data, my_approxmatvec<
                 MatrixTypeC,
		 MatrixTypeD,
                 MatrixTypeDtr,
		 MatrixTypeH,
		 MatrixTypeMpLp,
		 VectorType>);

  my_data->_dataAztec_s->aztecOptionsFromDataFile(options_s,params_s);

  // temporary vector z1
  Vector z1(dim);
  z1=0.;

  //initialisation of first recursion level for AZTEC memory manager
  options_s[AZ_recursion_level]=my_data->_recur;

  AZ_iterate(&z1[0], z, options_s, params_s, status_s,
	     proc_config, matS, NULL, NULL);

  //return to zero recursion level for AZTEC memory manager
  //--options_s[AZ_recursion_level];

  // ------------------------------
  // (ii) Solve the system Mp*z2= z
  // ------------------------------

  // AZTEC specifications
  int    proc_config_ii[AZ_PROC_SIZE];// Processor information:
  int    data_org_ii[AZ_COMM_SIZE];   // data organisation
  int    options_ii[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_ii[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_ii[AZ_STATUS_SIZE];   // Information returned from AZ_solve()

  AZ_set_proc_config(proc_config_ii, AZ_NOT_MPI);

  AZ_defaults(options_ii, params_ii);

  options_ii[AZ_output]=AZ_none;
  options_ii[AZ_precond] =  AZ_dom_decomp;
  options_ii[AZ_subdomain_solve] = AZ_ilut;
  options_ii[AZ_conv]     = AZ_rhs;
  options_ii[AZ_poly_ord] = 5;
  options_ii[AZ_kspace]   = 40;
  params_ii[AZ_tol]       = 1.00e-6;
  params_ii[AZ_drop]      = 1.00e-4;

  // data_org assigned "by hands" while no parallel computation is performed
  data_org_ii[AZ_N_internal]= N_eq;
  data_org_ii[AZ_N_border]= 0;
  data_org_ii[AZ_N_external]= 0;
  data_org_ii[AZ_N_neigh]= 0;
  data_org_ii[AZ_name]= 150;

  // AZTEC matrix for Mp
  AZ_MATRIX *Mp_az;
  // AZTEC preconditioner for Mp
  AZ_PRECOND *prec_Mp;

  Mp_az= AZ_matrix_create(N_eq);
  AZ_set_MSR(Mp_az, (int*) my_data->_Mp->Patt()->giveRaw_bindx(),
	     (double*) my_data->_Mp->giveRaw_value(),
	     data_org_ii, 0, NULL, AZ_LOCAL);

  prec_Mp= AZ_precond_create(Mp_az, AZ_precondition, NULL);

  //rhs v_muz = mu * z
  Vector v_muz(dim);
  point2Vector(z,v_muz);
  v_muz = my_data->_mu * v_muz;

  // keep factorisation for reusing it
  options_ii[AZ_keep_info]= 1;                        // keep information
  // reuse information if computed in former iteration
  if (*(my_data->_idFacto))
    options_ii[AZ_pre_calc]=AZ_reuse;

  //initialisation of first recursion level for AZTEC memory manager
  options_ii[AZ_recursion_level]=my_data->_recur;

  AZ_iterate(z, &v_muz[0], options_ii, params_ii, status_ii,
	     proc_config_ii, Mp_az, prec_Mp, NULL);

  //return to zero recursion level for AZTEC memory manager
  //--options_ii[AZ_recursion_level];

  // destroy matS
  AZ_matrix_destroy(&matS);
  AZ_matrix_destroy(&Mp_az);
  AZ_precond_destroy(&prec_Mp);

  // switch idFacto to 1 after first iteration
  *(my_data->_idFacto) = 1;

  for (UInt i=0; i< dim; i++)
    z[i]+=z1[i];

  if (my_data->_fullEssential)
    z[dim-1]=lastz; // diagonalisation of the last row.
}
}

#endif

