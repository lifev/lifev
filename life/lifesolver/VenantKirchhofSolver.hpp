//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief This file contains solvers for St. Venant-Kirchhof materials (linear for the moment)
 *
 *  @version 1.0
 *  @date 01-06-2003
 *  @author Miguel Angel Fernandez
 *
 *  @version 1.1
 *  @date 01-09-2010
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _VENANTKIRCHHOFSOLVER_H_
#define _VENANTKIRCHHOFSOLVER_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifecore/chrono.hpp>

#include <life/lifealg/dataNewton.hpp>
#include <life/lifealg/newton.hpp>
#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifesolver/dataElasticStructure.hpp>
#include <life/lifecore/displayer.hpp>

#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#include<boost/scoped_ptr.hpp>

namespace LifeV
{


/*!
  \class VenantKirchhofSolver
  \brief
  This class solves the linear elastodynamics equations for a (only linear right now)
  St. Venant-Kirchoff material

*/
template <typename Mesh,
typename SolverType = LifeV::SolverTrilinos >

class VenantKirchhofSolver
{
public:

    //!@name Type definitions
    //@{
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_type;

    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;


    typedef SolverType                    solver_type;

    typedef typename solver_type::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>         matrix_ptrtype;
    typedef typename solver_type::vector_type      vector_type;
    typedef boost::shared_ptr<vector_type>         vector_ptrtype;


    typedef typename SolverType::prec_raw_type     prec_raw_type;
    typedef typename SolverType::prec_type         prec_type;

    typedef DataElasticStructure                   data_type;

    typedef singleton<factory<VenantKirchhofSolver,  std::string> >  StructureSolverFactory;

    //@}


    //! @name Constructor
    //@{

    VenantKirchhofSolver();

    //@}

    //!@name Methods
    //@{

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
      \param BCh boundary conditions for the displacement
    */
    void setup( boost::shared_ptr<data_type> data,
                const boost::shared_ptr< FESpace<Mesh, EpetraMap> >&   FESpace,
                bchandler_type&       BCh,
                boost::shared_ptr<Epetra_Comm>&     comm
              );

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
    */

    void setup( boost::shared_ptr<data_type> data,
                const boost::shared_ptr< FESpace<Mesh, EpetraMap> >&   FESpace,
                boost::shared_ptr<Epetra_Comm>&     comm
              );

    virtual void setup( boost::shared_ptr<data_type> data,
                        const boost::shared_ptr< FESpace<Mesh, EpetraMap> >&   dFESpace,
                        boost::shared_ptr<Epetra_Comm>&     comm,
                        const boost::shared_ptr<const EpetraMap>&       monolithicMap,
                        UInt       offset=0
                                          //boost::shared_ptr<FESpace<Mesh, EpetraMap> >   uFESpace=0
                      );


    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    //    void updateSystem( source_type const& source );
    void updateSystem();

    virtual void updateSystem(matrix_ptrtype& stiff);

    //void buildSystem(matrix_type & bigMatrixStokes); // used for monolithic
    virtual void buildSystem(matrix_ptrtype matrix, const Real& factor=1.);
    void buildSystem();

    //! Solve the non-linear system

    void iterate( vector_type& sol );
    virtual void iterate( bchandler_type& bch );
    void iterateLin( bchandler_type& bch );

    //! Output
    //void showMe( std::ostream& c = std::cout ) const;

    virtual void updateJacobian( vector_type& sol, matrix_ptrtype& jac )=0;
  
    //! solves the tangent problem for newton iterations
    virtual void solveJac( vector_type&       step,
                           const vector_type& res,
                           Real&            linear_rel_tol)=0 ;
    //    void solveJac( const Vector& res, Real& linear_rel_tol, Vector &step);
    //! solves the tangent problem with custom BC

    virtual void solveJacobian( vector_type&       step,
                                const vector_type& res,
                                Real&            linear_rel_tol,
                                bchandler_type&    BCd )=0 ;


    //! evaluates residual for newton interations
    void evalResidual( vector_type &res, const vector_type& sol, int iter);

    void evalConstraintTensor();

    virtual void initialize( const Function& d0, const Function& w0, const Function& a0 = Function() );
    void initialize( vector_ptrtype d0,  vector_ptrtype w0 = vector_ptrtype(),  vector_ptrtype a0 = vector_ptrtype() );
    void initializeVel( const vector_type& w0);

    virtual void updateVel();

    void reduceSolution( Vector& disp, Vector& vel );

    void rescaleMatrices(); // used for monolithic

    /**
       in the linear case the solid matrix is constant, thus it does not need to be recomputed.
     */

    void computeMatrix( matrix_ptrtype& stiff, const vector_type& sol, Real const& factor )
    {
    }

    void computeMatrix( const vector_type& sol, Real const& factor )
    {
    }
  
    //void updateMatrix(matrix_type & bigMatrixStokes);// used for monolithic
    //void updateCoupling(matrix_type couplingMatrix);// used for monolithic

    //@}

    //! @name Set Methods
    //@{
    //!setters
    void setBC(bchandler_type& BCd)   {M_BCh = BCd;}
    void setSourceTerm( source_type const& __s ) { M_source = __s; }

    void resetPrec(bool reset = true) { if (reset) M_linearSolver.precReset(); }

    virtual void setDisp(const vector_type& disp) {*M_disp = disp;} // used for monolithic
    //! recur setter
    void setRecur(UInt recur) {_recur = recur;}

     void setDataFromGetPot( const GetPot& dataFile );

    //@}


    //! @name Get Methods
    //@{

    //! getters
    EpetraMap   const& getMap()       const { return *M_localMap; }

    Displayer   const& getDisplayer() const { return *M_Displayer; }

    matrix_ptrtype const getMassStiff() const {return M_massStiff; }

    matrix_ptrtype const getMass() const {return M_mass; }

    matrix_ptrtype const getLinearStiff() const {return M_linearStiff; }

    //! BCHandler getter and setter
//    LIFEV_DEPRECATED BCHandler const & BC_solid() const {return BCh_solid();}

    FESpace<Mesh, EpetraMap>& dFESpace() {return M_FESpace;}

    bchandler_type const & BChandler() const {return M_BCh;}

    //! residual getter
    vector_type& residual()             {return *M_residual_d;}

    source_type const& sourceTerm() const { return M_source; }

    vector_type& disp()        { return *M_disp; }
    vector_type& vel()         { return *M_vel; }
    vector_ptrtype& rhsWithoutBC() { return M_rhsNoBC; }

    //const Dof& dDof() const { return M_FESpace.dof(); }

    //const Mesh& mesh() const { return M_FESpace.mesh(); }

    //Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getRepeatedEpetra_Map(); }

    boost::shared_ptr<Epetra_Comm> const& comm()         const {return M_Displayer->comm();}

    Real rescaleFactor() {return M_rescaleFactor;}

    const UInt& offset() const { return M_offset; }

    /**
       Do nothing in the linear case: the matrix remains constant. Otherwise substitute the matrix with an updated one
     */
    void getSolidMatrix( matrix_ptrtype& matrix)
    {
    }

    // Physic constant
    const Real& thickness() const { return M_data->thickness(); }
    const Real& density()   const { return M_data->rho(); }
    const Real& young()     const { return M_data->young(); }
    const Real& poisson()   const { return M_data->poisson(); }
    const Real& rho()       const { return M_data->rho(); }

    //@}

protected:

    virtual void applyBoundaryConditions(matrix_type &matrix,
                                         vector_type &rhs,
                                         bchandler_type& BCh,
                                         UInt         offset=0);

    UInt dim() const { return M_FESpace->dim(); }


    //!Protected Members

    boost::shared_ptr<data_type>   M_data;

    boost::shared_ptr<FESpace<Mesh, EpetraMap> >      M_FESpace;

    boost::scoped_ptr<Displayer>   M_Displayer;

    int                            M_me;

    //! data for solving tangent problem with aztec
    boost::shared_ptr<solver_type>                    M_linearSolver;

    //! Elementary matrices and vectors
    boost::shared_ptr<ElemMat>                        M_elmatK; // stiffnes
    boost::shared_ptr<ElemMat>                        M_elmatM; // mass
    boost::shared_ptr<ElemMat>                        M_elmatC; // mass + stiffness
    //    ElemVec                        M_elvec;  // Elementary right hand side
    //    ElemVec                        M_dk_loc; // Local displacement

    //! linearized velocity

    vector_ptrtype                    M_disp;
    vector_ptrtype                    M_vel;

    //! right  hand  side displacement
    vector_ptrtype                    M_rhs;

    //! right  hand  side velocity
    vector_ptrtype                    M_rhsW;

    //! right  hand  side
    vector_ptrtype                    M_rhsNoBC;

    //! right  hand  side
    boost::shared_ptr<vector_type>                    M_f;

    //! residual
    boost::shared_ptr<vector_type>                    M_residual_d;

//    vector_type*                   M_sxx;

    vector_ptrtype                    M_sxx;
    vector_ptrtype                    M_syy;
    vector_ptrtype                    M_szz;

    //! files for lists of iterations and residuals per timestep
    std::ofstream                  M_out_iter;
    std::ofstream                  M_out_res;

    bchandler_type   M_BCh;

    boost::shared_ptr<const EpetraMap>                      M_localMap;


    //! Matrix M: mass
    matrix_ptrtype                    M_mass;

    //! Matrix C: mass + linear stiffness
    matrix_ptrtype                    M_massStiff;

    //! Matrix Knl: stiffness non-linear
    matrix_ptrtype                    M_stiff;

    //! Matrix Kl: stiffness linear
    matrix_ptrtype                    M_linearStiff;


    //! Matrix J: jacobian
    matrix_ptrtype                    M_jacobian;

    //! level of recursion for Aztec (has a sens with FSI coupling)
    UInt _recur;

    source_type                    M_source;


    int                            M_count;

    UInt                           M_offset;
    Real                            M_rescaleFactor;

    //
    //! methods
    //

    matrix_ptrtype                  M_matrFull;

};

}

#endif
