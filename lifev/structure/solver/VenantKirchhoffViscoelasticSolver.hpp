//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief VenantKirchhoffViscoelasticSolver - Class to solve second order problem as linear visco-elastic
 *  problem and waves problem.
 *
 *  @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 *
 *  @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 *
 *  @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */


#ifndef  VENANTKIRCHHOFFVISCOELASTICSOLVER_H_
#define  VENANTKIRCHHOFFVISCOELASTICSOLVER_H_ 1


#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#include<boost/scoped_ptr.hpp>


#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>

#include <lifev/structure/solver/VenantKirchhoffViscoelasticData.hpp>
#include <lifev/core/util/Displayer.hpp>


namespace LifeV
{
//! SecondOrderProblem this class solver second order problem,  waves equation and linear viscoelastic problem
/*!
  \class VenantKirchhoffViscoelasticSolver
  \brief
    This class solves the following waves equation:

 \f[M \frac{d^2 u}{d t^2} + D(u,\frac{d u}{d t} ) \frac{d u}{d t} + A(u, \frac{d u}{dt}) = f\f]

where \f$M\f$ is mass matrix, \f$A\f$  stiffness matrix and \f$D\f$ is  damping matrix given by  \f$ \beta M + \gamma A\f$.

 If \f$A\f$ and \f$D\f$ depend on \f$u\f$ and \f$dot{u}\f$ we linearize the problem using suitable extrapolations.
 we use the time advancing scheme  defined in TimeAdvanceBase class for details see TimeAdvanceBase.pdf notes.
*/

template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >

class VenantKirchhoffViscoelasticSolver
{
public:

    //! @name Public Types
    //@{

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > source_type;

    typedef BCHandler                                                        bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type>  bchandler_type;

    typedef SolverType                                               solver_type;

    typedef typename solver_type::matrix_type     matrix_type;
    typedef boost::shared_ptr<matrix_type>         matrix_ptrtype;
    typedef typename solver_type::vector_type     vector_type;
    typedef boost::shared_ptr<vector_type>          vector_ptrtype;

    typedef typename SolverType::prec_raw_type prec_raw_type;
    typedef typename SolverType::prec_type         prec_type;

    typedef VenantKirchhoffViscoelasticData                                    data_type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    VenantKirchhoffViscoelasticSolver( );

    //! virtual Destructor
    virtual ~VenantKirchhoffViscoelasticSolver() {};

    //@}

    //! @name Methods
    //@{

    //! class setup
    /*!
    @param data file
    @param feSpace finite element space
    @param comm communicator
    */
    void setup ( boost::shared_ptr<data_type> data,
                 const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   feSpace,
                 boost::shared_ptr<Epetra_Comm>&     comm
               );

    //!  class setup
    /*!
    @param data file
    @param feSpace finite element space
    @param BCh boundary condition
    @param comm communicator
    */
    void setup ( boost::shared_ptr<data_type> data,
                 const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   feSpace,
                 bchandler_type&       BCh,
                 boost::shared_ptr<Epetra_Comm>&     comm
               );


    //! class setup
    /*!
    @param data file
    @param feSpace finite element space
    @param comm communicator
    @param epetraMap  epetra vector
    @param offset
    */
    virtual void setup ( boost::shared_ptr<data_type> data,
                         const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   feSpace,
                         boost::shared_ptr<Epetra_Comm>&     comm,
                         const boost::shared_ptr<const MapEpetra>&      epetraMap,
                         UInt       offset = 0   );


    //! buildSystem
    /*!
     @param xi  is the coefficient of mass term defined as \f$\xi^0/dt^2\f$;
     @param alpha is the coefficient of damping term defined as \f$\alpha^0/dt\f$;
     @note by default \f$xi\f$ is equal to 1 and \f$alpha\f$ is equal to 0;
     */
    void buildSystem (const Real& xi = 1, const  Real& alpha = 0);

    //! buildSystem
    /*!
    @param matrix is the Mass \f$ + \f$  Damping \f$+\f$ Stiffness
    @param xi is the coefficient of mass term defined as \f$\xi^0/dt^2\f$;
    */
    void buildSystem (matrix_ptrtype matrix, const Real& xi);

    //! build Damping matrix
    /*!
    building damping matrix defined as \f$\gamma  A + \beta  M\f$
    @param damping matrix
    @param alpha is the coefficient of time advancing
    scheme defined as \f$\alpha^0/dt\f$
    @notes alpha can be 1 and we can introduce
    the time advancing coefficient in the updateSystem
    */
    void buildDamping (matrix_ptrtype damping, const Real& alpha);

    //!updateSystem with coefficients of time advance methods
    /*
    @param xi parameter of mass matrix \f$(\xi^0/dt^2)\f$;
    @param alpha parameter of damping matrix \f$(\alpha^0/dt)\f$;
    @param rhs  vector
    @note: this method must be used when we call buildSystem(1,1);
    in this case we must set the coefficients of time advancing scheme;
    */
    void updateSystem (const vector_type& rhs,
                       const Real& xi = 0,
                       const Real& alpha = 0 );

    //! compute the start solution
    /*! we consider the system \f$M w + D v + A u = f\f$;
    in general we knowns the \f$u^0\f$, \f$v^0\f$ but not  \f$w^0\f$.
    This method computes \f$w^0\f$  starting by \f$F = f^0 - D v^0 -A u^0\f$
    and boundary conditions associated to \f$w\f$;
    @param solution is the vector where save the solution;
    @param bch the boundary condition respect to \f$w\f$;
    @param rhs is \f$f^0 - D v^0 -A u^0\f$
    */
    void computeStartValue ( vector_type& solution, const  bchandler_raw_type& bch, const vector_type& rhs );

    //! Solve the system
    /*!
    solve the system
    @param  bch is boundary condition;
    */
    void iterate ( bchandler_raw_type& bch );

    //! update source term
    /*!
    @param source is the vector where we evaluate the source term;
    */
    void updateSourceTerm (const vector_type&  source);

    //! updateRHS
    /*!
    @param rhs is the right and side;
    */
    inline void updateRHS (const vector_type& rhs)
    {
        M_displayer->leaderPrint ("  P-  Updating right hand side... ");

        LifeChrono chrono;
        chrono.start();

        *M_rhsNoBC = rhs;
        M_rhsNoBC->globalAssemble();

        chrono.stop();
        M_displayer->leaderPrintMax ("done in ", chrono.diff() );
    }

    //! reset the Prec
    void resetPrec()
    {
        M_resetPrec = true;
    }

    //@}

    //@name Set Methods
    //@{

    //! Sets source term
    /*!
    @param source is function
    */
    void setSourceTerm ( source_type const& source )
    {
        M_source = source;
    }

    //! Set parameters
    /*!
    @param dataFile is the file contains the parameters of problem
    this method must be removed
    */
    void setDataFromGetPot ( const GetPot& dataFile );

    //@}

    //@name Get Methods
    //@{

    //! Return the  map
    /*!
     @returns epetraMap
     */
    MapEpetra   const& map()       const
    {
        return M_localMap;
    }

    //! Return the  map
    /*!
     @returns  displayer
     */
    Displayer   const& displayer() const
    {
        return M_displayer;
    }

    //!Return the mass matrix
    /*!
     @returns the mass matrix (it doesn't divided by \f$\xi/dt^2\f$)
     */
    matrix_type const matrMass()   const
    {
        return *M_matrMass;
    }

    //!Return the damping matrix
    /*!
    @returns the damping matrix (it doesn't divided by \f$\alpha/dt\f$)
    */
    matrix_ptrtype const matrDamping()  const
    {
        return M_matrDamping;
    }

    //!Return the damping matrix
    /*!
    @returns the system matrix given by Stiffness + Damping \f$\alpha/dt\f$ + Mass \f$\xi/dt^2\f$
    */
    matrix_ptrtype const matrSystem()   const
    {
        return M_matrSystem;
    }

    //!Return the Stiffness matrix
    /*!
    @returns the stiffness matrix
    */
    matrix_ptrtype const matrLinearStiff() const
    {
        return M_matrLinearStiffness;
    }

    //! Return FESpace
    FESpace<Mesh, MapEpetra>& feSpace()
    {
        return M_FESpace;
    }

    //! BCHandler getter and setter
    /*!
     @returns the boundary conditions
     */
    BCHandler const& BChandler() const
    {
        return M_BCh;
    }

    //! Return the solution
    /*!
     @returns the solution
     */
    vector_ptrtype&  solution()
    {
        return M_solution;
    }

    //! Return residual
    /*
     @returns the residual
     */
    vector_ptrtype& residual()
    {
        return M_residual;
    }

    //!Return right hand side without boundary conditions
    /*
    @returns the right hand side without boundary conditions;
    */
    vector_ptrtype& rhsContributionSecondDerivativeithoutBC()
    {
        return M_rhsNoBC;
    }

    //!return the communicator
    /*!
    @returns the communicator
    */
    boost::shared_ptr<Epetra_Comm> const& comm()   const
    {
        return M_displayer->comm();
    }

    //! Return offset;
    /*
    @returns offset;
    */
    const UInt& offset() const
    {
        return M_offset;
    }

    //!thickness
    /*!
    @returns the thickness
    */
    inline const Real& thickness() const
    {
        return M_data->thickness();
    }
    //!density
    /*!
    @returns the density
     */
    inline const Real& density()   const
    {
        return M_data->rho();
    }

    //!young
    /*!
    @returns the young
    */
    inline const Real& young()     const
    {
        return M_data->young();
    }
    //!poisson
    /*!
    @returns the poisson
    */
    inline const Real& poisson()   const
    {
        return M_data->poisson();
    }

    //@}

private:
    //@ Private Methods

    //! apply boundary condition
    /*!
    @param matrix of system;
    @param rhs right hand side without boundaryCondition;
    @param BCh boundary condition;
    @param offset
     */
    void applyBoundaryConditions (matrix_type& matrix,
                                  vector_type& rhs,
                                  bchandler_raw_type& BCh,
                                  UInt         offset = 0);

    //@name Attributes

protected :

    //! data of problem
    boost::shared_ptr<data_type>   M_data;

    //! feSpace
    boost::shared_ptr<FESpace<Mesh, MapEpetra> >      M_FESpace;

    //! displayer
    boost::scoped_ptr<Displayer>   M_displayer;

    //! current rank
    Int                            M_me;

    //! boundary condition
    bchandler_type   M_BCh;

    //! Epetra map need to define the VectorEpetra;
    boost::shared_ptr<const MapEpetra>       M_localMap;

    //! Matrix  mass
    matrix_ptrtype                          M_matrMass;

    //! Matrix  stiffness linear
    matrix_ptrtype               M_matrLinearStiffness;

    //! Matrix System
    matrix_ptrtype                 M_matrSystem;

    //! Matrix Damping \f$\gamma  A + \beta  M\f$
    matrix_ptrtype                 M_matrDamping;

    //!@name Elementary matrices and vectors:
    //@{
    //! linear stiffness
    boost::shared_ptr<MatrixElemental>                       M_elmatK;

    //! mass
    boost::shared_ptr<MatrixElemental>                       M_elmatM;
    //!mass+ linear stiffness
    boost::shared_ptr<MatrixElemental>                       M_elmatC;

    //!damping
    boost::shared_ptr<MatrixElemental>                       M_elmatD;
    //@}

    //! unknowns vector
    vector_ptrtype                    M_solution;

    //! right hand  side
    vector_ptrtype                 M_rhs;

    //! right  hand  side without boundary condition
    vector_ptrtype                 M_rhsNoBC;

    //! residual
    boost::shared_ptr<vector_type> M_residual;

    //! source term
    source_type                M_source;

    //! data for solving tangent problem with aztec
    boost::shared_ptr<solver_type> M_linearSolver;

    //! true if reuse the preconditonar
    bool                            M_reusePrec;

    //! max number of iteration before to update preconditonator
    UInt                            M_maxIterForReuse;

    //! reset preconditionator
    bool                            M_resetPrec;

    //! maximum number of iteration for solver method
    UInt                            M_maxIterSolver;

    //! offset
    UInt                            M_offset;

};

// ==============================================================
// Constructors & Destructors
// ==============================================================

//! Empty Constructor

template <typename Mesh, typename SolverType>
VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
VenantKirchhoffViscoelasticSolver( ) :
    M_data                             (    ),
    M_FESpace                          (    ),
    M_displayer                        (    ),
    M_me                               ( 0  ),
    M_BCh                              (    ),
    M_localMap                         (    ),
    M_matrMass                         (    ),
    M_matrLinearStiffness              (    ),
    M_matrSystem                       (    ),
    M_matrDamping                      (    ),
    M_elmatK                           (    ),
    M_elmatM                           (    ),
    M_elmatC                           (    ),
    M_elmatD                           (    ),
    M_solution                         (    ),
    M_rhs                              (    ),
    M_rhsNoBC                          (    ),
    M_residual                         (    ),
    M_source                           (    ),
    M_linearSolver                     (    ),
    M_reusePrec                        (    ),
    M_maxIterForReuse                  (    ),
    M_resetPrec                        (    ),
    M_maxIterSolver                    (    ),
    M_offset                           (    )

{
}

// ==============================================================
// Methods
// ==============================================================

template <typename Mesh, typename SolverType>
void
VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::setup (
    boost::shared_ptr<data_type>        data,
    const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& feSpace,
    boost::shared_ptr<Epetra_Comm>&     comm,
    const boost::shared_ptr<const MapEpetra>&  epetraMap,
    UInt                                offset
)
{
    M_data                                        = data;
    M_FESpace                                 = feSpace;
    M_displayer.reset                      (new Displayer (comm) );
    M_me                                           = comm->MyPID();
    M_elmatK.reset                           ( new MatrixElemental ( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
    M_elmatM.reset                           ( new MatrixElemental ( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
    M_elmatC.reset                           ( new MatrixElemental ( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
    M_elmatD.reset                           ( new MatrixElemental ( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
    M_localMap                                   = epetraMap;
    M_solution.reset                          (new vector_type (*M_localMap) );
    M_rhsNoBC.reset                         ( new vector_type (*M_localMap) );
    M_matrMass.reset                       (new matrix_type (*M_localMap) );
    M_matrLinearStiffness.reset     (new matrix_type (*M_localMap) );
    M_matrDamping.reset               (new matrix_type (*M_localMap) );
    M_matrSystem.reset                  (new matrix_type (*M_localMap) );
    M_offset                                        = offset;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::setup (
    boost::shared_ptr<data_type>        data,
    const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& feSpace,
    boost::shared_ptr<Epetra_Comm>&     comm
)
{
    setup ( data, feSpace, comm, feSpace->mapPtr(), (UInt) 0 );
    M_rhs.reset                              ( new vector_type (*M_localMap) );
    M_residual.reset                     ( new vector_type (*M_localMap) );
    M_linearSolver.reset              ( new SolverType ( comm ) );

}

template <typename Mesh, typename SolverType>
void
VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::setup (
    boost::shared_ptr<data_type>          data,
    const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& feSpace,
    bchandler_type&                BCh,
    boost::shared_ptr<Epetra_Comm>&              comm
)
{
    setup (data, feSpace, comm);
    M_BCh = BCh;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::setDataFromGetPot ( const GetPot& dataFile )
{
    M_linearSolver->setDataFromGetPot ( dataFile, "problem/solver" );
    M_linearSolver->setupPreconditioner (dataFile, "problem/prec");
}

template <typename Mesh, typename SolverType>
void
VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
buildSystem (const Real& xi, const Real& alpha)
{
    M_displayer->leaderPrint ("  P-  Computing constant matrices ...          ");
    LifeChrono chrono;
    chrono.start();

    // these lines must be removed next week
    this->buildSystem (M_matrSystem, xi);

    if (M_data->damping() )
    {
        M_matrDamping.reset (new matrix_type (*M_localMap) );
        this->buildDamping (M_matrSystem, alpha);
    }
    M_matrSystem->globalAssemble();

    chrono.stop();
    M_displayer->leaderPrintMax ( "done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void
VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
buildSystem (matrix_ptrtype matrSystem, const Real& xi)
{
    M_displayer->leaderPrint ( "P-  Building the system             ... ");

    LifeChrono chrono;
    chrono.start();

    // Elementary computation and matrix assembling
    // Loop on elements

    for ( UInt iVol = 0; iVol < this->M_FESpace->mesh()->numVolumes(); ++iVol )
    {
        this->M_FESpace->fe().updateFirstDeriv ( this->M_FESpace->mesh()->element ( iVol ) );

        Int marker    = this->M_FESpace->mesh()->volumeList ( iVol ).markerID();

        this->M_elmatK->zero();
        this->M_elmatM->zero();

        // building stiffness matrix
        stiff_strain (   2 * M_data->mu (marker), *this->M_elmatK, this->M_FESpace->fe() );
        stiff_div   ( M_data->lambda (marker), *this->M_elmatK, M_FESpace->fe() );

        this->M_elmatC->mat() = this->M_elmatK->mat();

        // mass*xi to compute mass+stiff

        mass ( xi * M_data->rho(), *this->M_elmatM, this->M_FESpace->fe(), 0, 0, nDimensions );

        this->M_elmatC->mat() += this->M_elmatM->mat();

        //mass
        this->M_elmatM->zero();
        mass (M_data->rho(), *this->M_elmatM, this->M_FESpace->fe(), 0, 0, nDimensions );

        assembleMatrix ( *M_matrLinearStiffness,
                         *this->M_elmatC,
                         this->M_FESpace->fe(),
                         this->M_FESpace->fe(),
                         this->M_FESpace->dof(),
                         this->M_FESpace->dof(),
                         0, 0, 0, 0);

        assembleMatrix ( *matrSystem,
                         *this->M_elmatC,
                         this->M_FESpace->fe(),
                         this->M_FESpace->fe(),
                         this->M_FESpace->dof(),
                         this->M_FESpace->dof(),
                         0, 0, 0, 0);

        assembleMatrix ( *this->M_matrMass,
                         *this->M_elmatM,
                         this->M_FESpace->fe(),
                         this->M_FESpace->fe(),
                         this->M_FESpace->dof(),
                         this->M_FESpace->dof(),
                         0, 0, 0, 0);
    }

    M_matrMass->globalAssemble();
    M_matrLinearStiffness->globalAssemble();
    chrono.stop();
    M_displayer->leaderPrintMax ("done in ", chrono.diff() );

}

template <typename Mesh, typename SolverType>
void VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
buildDamping (matrix_ptrtype damping, const Real& alpha)
{
    LifeChrono chrono;
    chrono.start();

    M_displayer->leaderPrint ( "P-  Building the system   Damping matrix          ... ");

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt iVol = 0; iVol < this->M_FESpace->mesh()->numVolumes(); ++iVol )
    {
        this->M_FESpace->fe().updateFirstDeriv ( this->M_FESpace->mesh()->element ( iVol ) );

        Int marker    = this->M_FESpace->mesh()->volumeList ( iVol ).markerID();

        Real gamma = M_data->gamma (marker);
        Real beta  = M_data->beta (marker);

        //building damping matrix

        this->M_elmatD->zero();

        stiff_strain ( 2.0 * M_data->mu (marker) * gamma , *this->M_elmatD, this->M_FESpace->fe() );
        stiff_div   ( M_data->lambda (marker) * gamma, *this->M_elmatD, this->M_FESpace->fe() );
        mass ( beta * M_data->rho(), *this->M_elmatD, this->M_FESpace->fe(), 0, 0);

        assembleMatrix ( *damping,
                         *this->M_elmatD,
                         this->M_FESpace->fe(),
                         this->M_FESpace->fe(),
                         this->M_FESpace->dof(),
                         this->M_FESpace->dof(),
                         0, 0, 0, 0);
    }

    *M_matrSystem += *damping * alpha;
    *M_matrDamping = *damping;

    M_matrDamping->globalAssemble();
    M_displayer->leaderPrintMax ( " done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
updateSystem (const vector_type& rhs,
              const Real&        xi,
              const Real&        alpha)
{

    LifeChrono chrono;

    updateRHS (rhs);

    if (xi == 0 & alpha == 0 )
    {
        M_displayer->leaderPrint ("P - use the same System matrix ...  \n ");
    }
    else
    {
        M_displayer->leaderPrint ("P - updating the System matrix ....      ");
        chrono.start();
        M_matrSystem = M_matrLinearStiffness;

        if (xi != 0. )
        {
            *M_matrSystem += *M_matrMass * xi;
        }

        if (alpha != 0)
        {
            *M_matrSystem += *M_matrDamping * alpha;
        }

        M_matrSystem->GlobalAssemble();

        chrono.stop();
        M_displayer->leaderPrintMax ("done in ", chrono.diff() );
    }

}

template <typename Mesh, typename SolverType>
void VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
computeStartValue ( vector_type& solution, const  bchandler_raw_type& bch, const vector_type& rhs )
{
    vector_type rhsFull (rhs);

    prec_type prec;

    std::string precType = "Ifpack" ;

    prec.reset ( PRECFactory::instance().createObject (precType) );

    prec->buildPreconditioner (M_matrMass);

    Real condest = prec->Condest();

    M_linearSolver->setPreconditioner (prec);

    M_linearSolver->setMatrix (*M_matrMass);

    Int numIter =  M_linearSolver->solve (solution, rhs);
}

template <typename Mesh, typename SolverType>
void VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
updateSourceTerm (const  vector_type&  source)
{
    M_displayer->leaderPrint ("P - updating the Source Term....      ");
    LifeChrono chrono;
    chrono.start();

    for ( UInt iVol = 0; iVol < this->M_FESpace->mesh()->numVolumes(); ++iVol )
    {
        //  M_elvecSource.zero();
        Real f, x, y, z;

        UInt i, inod, ig, ic;
        UInt eleID = this->M_FESpace->fe().currentLocalId();
        Real u_ig;
        for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
        {
            for ( ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ++ig )
            {
                this->M_FESpace->fe().coorQuadPt ( x, y, z, ig );
                f = M_source (M_data->dataTime()->time(), x, y, z, iComp + 1 );
                u_ig = 0.;

                for ( i = 0; i < M_FESpace->fe().nbFEDof(); ++i )
                {
                    inod = this->M_FESpace->dof().localToGlobalMap ( eleID, i ) + iComp * this->M_FESpace->dof().numTotalDof();
                    u_ig = f * this->M_FESpace->fe().phi ( i, ig );
                    source.sumIntoGlobalValues (inod, u_ig * this->M_FESpace->fe().weightDet ( ig ) );
                }
            }
        }
    }
    source.GlobalAssemble();

    chrono.stop();
    M_displayer->leaderPrintMax ("done on ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
iterate ( bchandler_raw_type& bch )
{
    LifeChrono chrono;

    // matrix and vector assembling communication
    M_displayer->leaderPrint ("  P-  Solving the system ... \n");

    chrono.start();

    matrix_ptrtype matrFull ( new matrix_type (*M_localMap, M_matrSystem->meanNumEntries() ) );
    *matrFull += *M_matrSystem;

    M_rhsNoBC->globalAssemble();

    vector_type rhsFull (*M_rhsNoBC);

    // boundary conditions update
    M_displayer->leaderPrint ("  P-  Applying boundary conditions ...         ");

    chrono.start();
    this->applyBoundaryConditions (*matrFull, rhsFull, bch);

    chrono.stop();

    M_displayer->leaderPrintMax ("done in " , chrono.diff() );

    M_linearSolver->resetPreconditioner();
    M_resetPrec = false;

    // solving the system
    M_linearSolver->setMatrix (*matrFull);
    Real numIter = M_linearSolver->solveSystem ( rhsFull, *M_solution, matrFull);

    numIter = std::abs (numIter);

    if (numIter >= M_maxIterForReuse || numIter >= M_maxIterSolver)
    {
        resetPrec();
    }

    *M_residual =  *M_matrSystem * ( *M_solution);
    *M_residual -= *M_rhsNoBC;

} // iterate()

//=========================================================
// Private Methods
//=========================================================

template<typename Mesh, typename SolverType>
void VenantKirchhoffViscoelasticSolver<Mesh, SolverType>::
applyBoundaryConditions (matrix_type&        matrix,
                         vector_type&        rhs,
                         bchandler_raw_type& BCh,
                         UInt                offset)
{
    if (offset)
    {
        BCh.setOffset (offset);
    }

    if ( !BCh.bcUpdateDone() )
    {
        BCh.bcUpdate ( *this->M_FESpace->mesh(), this->M_FESpace->feBd(), this->M_FESpace->dof() );
    }

    vector_type rhsFull (rhs, Unique); // bcManages now manages the also repeated parts

    bcManage ( matrix, rhsFull, *this->M_FESpace->mesh(), this->M_FESpace->dof(), BCh, this->M_FESpace->feBd(), 1., M_data->dataTime()->time() );

    // matrix should be GlobalAssembled by  bcManage

    rhs = rhsFull;

} // applyBoundaryCondition

} // namespace
#endif /* VENANTKIRCHHOFFVISCOELASTICSOLVER */
