/* -*- mode: c++ -*-
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
/**
 * @file monolithic.hpp
 * @author Paolo Crosetto
 * @date 13-09-2008
 * Class handling the monolithic solver for FSI problems. The block structure of the matrix can be
 \f$\left(\begin{array}{cc}
 C&B\\
 D&N
 \end{array}\right)\f$
 if the time discretization at hand is the Geometry-Explicit one (implemented in monolithicGE.hpp), or
 \f$\left(\begin{array}{ccc}
 C&B&D\\
 D&N&0\\
 0&I&H
 \end{array}\right)\f$
 if the time discretization at hand is the Geometry-Implicit one (implemented in monolithicGI.hpp)
 */
#ifndef _MONOLITHIC_HPP
#define _MONOLITHIC_HPP

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <EpetraExt_MatrixMatrix.h>
//#include <EpetraExt_Reindex_MultiVector.h>
//#include <EpetraExt_Reindex_CrsMatrix.h>
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <life/lifecore/chrono.hpp>

#include <life/lifefem/FESpace.hpp>

#include <life/lifesolver/FSIOperator.hpp>

#include <lifemc/lifealg/ComposedPreconditioner.hpp>
#include <lifemc/lifealg/ComposedOperator.hpp>
#ifdef HAVE_TRILINOS_ANASAZI
#include <lifemc/lifealg/eigenSolver.hpp>
#endif

#include <lifemc/lifesolver/BlockMatrix.hpp>




//#include <Epetra_IntVector.h>

namespace LifeV
{
//! Monolithic.hpp pure virtual class containing the core methods of the monolithic FSI solver
/*!
 * Class handling the monolithic solver for FSI problems. The block structure of the matrix can be
 \f$\left(\begin{array}{cc}
 C&B\\
 D&N
 \end{array}\right)\f$
 if the time discretization at hand is the Geometry-Explicit one (implemented in monolithicGE.hpp), or
 \f$\left(\begin{array}{ccc}
 C&B&D\\
 D&N&0\\
 0&I&H
 \end{array}\right)\f$
 if the time discretization at hand is the Geometry-Implicit one (implemented in monolithicGI.hpp),
 * where \f$N\f$ represents the solid block, \f$C\f$ the fluid block, \f$H\f$ the harmonic extension block, while the extra
 * diagonal blocks represent the coupling. The implementation of the stress continuity coupling condition
 * is obtained by means of an augmented formulation.
 * Different possible preconditioners are implemented.
 * The flag semiImplicit in the data file is used to distinguish between the GCE and CE (with quasi Newton) time discretizations.
 * Exact Newton method and full implicit time discretization are implemented in the MonolithicGI class.
 */

class WRONG_PREC_EXCEPTION;

class Monolithic : public FSIOperator
{
public:

    //!@name Typedefs
    //@{
    typedef FSIOperator                                        super;
    typedef FSIOperator::fluid_type::value_type::matrix_type/*matrix_Type*/   matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                     matrixPtr_Type;
    typedef BlockInterface                                     prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>                   prec_ptrtype;

    typedef BlockMatrix                                        block_matrix_raw_type;
    typedef boost::shared_ptr<block_matrix_raw_type>           block_matrix_ptrtype;
    typedef SolverTrilinos                                     solver_type;
    //@}

    //! OBSOLETE typedefs
    // typedef BlockInterface                                     prec_raw_type;
//     typedef boost::shared_ptr<prec_raw_type>                   prec_ptrtype;

//     typedef BlockMatrix                                        block_matrix_raw_type;
//     typedef boost::shared_ptr<block_matrix_raw_type>           block_matrix_ptrtype;
//     typedef SolverTrilinos                                     solver_type;
    // END of OBSOLETE typedefs

    // constructors

//!@name Constructors, Destructor
//!@{
    Monolithic();
    ~Monolithic();
//!@}

//!@name Public Setup Methods
//!@{

    /**
       create FEspace
    */
    void setupFEspace();


    /**
       sets the interface map between the fluid and solid meshes
       (non scalable, do not use for massively parallel simulations)
    */
    virtual void setupDOF( void );

    /**
       reads the interface map between the fluid and solid meshes from file.
    */
    void setupDOF( mesh_filtertype& filterMesh );

    //! sets the parameters from data file
    /**
       Calls the setup of the fluid problem and the setUp method.
     */
    virtual void setupSystem( );

//     /** stores the data file into a member */
//     virtual void setDataFile( const GetPot& dataFile );

//! setup method for the Monolithic solver
/**
   sets some parameters specific to the Monolithic class
*/
    virtual void setUp( const GetPot& dataFile );

    //! builds the global Epetra map
    /**
       assigns each mesh partition to the corresponding processor, builds the monolithic map
    */
    virtual void setupFluidSolid();

    //! builds the global Epetra map
    /**
       assigns each mesh partition to the corresponding processor, builds the monolithic map with a number of fluxes
       specified from input
    */
    virtual void setupFluidSolid( UInt const fluxes );

    //!@}

    //!@name Public Methods
    //!@{

    //! Transfer a vector to the interface
    /**
       restricts a vector with a monolithic map on the solid interface map
       \param _lambdasolid: vector on the solid interface
       \param _disp: monolithic vector
    */
    void monolithicToInterface(    vector_type& _lambdaSolid, const vector_type& _sol) ;


    //!Transfer a vector to a subdomain
    /**
       restricts a vector with a monolithic map on another map that must
       have a sequential numbering (not the interface map)
       \param _dispFluid: vector on the fluid domain
       \param _disp: monolithic vector
    */
    void monolithicToX(const vector_type& _disp, vector_type& _dispFluid, EpetraMap& map, UInt offset=(UInt)0);

    /**
       sets the vector M_solid->dispSolid() to the monolithic solution M_solid->disp() in the solid nodes, to 0 in the fluid nodes
    */
    void setDispSolid(const vector_type& sol);

    /**
       builds the constant part of the monolithic matrix
    */
    void buildSystem();

    //! Initializer for the solution M_un
    /**
       \small initializes the current solution vector. Note: this is not sufficient for the correct initialization
       of bdf!
    */
    void initialize( vector_ptrtype u0)
    {
        M_un=u0;

//         M_BCh_u->merge(*M_BCh_flux);
//         M_BCh_flux.reset();
//         M_BCh_d->merge(*M_BCh_Robin);
//         M_BCh_Robin.reset();

    }


    void mergeBCHandlers()
    {
        M_BCh_u->merge(*M_BCh_flux);
        M_BCh_flux.reset();
        M_BCh_d->merge(*M_BCh_Robin);
        M_BCh_Robin.reset();
    }

#ifdef HAVE_TRILINOS_ANASAZI
    //! Computes the maximum singular value
    /**
       \small Computes the maximum singular value of the preconditioned system \f$P^-1A\f$ where \f$P\f$ is an
       instance of ComposedOperator and \f$A\f$ is the system matrix.
    */
    LifeV::Real& computeMaxSingularValue();
#endif
    //! Computes the normals to the fluid domain
    /**
       \small Computes the normals to the fluid domain in the nodes. It is an example of how
       to use the boundary conditions methods to compute the normal field on
       a surface.
    */
    void computeFNormals( vector_type& normals);


    //!Evaluates the nonlinear residual
    /**
       This class is pure virtual, it depends on which type of monolithic solver is used
       \param res: output
       \param _sol: monolithic solution
       \param iter: current nonLinRichardson (Newton) iteration
    */
    virtual void   evalResidual(vector_type&        res,
                        const vector_type& _sol,
                        const UInt          _iter) = 0 ;

    /**
       solves the Jacobian system
       \param _muk: output, solution at the current Newton step
       \param _res: nonlinear residual
       \param _linearRelTol: not used

       \small The preconditioner type is usually an algebraic additive Schwarz. The following values
       assigned to the field DDBlockPrec in the data file correspond to different variants:

       - DDBlockPrec = AdditiveSchwarz is AAS on a the system matrix
      Only for the Monolithic Geometry Explicit:
       - DDBlockPrec = ComposedDN is AAS on a Dirichlet-Neumann preconditioner using the ComposedOperator strategy
       - DDBlockPrec = ComposedDN2 is AAS on an alternative Dirichlet-Neumann preconditioner using the ComposedOperator strategy
       - DDBlockPrec = ComposedNN is AAS on a Neumann-Neumann preconditioner using the ComposedOperator strategy
       - DDBlockPrec = ComposedDNND is AAS on a Dirichler-Neumamm -- Neumann-Dirichlet preconditioner using the ComposedOperator strategy

       Only for the Geometry Implicit:
       - DDBlockPrec = 9 is AAS on the quasi-newton matrix obtained with the ComposedOperator strategy
       - DDBlockPrec = 10 is AAS on an alternative matrix obtained with the ComposedOperator strategy
       - DDBlockPrec = 11 is AAS on the quasi-newton matrix obtained with the ComposedOperator strategy, composing
       3 preconditioners
       - DDBlockPrec = 12 is AAS on an alternative matrix obtained with the ComposedOperator strategy, composing
       3 preconditioners
    */
    virtual void   solveJac(vector_type&       _muk,
                            const vector_type& _res,
                            const Real       _linearRelTol);

    /**
       updates the meshmotion, advances of a time step
       \param _sol: solution
    */
    virtual void updateSystem();


    /**
       \small initialize with functions
    */
    virtual void initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                             FSIOperator::solid_type::value_type::Function const& p0,
                             FSIOperator::solid_type::value_type::Function const& d0,
                             FSIOperator::solid_type::value_type::Function const& w0,
                             FSIOperator::solid_type::value_type::Function const& df0);

    /**
       \small initialize the mesh displacement
    */
    virtual void initializeMesh(vector_ptrtype fluid_dispOld);

    //@}

    //!@name Set Methods
    //@{
    //!Sets the restart flag
    void  setRestarts( bool restarts ){ M_restarts = restarts; }

    //! returns a non-const pointer to the preconditioner. Can be used either as a setter or a getter.
    prec_ptrtype& precPtr(){ return M_precPtr; }

    //! returns a non-const pointer to the preconditioner. Can be used either as a setter or a getter.
    block_matrix_ptrtype& operatorPtr(){ return M_monolithicMatrix; }

    /**
       \small sets the fluid BCHandler and merges it with the flux BCHandler
    */
    virtual void setFluidBC     ( const fluid_bchandler_type& bc_fluid )
    {
        super::setFluidBC(bc_fluid);
        //bc_fluid->merge(*M_BCh_flux);
    }

    /**
       \small sets the solid BCHandler and merges it with the Robin BCHandler
    */
    virtual void setSolidBC     ( const fluid_bchandler_type& bc_solid )
    {
        super::setSolidBC(bc_solid);
        //bc_solid->merge(*M_BCh_Robin);
    }

    //! initializes the solution by reference (through a shared_ptr)
    /*!
      \param sol: input pointer
    */
    virtual void setSolutionPtr                     ( const vector_ptrtype& sol)=0;

    //! set the solution
    virtual void setSolution( const vector_type& solution ) = 0;

    //!@}
    //!@name Get Methods
    //!@{

    //! Returns the wall stress
    //!@note still not fixed in parallel
//     vector_ptrtype getWS( )
//     {
//         return M_wss;
//     }

    //! returns a boost shared pointer to the preconditioner
    //prec_raw_type & getPrec(){return M_prec.getPrec();}

#ifdef OBSOLETE
    /** get the shape derivatives vector*/
    vector_type getRhsShapeDerivatives(){return *M_rhsShapeDerivatives;}
#endif
    //    const boost::shared_ptr<EpetraMap>& monolithicMap() {return M_monolithicMap;}

    //!get the total dimension of the FS interface
    const UInt getDimInterface() const {return nDimensions*M_monolithicMatrix->getInterface();}

    //! Returns the solution at the previous time step
    vector_ptrtype const&       un(){return M_un;}

    //! Returns true if CE of FI methods are used, false otherwise (GCE)
    //bool const isFullMonolithic(){return M_fullMonolithic;}

    //! Returns the offset assigned to the solid block
    UInt const getOffset(){return M_offset;}


    //!Get the solid displacement from the solution M_un
    /*!
      \param soliddisp: input vector
    */
    void getSolidDisp(vector_type& soliddisp)
    {
        soliddisp.subset(*un(), M_offset);
        soliddisp *= dataFluid()->dataTime()->getTimeStep()*M_solid->rescaleFactor();
    }

    //!Get the solid velocity
    /*!
       fills an input vector with the solid displacement trom the solution M_un.
      \param solidvel: input vector (output solid velocity)
     */
    void getSolidVel(vector_type& solidvel)
    {
        solidvel.subset(M_solid->vel(), M_offset);
        solidvel *= dataFluid()->dataTime()->getTimeStep()*M_solid->rescaleFactor();
    }

    //! Gets the fluid and pressure
    /**
       fills an input vector with the fluid and pressure trom the solution M_un.
       It performs an Import. Thus it works also for the velocity, depending on the map of the input vector
       \param sol: input vector
    */
    void getFluidVelAndPres(vector_type& sol)
    {
        sol = *un();
    }


    //! Returns the monolithic map
    virtual    boost::shared_ptr<EpetraMap>& getCouplingVariableMap(){return M_monolithicMap;}

    //! get the solution vector
    virtual const vector_type& getSolution() const = 0;

    //virtual vector_ptrtype& solutionPtr() = 0;
    //@}


protected:


    //!@name Protected methods
    //!@{

    //! pure virtual: creates the operator (either of type MonolithicGI or MonolithicGE)
    virtual void createOperator( std::string& operType )=0;

    /**
       \small solves the monolithic system, once a solver, a preconditioner and a rhs have been defined.
    */
    void iterateMonolithic(const vector_type& rhs, vector_type& step);

    /**
       \small adds the part due to coupling to the rhs
       \param rhs: right hand side
       \param un: current solution
    */
    void couplingRhs( vector_ptrtype rhs , vector_ptrtype un );


    /**\small evaluates the linear residual
       \param sol: solution vector
       \param rhs: right-hand side
       \param res: the output residual
       \param diagonalScaling: flag stating wether to perform diagonal scaling
    */
    void evalResidual( const vector_type& sol, const vector_ptrtype& rhs,  vector_type& res, bool diagonalScaling=false);

    //!\small says if the preconditioner will be recomputed
    bool recomputePrec() {return(!M_reusePrec || M_resetPrec);}

    //!\small updates the rhs of the solid block.
    void updateSolidSystem(vector_ptrtype& rhsFluidCoupling);

    /**\small scales matrix and rhs
       \param rhs: the output rhs
       \param matrFull: the output matrix*/
    void    diagonalScale(vector_type& rhs, matrixPtr_Type matrFull);

    //!Empty method kept for compatibility with FSIOperator.
    void shiftSolution(){}

    //!Empty method kept for compatibility with FSIOperator.
    void couplingVariableExtrap(vector_ptrtype& /*lambda*/, vector_ptrtype& /*lambdaDot*/, bool& /*firstIter*/)
    { }

    //! Constructs the solid FESpace
    /**
       Creates the solid FESpace with an unpartitioned mesh, necessary step to create the dof interconnections
       at the interface. The solid FESpace will be reset in variablesInit using the partitioned mesh.
       If the interface map is created offline this method is never called.
       \param dOrder: discretization order
     */
    void solidInit(std::string const& dOrder);

    //! Constructs the solid FESpace and initializes the coupling variables at the interface
    /**
       If the mesh is partitioned online the previous FESpace constructed with the unpartitioned mesh is discarded and
       replaced with one using a partitioned mesh.
     */
    void variablesInit(std::string const& dOrder);

    //!
    int  setupBlockPrec( );

#ifdef OBSOLETE
    void setOperator(Epetra_Operator& epetraOperator) {M_linearSolver->setOperator(epetraOperator);}
#endif

    //! assembles the solid problem (the matrix and the rhs due to the time derivative)
    /*
      \param iter: current nonlinear iteration: used as flag to distinguish the first nonlinear iteration from the others
     */
    void assembleSolidBlock(UInt iter, vector_ptrtype& solution);


    //! assembles the fluid problem (the matrix and the rhs due to the time derivative)
    /*
      \param iter: current nonlinear iteration: used as flag to distinguish the first nonlinear iteration from the others
     */
    void assembleFluidBlock(UInt iter, vector_ptrtype& solution);

    //! Updates the right hand side
    /**
       Adds to the rhs the fluid time discretization terms
       \todo this should be handled externally
     */
    void updateRHS();
    //!@}


    //!@name Protected attributes
    //@{
    boost::shared_ptr<EpetraMap>                      M_monolithicMap;
    boost::shared_ptr< EpetraMap >                    M_interfaceMap;///the solid interface map
    boost::shared_ptr<vector_type>                    M_beta;
    boost::shared_ptr<BlockMatrix >                   M_monolithicMatrix;
    //matrixPtr_Type                                    M_precMatrPtr;
    prec_ptrtype                                      M_precPtr;
    boost::shared_ptr<vector_type>                    M_rhsFull;

    fluid_bchandler_type                              M_BCh_flux;
    solid_bchandler_type                              M_BCh_Robin;
    //UInt                                              M_fluxes;
    solid_bchandler_type                              M_BChWSS;
    BCFunctionMixte                                   M_bcfWss;
    //    matrixPtr_Type                                    M_robinCoupling;
    UInt                                              M_offset;
    UInt                                              M_solidAndFluidDim;
    matrixPtr_Type                                    M_fluidBlock;
    matrixPtr_Type                                    M_solidBlock;
    matrixPtr_Type                                    M_solidBlockPrec;
    matrixPtr_Type                                    M_robinCoupling; //uninitialized if not needed
    boost::shared_ptr<solver_type>                    M_linearSolver;
    boost::shared_ptr<vector_type>                    M_numerationInterface;
    std::vector<fluid_bchandler_type>                 M_BChs;
    std::vector<boost::shared_ptr<FESpace<mesh_type, EpetraMap> > >          M_FESpaces;
    bool                                              M_diagonalScale;
    bool                                              M_reusePrec;//!\todo to move to private
    bool                                              M_resetPrec;//!\todo to move to private
    UInt                                              M_maxIterSolver;//!\todo to move to private
    bool                                              M_restarts;
    //@}

private:
    //!@name Private methods
    //!@{
    void initialize( vector_type const& u0, vector_type const& p0, vector_type const& d0);
    //!@}
    //!@name Private attributes
    //!@{
    //! operator \f$P^{-1}AA^TP^{-T}\f$, where P is the preconditioner and A is the monolithic matrix
    boost::shared_ptr<ComposedOperator<Epetra_Operator> > M_preconditionedSymmetrizedMatrix;
#ifdef OBSOLETE
    boost::shared_ptr<vector_type>                    M_rhsShapeDerivatives;
#endif
    static bool                                       reg;
    //@}
};

class WRONG_PREC_EXCEPTION
{
public:

    //! Exception thrown when a wrong preconditioning strategy is set
    WRONG_PREC_EXCEPTION(){}
    virtual ~WRONG_PREC_EXCEPTION(){}

};

}
#endif
