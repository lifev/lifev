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

#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/fsi/solver/FSIOperator.hpp>

#include <lifev/core/algorithm/PreconditionerComposed.hpp>
#include <lifev/core/algorithm/ComposedOperator.hpp>
#ifdef HAVE_TRILINOS_ANASAZI
#include <lifev/core/algorithm/EigenSolver.hpp>
#endif

#include <lifev/fsi/solver/MonolithicBlockMatrix.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposed.hpp>
#ifdef HAVE_NS_PREC
#include <life/lifealg/PreconditionerPCD.hpp>
#endif

namespace LifeV
{

class WRONG_PREC_EXCEPTION;

//! FSIMonolithic.hpp pure virtual class containing the core methods of the FSIMonolithic FSI solver
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
 * Exact Newton method and full implicit time discretization are implemented in the FSIMonolithicGI class.
 */

class FSIMonolithic : public FSIOperator
{
public:

    //!@name Typedefs
    //@{

    typedef FSIOperator                                               super_Type;
    typedef FSIOperator::fluidPtr_Type::value_type::matrix_Type       matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                            matrixPtr_Type;
    typedef MonolithicBlock                                           prec_Type;
    typedef boost::shared_ptr<prec_Type>                              precPtr_Type;
    typedef MonolithicBlockMatrix                                     blockMatrix_Type;
    typedef boost::shared_ptr<blockMatrix_Type>                       blockMatrixPtr_Type;
    typedef FactorySingleton< Factory< FSIMonolithic, std::string > > factory_Type;
    typedef SolverAztecOO                                             solver_Type;

    //@}

    // constructors

    //!@name Constructors, Destructor
    //!@{
    FSIMonolithic();

    ~FSIMonolithic();
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

#ifdef HAVE_HDF5
    /**
       reads the interface map between the fluid and solid meshes from file.
     */
    void setupDOF( meshFilter_Type& filterMesh );
#endif

    //! sets the parameters from data file
    /**
       Calls the setup of the fluid problem and the setUp method.
     */
    virtual void setupSystem( );

    //     /** stores the data file into a member */
    //     virtual void setDataFile( const GetPot& dataFile );

    //! setup method for the FSIMonolithic solver
    /**
       sets some parameters specific to the FSIMonolithic class
     */
    virtual void setup( const GetPot& dataFile );

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

    //! Transfers a vector to the interface
    /**
       restricts a vector with a monolithic map on the solid interface map
       \param lambdaSolid: vector on the solid interface
       \param disp: monolithic vector
     */
    void monolithicToInterface(    vector_Type& lambdaSolid, const vector_Type& sol) ;


    //!Transfers a vector to a subdomain
    /**
       restricts a vector with a monolithic map on another map that must
       have a sequential numbering (not the interface map)
       \param disp: monolithic vector
       \param dispFluid: vector on the fluid domain
       \param map: MapEpetra of the part of vector that we want to transfer
       \param offset: offset for the monolithic vector (also alpplied to the input map)
     */
    void monolithicToX(const vector_Type& disp, vector_Type& dispFluid, MapEpetra& map, UInt offset=(UInt)0);


    /**
       builds the constant part of the monolithic matrix
     */
    void buildSystem();


    //! Merges the flux boundary conditions into the fluid BCHandler
    /*!two separate BCHandlers are initially created for the flux-type boundary conditions, these are later merged with the fluid BCHandler
      automatically, using this method
     */
    void mergeBCHandlers()
    {
        if(M_BCh_u.get())
            M_BCh_u->merge(*M_BCh_flux);
        else // in this case only fluxes are imposed on the fluid
            M_BCh_u=M_BCh_flux;
        M_BCh_flux.reset();
    }

#ifdef HAVE_TRILINOS_ANASAZI
    //! Computes the maximum singular value
    /**
       \small Computes the maximum singular value of the preconditioned system \f$P^{-1}A\f$ where \f$P\f$ is an
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
    void computeFluidNormals( vector_Type& normals);


    //!Evaluates the nonlinear residual
    /**
       This class is pure virtual, it depends on which type of monolithic solver is used
       \param res: output
       \param _sol: monolithic solution
       \param iter: current NonLinearRichardson (Newton) iteration
     */
    virtual void   evalResidual(vector_Type&        res,
                                const vector_Type& sol,
                                const UInt         iter) = 0 ;

    /**
       solves the Jacobian system
       \param muk: output, solution at the current Newton step
       \param res: nonlinear residual
       \param linearRelTol: not used

       \small The preconditioner type is usually an algebraic additive Schwarz. The following values
       assigned to the field DDBlockPrec in the data file correspond to different variants:

       Only for the FSIMonolithic Geometry Explicit:
       - DDBlockPrec = AdditiveSchwarz is AAS on a the system matrix
       - DDBlockPrec = MonolithicBlockComposedDN is AAS on a Dirichlet-Neumann preconditioner using the ComposedOperator strategy
       - DDBlockPrec = ComposedDN2 is AAS on an alternative Dirichlet-Neumann preconditioner using the ComposedOperator strategy
       - DDBlockPrec = MonolithicBlockComposedNN is AAS on a Neumann-Neumann preconditioner using the ComposedOperator strategy
       - DDBlockPrec = MonolithicBlockComposedDNND is AAS on a Dirichler-Neumamm -- Neumann-Dirichlet preconditioner using the ComposedOperator strategy

       Only for the Geometry Implicit:
       - DDBlockPrec = AdditiveSchwarzGI is AAS on the whole matrix.
       - DDBlockPrec = MonolithicBlockComposedDNDGI is AAS on the quasi-newton matrix obtained with the ComposedOperator strategy. Split in 3 factors.
       - DDBlockPrec = ComposedDND2GI is AAS on the quasi-newton matrix obtained with the ComposedOperator strategy. Split in 3 factors.
       - DDBlockPrec = MonolithicBlockComposedDNGI is AAS on the full Jacobian matrix obtained with the ComposedOperator strategy by neglecting part
       of the fluid--structure coupling, split in 3 factors
       - DDBlockPrec = MonolithicBlockComposedDN2GI is AAS on an alternative matrix obtained with the previous strategy
     */
    virtual void   solveJac(vector_Type&       muk,
                            const vector_Type& res,
                            const Real       linearRelTol);

    /**
       updates the meshmotion, advances of a time step
       \param _sol: solution
     */
    virtual void updateSystem();

    //! activates the computation of the wall stress on the boundary with a specified flag.
    /**
       Notice that the specified flag must be in the coupling fluid-structure interface
     */
    void enableStressComputation(UInt  flag);

    /**
       Computes the stress on the coupling boundary (the traction vector)
     */
    vectorPtr_Type computeStress();


    //@}

    //!@name Set Methods
    //@{

    //! returns a non-const pointer to the preconditioner. Can be used either as a setter or a getter.
    precPtr_Type& precPtrView(){ return M_precPtr; }

    //! returns a non-const pointer to the preconditioner. Can be used either as a setter or a getter.
    blockMatrixPtr_Type& operatorPtrView(){ return M_monolithicMatrix; }

    /**
       \small sets the solid BCHandle
     */
    virtual void setSolidBC     ( const fluidBchandlerPtr_Type& bc_solid )
    {
        super_Type::setSolidBC(bc_solid);
        //bc_solid->merge(*M_BCh_Robin);
    }

    //! initializes the solution by reference (through a shared_ptr)
    /*!
      \param sol: input pointer
     */

    void setFluidBC     ( const fluidBchandlerPtr_Type& bc_fluid )
    {
        super_Type::setFluidBC(bc_fluid);
        if(M_BChs.size())
        {
            UInt nfluxes(M_BChs[1]->numberOfBCWithType(Flux));
            //M_substituteFlux.reset(new std::vector<bool>(nfluxes))
            M_fluxOffset.resize(nfluxes);
            M_BCFluxNames = M_BChs[1]->findAllBCWithType(Flux);
            for (UInt i=0; i<nfluxes; ++i)
            {
                const BCBase* bc = M_BChs[1]->findBCWithName(M_BCFluxNames[i]);
                M_fluxOffset[i]=bc->offset();
            }
            M_BChs[1]=bc_fluid;
            M_monolithicMatrix->setConditions(M_BChs);
            M_precPtr->setConditions(M_BChs);
        }
    }


    //!get the total dimension of the FS interface
    UInt dimInterface() const {return nDimensions*M_monolithicMatrix->interface();}

    //! Returns true if CE of FI methods are used, false otherwise (GCE)
    //bool const isFullMonolithic(){return M_fullMonolithic;}

    //! Returns the offset assigned to the solid block
    UInt  offset() const {return M_offset;}

    //!Get the solid displacement from the solution
    /*!
      \param solidDisplacement: input vector
     */
    void exportSolidDisplacement( vector_Type& solidDisplacement )
    {
        solidDisplacement.subset( M_fluidTimeAdvance->singleElement(0), M_offset );
        solidDisplacement *= M_solid->rescaleFactor();
    }

    //!Get the solid velocity
    /*!
      fills an input vector with the solid displacement from the solution.
      \param solidVelocity: input vector (output solid velocity)
     */
    void exportSolidVelocity( vector_Type& solidVelocity )
    {
        solidVelocity = M_solidTimeAdvance->firstDerivative();
        solidVelocity *= M_solid->rescaleFactor();
    }

    //!Get the solid accelration
    /*!
      fills an input vector with the solid displacement from the solution.
      \param solidVelocity: input vector (output solid acceleration)
     */
    void exportSolidAcceleration( vector_Type& solidAcceleration )
    {
        solidAcceleration = M_solidTimeAdvance->secondDerivative();
        solidAcceleration *= M_solid->rescaleFactor();
    }

    //! Export the fluid velocity by copying it to an external vector
    /*!
     * @param fluidVelocity vector to be filled with the fluid velocity
     */
    void exportFluidVelocity( vector_Type& fluidVelocity ) { fluidVelocity.subset( M_fluidTimeAdvance->singleElement(0), 0 ); }

    //! Export the fluid pressure by copying it to an external vector
    /*!
     * @param fluidPressure vector to be filled with the fluid pressure
     */
    void exportFluidPressure( vector_Type& fluidPressure ) { fluidPressure.subset( M_fluidTimeAdvance->singleElement(0), static_cast<UInt> (3 * M_uFESpace->dof().numTotalDof() ) ); }

    //! Gets the fluid and pressure
    /**
       fills an input vector with the fluid and pressure from the solution M_un.
       It performs a trilinos import. Thus it works also for the velocity, depending on the map of the input vector
       \param fluidVelocityandPressure: input vector
     */
    void exportFluidVelocityAndPressure( vector_Type& fluidVelocityAndPressure )
    {
        fluidVelocityAndPressure = M_fluidTimeAdvance->singleElement(0);
    }

    //! Returns the monolithic map
    virtual boost::shared_ptr<MapEpetra>& couplingVariableMap() { return M_monolithicMap; }

    //! get the solution vector
    virtual const vector_Type& solution() const = 0;

    //! Update the solution after NonLinearRichardson is called.
    /*!
     *  Here it is used also to update the velocity for the post-processing.
     */
    virtual void updateSolution( const vector_Type& solution )
    {
        this->M_fluidTimeAdvance->shiftRight(solution);
        if(M_data->dataFluid()->conservativeFormulation())
            this->M_fluidMassTimeAdvance->shiftRight(M_fluid->matrixMass()*solution);
        this->M_solidTimeAdvance->shiftRight(solution);
    }

    //! Updates the right hand side
    /**
       Adds to the rhs the fluid time discretization terms
       \todo this should be handled externally
     */
    void updateRHS();

    //! Set vectors for restart
    /*!
     *  Set vectors for restart
     */
    void setVectorInStencils( const vectorPtr_Type& vel,
                              const vectorPtr_Type& pressure,
                              const vectorPtr_Type& solidDisp,
                              //const vectorPtr_Type& fluidDisp,
                              const UInt iter);

    void setFluidVectorInStencil( const vectorPtr_Type& vel, const vectorPtr_Type& pressure, const UInt iter);

    void setSolidVectorInStencil( const vectorPtr_Type& solidDisp, const UInt iter);

    virtual void setALEVectorInStencil( const vectorPtr_Type& fluidDisp, const UInt iter,const bool lastVector) = 0;

    void finalizeRestart();

    //@}


protected:


    //!@name Protected methods
    //!@{

    //! pure virtual: creates the operator (either of type FSIMonolithicGI or FSIMonolithicGE)
    virtual void createOperator( std::string& operType )=0;

    /**
       \small solves the monolithic system, once a solver, a preconditioner and a rhs have been defined.
     */
    void iterateMonolithic(const vector_Type& rhs, vector_Type& step);

    /**
       \small adds the part due to coupling to the rhs
       \param rhs: right hand side
       \param un: current solution
     */
    void couplingRhs( vectorPtr_Type rhs);


    /**\small evaluates the linear residual
       \param sol: solution vector
       \param rhs: right-hand side
       \param res: the output residual
       \param diagonalScaling: flag stating wether to perform diagonal scaling
     */
    void evalResidual( const vector_Type& sol, const vectorPtr_Type& rhs,  vector_Type& res, bool diagonalScaling=false);

    //!\small says if the preconditioner will be recomputed
    bool recomputePrec() {return(!M_reusePrec || M_resetPrec);}

    //!\small updates the rhs of the solid block.
    void updateSolidSystem(vectorPtr_Type& rhsFluidCoupling);

    /**\small scales matrix and rhs
       \param rhs: the output rhs
       \param matrFull: the output matrix*/
    void diagonalScale(vector_Type& rhs, matrixPtr_Type matrFull);


    //! Constructs the solid FESpace
    /**
       Creates the solid FESpace with an unpartitioned mesh, necessary step to create the dof interconnections
       at the interface. The solid FESpace will be reset in variablesInit using the partitioned mesh.export
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
    virtual void setupBlockPrec();


    //! assembles the solid problem (the matrix and the rhs due to the time derivative)
    /*
      \param iter: current nonlinear iteration: used as flag to distinguish the first nonlinear iteration from the others
      \param solution: current solution, used for the time advance implementation, and thus for the update of the right hand side
     */
    void assembleSolidBlock(UInt iter, const vector_Type& solution);


    //! assembles the fluid problem (the matrix and the rhs due to the time derivative)
    /*
      \param iter: current nonlinear iteration: used as flag to distinguish the first nonlinear iteration from the others
      \param solution: current solution, used for the time advance implementation, and thus for the update of the right hand side
     */
    void assembleFluidBlock(UInt iter, const vector_Type& solution);

    //! Checks if the flux bcs changed during the simulation, e.g. if a flux b.c. has been changed to Natural
    //! (this can be useful when modeling valves)
    /**
       When the fluxes bcs changed a '1' is added in the line corresponding to the Lagrange multiplier. This method must
       be called for both operator and preconditioner
     */
    void checkIfChangedFluxBC( precPtr_Type oper );

    //!@}


    //!@name Protected attributes
    //@{
    boost::shared_ptr<MapEpetra>                      M_monolithicMap;
    boost::shared_ptr< MapEpetra >                    M_interfaceMap;///the solid interface map
    boost::shared_ptr<vector_Type>                    M_beta;
    boost::shared_ptr<MonolithicBlockMatrix >                   M_monolithicMatrix;
    precPtr_Type                                      M_precPtr;
    boost::shared_ptr<vector_Type>                    M_rhsFull;

    fluidBchandlerPtr_Type                            M_BCh_flux;
    solidBchandlerPtr_Type                            M_BChWS;
    BCFunctionRobin                                   M_bcfWs;
    UInt                                              M_offset;
    UInt                                              M_solidAndFluidDim;
    FSIOperator::fluidPtr_Type::value_type::matrixPtr_Type M_fluidBlock;
    matrixPtr_Type                                    M_solidBlockPrec;
    matrixPtr_Type                                    M_robinCoupling; //uninitialized if not needed
    matrixPtr_Type                                    M_boundaryMass;
    boost::shared_ptr<solver_Type>                    M_linearSolver;
    boost::shared_ptr<vector_Type>                    M_numerationInterface;
    std::vector<fluidBchandlerPtr_Type>                 M_BChs;
    std::vector<boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > >          M_FESpaces;
    bool                                              M_diagonalScale;
    bool                                              M_reusePrec;//!\todo to move to private
    bool                                              M_resetPrec;//!\todo to move to private
    Int                                               M_maxIterSolver;//!\todo to move to private
    bool                                              M_restarts;
    //@}

private:
    //!@name Private attributes
    //!@{
    //! operator \f$P^{-1}AA^TP^{-T}\f$, where P is the preconditioner and A is the monolithic matrix
    boost::shared_ptr<ComposedOperator<Epetra_Operator> > M_preconditionedSymmetrizedMatrix;
    boost::shared_ptr<vector_Type>                    M_stress;
    UInt                                              M_fluxes;
    std::vector<bcName_Type>                          M_BCFluxNames;
    std::vector<UInt>                                 M_fluxOffset;

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
