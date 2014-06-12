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
    @file
    @brief File containing the implementation of the file OseenSolver.hpp

    @author Davide Forti <davide.forti@epfl.ch>

    @date 28-01-2014

 */

#include <lifev/core/LifeV.hpp>

namespace LifeV
{
    
class buildVector
{
public:
    typedef VectorSmall<3> return_Type;
    
    inline return_Type operator() (Real a, Real b, Real c)
    {
        VectorSmall<3> o;
        o[0] = a;
        o[1] = b;
        o[2] = c;
        return o;
    }
    
    buildVector() {}
    buildVector (const buildVector&) {}
    ~buildVector() {}
};
    
// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              boost::shared_ptr<Epetra_Comm>& communicator,
              const Int                       lagrangeMultiplier ) :
	  M_oseenData       ( dataType ),
	  M_velocityFESpace        ( velocityFESpace ),
	  M_pressureFESpace        ( pressureFESpace ),
	  M_Displayer              ( communicator ),
	  M_fluxMap                ( lagrangeMultiplier, communicator),
	  M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultiplier),
	  M_velocityMatrixMass     ( ),
	  M_pressureMatrixMass     ( ),
	  M_matrixStokes           ( ),
	  M_matrixNoBC             ( ),
	  M_matrixStabilization    ( ),
	  M_rightHandSideNoBC      ( ),
	  M_solution               ( new vector_Type ( M_localMap ) ),
	  M_residual               ( new vector_Type (M_localMap ) ),
	  M_linearSolver           ( new linearSolver_Type (communicator) ),
	  M_steady                 ( ),
	  M_postProcessing         ( new PostProcessingBoundary<mesh_Type> ( M_velocityFESpace.mesh(),
			  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 &M_velocityFESpace.feBd(),
			  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 &M_velocityFESpace.dof(),
			  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 &M_pressureFESpace.feBd(),
			  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 &M_pressureFESpace.dof(),
			  	  	  	  	  	  	  	  	  	  	  	  	  	  	  	 M_localMap ) ),
	  M_stabilization          ( false ),
	  M_reuseStabilization     ( false ),
	  M_resetStabilization     ( false ),
	  M_iterReuseStabilization ( -1 ),
	  //        M_ipStabilization        ( M_velocityFESpace.mesh(),
	  //                                   M_velocityFESpace.dof(),
	  //                                   M_velocityFESpace.refFE(),
	  //                                   M_velocityFESpace.feBd(),
	  //                                   M_velocityFESpace.qr(),
	  //                                   0., 0., 0.,
	  //                                   M_oseenData->viscosity() ),
	  M_betaFunction           ( 0 ),
	  M_divBetaUv              ( false ),
	  M_stiffStrain            ( false ),
	  M_diagonalize            ( false ),
	  M_count                  ( 0 ),
	  M_recomputeMatrix        ( false ),
	  M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), velocityFESpace.fieldDim() ),
	  M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), velocityFESpace.fieldDim() ),
	  M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
	  M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
			  	  	  	  	  	  M_velocityFESpace.fe().nbFEDof(), 0, velocityFESpace.fieldDim() ),
	  M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), 0,
	  M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
	  M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_un                     ( new vector_Type (M_localMap) ),
      M_velocityPreviousTimestep (new vector_Type (M_velocityFESpace.map()) ),
      M_pressurePreviousTimestep (new vector_Type (M_pressureFESpace.map()) ),
      M_fespaceUETA            ( new ETFESpace_velocity(M_velocityFESpace.mesh(), &(M_velocityFESpace.refFE()), communicator)),
      M_fespacePETA            ( new ETFESpace_pressure(M_pressureFESpace.mesh(), &(M_pressureFESpace.refFE()), communicator)),
      M_supgStabilization      (new StabilizationSUPG<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      M_VMSLESStabilization    (new StabilizationVMSLES<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace))
{
    *M_solution *= 0;
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              boost::shared_ptr<Epetra_Comm>& communicator,
              MapEpetra                       monolithicMap,
              UInt                            offset ) :
	  M_oseenData              ( dataType ),
	  M_velocityFESpace        ( velocityFESpace ),
	  M_pressureFESpace        ( pressureFESpace ),
	  M_Displayer              ( communicator ),
      M_fluxMap                ( offset, communicator),
	  M_localMap               ( monolithicMap ),
	  M_velocityMatrixMass     ( ),
	  M_matrixStokes           ( ),
	  M_matrixNoBC             ( ),
	  M_matrixStabilization    ( ),
	  M_rightHandSideNoBC      ( ),
	  M_solution               ( ),
	  M_residual               (  ),
	  M_linearSolver           ( ),
	  M_postProcessing         ( new PostProcessingBoundary<mesh_Type> (M_velocityFESpace.mesh(),
																		&M_velocityFESpace.feBd(),
																		&M_velocityFESpace.dof(),
																		&M_pressureFESpace.feBd(),
																		&M_pressureFESpace.dof(),
																		M_localMap ) ),
	  M_stabilization          ( false ),
	  M_reuseStabilization     ( false ),
	  M_resetStabilization     ( false ),
	  M_iterReuseStabilization ( -1 ),
	  M_betaFunction           ( 0 ),
	  M_divBetaUv              ( false ),
	  M_stiffStrain            ( false ),
	  M_diagonalize            ( false ),
	  M_count                  ( 0 ),
	  M_recomputeMatrix        ( false ),
	  M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
	  M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
	  M_elementMatrixPreconditioner                 ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
	  M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
								  M_velocityFESpace.fe().nbFEDof(), 0, M_velocityFESpace.fieldDim() ),
	  M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), 0,
								 M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
	  M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
	  M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
	  M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
	  M_un                     ( /*new vector_Type(M_localMap)*/ ),
      M_velocityPreviousTimestep (new vector_Type (M_velocityFESpace.map()) ),
      M_pressurePreviousTimestep (new vector_Type (M_pressureFESpace.map()) ),
      M_fespaceUETA            ( new ETFESpace_velocity(M_velocityFESpace.mesh(), &(M_velocityFESpace.refFE()), communicator)),
      M_fespacePETA            ( new ETFESpace_pressure(M_pressureFESpace.mesh(), &(M_pressureFESpace.refFE()), communicator)),
      M_supgStabilization      (new StabilizationSUPG<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      M_VMSLESStabilization    (new StabilizationVMSLES<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace))
{
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              const std::vector<Int>&         lagrangeMultipliers,
              boost::shared_ptr<Epetra_Comm>& communicator ) :
	  M_oseenData              ( dataType ),
	  M_velocityFESpace        ( velocityFESpace ),
	  M_pressureFESpace        ( pressureFESpace ),
	  M_Displayer              ( communicator ),
	  M_fluxMap                ( lagrangeMultipliers, communicator),
	  M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultiplier),
	  M_velocityMatrixMass     ( ),
	  M_matrixStokes           ( ),
	  M_matrixNoBC             ( ),
	  M_matrixStabilization    ( ),
	  M_rightHandSideNoBC      ( ),
	  M_solution               ( new vector_Type ( M_localMap ) ),
	  M_residual               (  ),
	  M_linearSolver           ( new linearSolver_Type (communicator) ),
	  M_postProcessing         ( new PostProcessingBoundary<mesh_Type> (M_velocityFESpace.mesh(),
																		&M_velocityFESpace.feBd(),
																		&M_velocityFESpace.dof(),
																		&M_pressureFESpace.feBd(),
																		&M_pressureFESpace.dof(),
																		M_localMap ) ),
	  M_stabilization          ( false ),
	  M_reuseStabilization     ( false ),
	  M_resetStabilization     ( false ),
	  M_iterReuseStabilization ( -1 ),
	  M_betaFunction           ( 0 ),
	  M_divBetaUv              ( false ),
	  M_stiffStrain            ( false ),
	  M_diagonalize            ( false ),
	  M_count                  ( 0 ),
	  M_recomputeMatrix        ( false ),
	  M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
	  M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
	  M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
	  M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
								  M_velocityFESpace.fe().nbFEDof(), 0, M_velocityFESpace.fieldDim() ),
	  M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), 0,
								 M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
	  M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_un                     ( new vector_Type (M_localMap) ),
      M_velocityPreviousTimestep (new vector_Type (M_velocityFESpace.map()) ),
      M_pressurePreviousTimestep (new vector_Type (M_pressureFESpace.map()) ),
      M_fespaceUETA            ( new ETFESpace_velocity(M_velocityFESpace.mesh(), &(M_velocityFESpace.refFE()), communicator)),
      M_fespacePETA            ( new ETFESpace_pressure(M_pressureFESpace.mesh(), &(M_pressureFESpace.refFE()), communicator)),
      M_supgStabilization      (new StabilizationSUPG<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      M_VMSLESStabilization    (new StabilizationVMSLES<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace))
{
    *M_solution *= 0;
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
~OseenSolver()
{

}


// ===================================================
// Methods
// ===================================================

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::setUp ( const GetPot& dataFile )
{
    if (M_linearSolver.get() )
    {
        M_linearSolver->setupPreconditioner ( dataFile, "fluid/prec" );
        M_linearSolver->setDataFromGetPot ( dataFile, "fluid/solver" );
    }

    M_stabilization = dataFile ( "fluid/stabilization/use",  /*(&M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ) ,*/ false);

    // If using P1-P1 the use of the stabilization is necessary
    //if(&M_velocityFESpace.refFE() == &M_pressureFESpace.refFE())
    //	M_stabilization = true;

    M_steady        = dataFile ( "fluid/miscellaneous/steady", 0 );
    
    if (M_stabilization)
    {
    	if(M_oseenData->stabilizationType() == "IP")
    	{
    		M_gammaBeta     = dataFile ( "fluid/ipstab/gammaBeta",  0. );
    		M_gammaDiv      = dataFile ( "fluid/ipstab/gammaDiv",   0. );
    		M_gammaPress    = dataFile ( "fluid/ipstab/gammaPress", 0. );
    		M_reuseStabilization     = dataFile ( "fluid/ipstab/reuse", false );

    		M_ipStabilization.setGammaBeta ( M_gammaBeta );
    		M_ipStabilization.setGammaDiv  ( M_gammaDiv );
    		M_ipStabilization.setGammaPress ( M_gammaPress );

    		if (M_linearSolver.get() )
    			M_iterReuseStabilization = dataFile ( "fluid/ipstab/max_iter_reuse", static_cast<Int> ( M_linearSolver->maxNumIterations() * 8. / 10. ) );
    	}
    	else if (M_oseenData->stabilizationType() == "SUPG")
    	{
    		int vel_order = dataFile ( "fluid/space_discretization/vel_order", 1 );
    		M_supgStabilization->setConstant ( vel_order );
    		M_supgStabilization->setETvelocitySpace(M_fespaceUETA);
    		M_supgStabilization->setETpressureSpace(M_fespacePETA);
    		M_supgStabilization->setCommunicator(M_velocityFESpace.map().commPtr());
    		M_supgStabilization->setDensity(M_oseenData->density());
    		M_supgStabilization->setViscosity(M_oseenData->viscosity());
    		M_supgStabilization->setTimeStep(M_oseenData->dataTime()->timeStep());

    	}
        else if (M_oseenData->stabilizationType() == "VMSLES")
    	{
            *M_velocityPreviousTimestep *= 0;
            *M_pressurePreviousTimestep *= 0;
    		int vel_order = dataFile ( "fluid/space_discretization/vel_order", 1 );
    		M_VMSLESStabilization->setConstant ( vel_order );
    		M_VMSLESStabilization->setETvelocitySpace(M_fespaceUETA);
    		M_VMSLESStabilization->setETpressureSpace(M_fespacePETA);
    		M_VMSLESStabilization->setCommunicator(M_velocityFESpace.map().commPtr());
    		M_VMSLESStabilization->setDensity(M_oseenData->density());
    		M_VMSLESStabilization->setViscosity(M_oseenData->viscosity());
    		M_VMSLESStabilization->setTimeStep(M_oseenData->dataTime()->timeStep());
    	}
	}
    // Energetic stabilization term
    M_divBetaUv   = dataFile ( "fluid/space_discretization/div_beta_u_v", false);
    // Enable grad( u )^T in stress tensor
    M_stiffStrain = dataFile ( "fluid/space_discretization/stiff_strain", false);
    M_diagonalize = dataFile ( "fluid/space_discretization/diagonalize", 1. );

    //    M_linearSolver.setAztecooPreconditioner( dataFile, "fluid/solver" );

}


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
initialize ( const function_Type& velocityFunction, const function_Type& pressureFunction )
{
    vector_Type velocityInitialGuess ( M_velocityFESpace.map() );
    M_velocityFESpace.interpolate ( velocityFunction,
                                    velocityInitialGuess,
                                    M_oseenData->dataTime()->time() );

    vector_Type pressureInitialGuess ( M_pressureFESpace.map() );
    M_pressureFESpace.interpolate ( pressureFunction,
                                    pressureInitialGuess,
                                    M_oseenData->dataTime()->time() );

    initialize ( velocityInitialGuess, pressureInitialGuess );
}


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
initialize ( const vector_Type& velocityInitialGuess, const vector_Type& pressureInitialGuess )
{

    *M_solution = velocityInitialGuess;
    *M_un = velocityInitialGuess;
    M_solution->add ( pressureInitialGuess, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof() );
    M_un->add ( pressureInitialGuess, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof() );

}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
initialize ( const vector_Type& velocityAndPressure )
{

    *M_un = velocityAndPressure;
    if ( M_solution.get() )
    {
        *M_solution = velocityAndPressure;
    }

}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::buildSystem()
{
	// New ETA PART
	/*
	 *  TESTING THE ASSEMBLY USING ETA
	 *
	 */

	M_Displayer.leaderPrint ( "  F-  Computing constant matrices ...          " );
	LifeChrono chrono;
	chrono.start();

	M_velocityMatrixMass.reset  ( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
	M_matrixStokes.reset(new matrix_block_Type( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ));
    
    *M_velocityMatrixMass *= 0;
	*M_matrixStokes *= 0;

	{
		using namespace ExpressionAssembly;
		
        if ( !M_steady )
		{
			integrate(
					elements(M_fespaceUETA->mesh()),
					M_velocityFESpace.qr(),
					M_fespaceUETA,
					M_fespaceUETA,
					M_oseenData->density() * dot(phi_i, phi_j)
			) >> M_velocityMatrixMass->block(0,0);
			
            M_velocityMatrixMass->globalAssemble();
            
            M_matrixMass.reset(new matrix_Type(M_localMap));
            *M_matrixMass += *M_velocityMatrixMass;
            M_matrixMass->globalAssemble();
		}

        if ( M_stiffStrain )
        {
            integrate(
                    elements(M_fespaceUETA->mesh()),
                    M_velocityFESpace.qr(),
                    M_fespaceUETA,
                    M_fespaceUETA,
                    M_oseenData->viscosity() * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
            ) >> M_matrixStokes->block(0,0);
        }
        else
        {
            integrate(
                    elements(M_fespaceUETA->mesh()),
                    M_velocityFESpace.qr(),
                    M_fespaceUETA,
                    M_fespaceUETA,
                    M_oseenData->viscosity() * dot( grad(phi_i), grad(phi_j) )
            ) >> M_matrixStokes->block(0,0);
        }
		

		integrate(
				elements(M_fespaceUETA->mesh()),
				M_velocityFESpace.qr(),
				M_fespaceUETA,
				M_fespacePETA,
				value(-1.0) * phi_j * div(phi_i)
		) >> M_matrixStokes->block(0,1);

		integrate(
				elements(M_fespaceUETA->mesh()),
				M_velocityFESpace.qr(),
				M_fespacePETA,
				M_fespaceUETA,
				phi_i * div(phi_j)
		) >> M_matrixStokes->block(1,0);
	}
    
    M_matrixStokes->globalAssemble();
    
	chrono.stop();
	M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
	
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
updateSystem ( const Real         alpha,
               const vector_Type& u_star,
               const vector_Type& sourceVector )
{
    if ( M_matrixNoBC.get() )
    {
        M_matrixNoBC.reset ( new matrix_Type ( M_localMap, M_matrixNoBC->meanNumEntries() ) );
    }
    else
    {
        M_matrixNoBC.reset ( new matrix_Type ( M_localMap ) );
    }

    updateSystem ( alpha, u_star, sourceVector, M_matrixNoBC, *M_un );
    if ( alpha != 0. )
    {
        M_matrixNoBC->globalAssemble();
    }

}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
updateSystem ( const Real          alpha,
		       const vector_Type&  u_star,
               const vector_Type&  rightHandSide,
               matrixPtr_Type      matrixNoBC,
               const vector_Type&  un )
{
    M_Displayer.leaderPrint ( "  F-  Updating mass term on right hand side... " );
    LifeChrono chrono;
    chrono.start();
    
    updateRightHandSide ( rightHandSide );
    
    chrono.stop();
    M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );
    
    M_matrixNoBC_block.reset(new matrix_block_Type(M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap));
    
    if ( M_recomputeMatrix )
        buildSystem();
    
    //! managing the convective term : semi-implicit approximation of the convective term
    
    Real normInf;
    u_star.normInf ( &normInf );
    
    MatrixSmall<3, 3> Eye;
    Eye *= 0.0;
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;
    
    // u_star: extrapolation of the velocity
    // NOTE:   for ALE formulation it has to be already: extrapolation of fluid velocity - extrapolation of fluid mesh velocity

    vector_Type u_starRepeated ( u_star, Repeated );
    vector_Type unRepeated ( un, Repeated );

    vector_Type wRepeated ( un, Repeated );
    wRepeated -= u_starRepeated;

    M_convectiveMatrix.reset  ( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
    *M_convectiveMatrix *= 0;
    
    // only if the extrapolation of the velocity differs from zero we update the convective term
    if ( normInf != 0. )
    {
        M_convectiveMatrix.reset  ( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
        *M_convectiveMatrix *= 0;
        
        if(M_oseenData->conservativeFormulation())
        {
            using namespace ExpressionAssembly;
            integrate(
                elements(M_fespaceUETA->mesh()), // Mesh
                M_velocityFESpace.qr(),          // QR
                M_fespaceUETA,
                M_fespaceUETA,
                dot(grad(phi_j) * M_oseenData->density() * value(M_fespaceUETA, u_starRepeated), phi_i)                              // semi-implicit convective term
                + 0.5 * M_oseenData->density() * dot ( value ( Eye ) , grad(M_fespaceUETA, u_starRepeated) ) * dot( phi_j , phi_i )  // consistency term
                - M_oseenData->density() * dot ( value ( Eye ) , grad(M_fespaceUETA, wRepeated) ) * dot( phi_j , phi_i )             // conservative formulation
                      )
            >> M_convectiveMatrix->block(0,0);
        }
        else
        {
            using namespace ExpressionAssembly;
            integrate(
                elements(M_fespaceUETA->mesh()), // Mesh
                M_velocityFESpace.qr(),          // QR
                M_fespaceUETA,
                M_fespaceUETA,
                dot( M_oseenData->density() * value(M_fespaceUETA, u_starRepeated)*grad(phi_j), phi_i) // semi-implicit convective term
                //+ 0.5 * M_oseenData->density() * dot ( value ( Eye ) , grad(M_fespaceUETA, u_starRepeated) ) * dot( phi_j , phi_i )   // consistency term
                      )
            >> M_convectiveMatrix->block(0,0);
        }
        
        // Stabilising term: -rho div(u^n) u v
        if ( M_divBetaUv )
        {
            ASSERT(0!=0,"If really needed, please code it using ET");
            // OLD-STYLE assembly:
            // mass_divw ( -0.5 * M_oseenData->density(),
            //           M_uLoc,
            //           M_elementMatrixStiff,
            //           M_velocityFESpace.fe(), 0, 0, numVelocityComponent );
        }
    }
    
    if(M_stabilization)
        computeStabilization(u_starRepeated, alpha);

    if ( alpha != 0. )
    {
        *M_matrixNoBC_block += (*M_velocityMatrixMass) * alpha;
    }
    
    *M_convectiveMatrix += *M_matrixStokes;
    M_convectiveMatrix->globalAssemble();

    *M_matrixNoBC_block += *M_convectiveMatrix;
    M_matrixNoBC_block->globalAssemble();
    MapEpetra fullMap ( M_velocityFESpace.map() + M_pressureFESpace.map() + M_fluxMap );
    M_matrixNoBC.reset(new matrix_Type(fullMap));
    *M_matrixNoBC += *M_matrixNoBC_block;
    M_matrixNoBC->globalAssemble();
    
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::computeStabilization ( const vector_Type& u_star, const Real& alpha )
{
	Real normInf;
	u_star.normInf ( &normInf );

	LifeChrono chrono;

	if ( normInf != 0. )
	{
		if( M_oseenData->stabilizationType() == "IP" /*&& ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) )*/ )
		{
			vector_Type u_starRepeated ( u_star, Repeated );
			M_Displayer.leaderPrint ( "  F-  Updating the IP stabilization terms ... " );
			chrono.start();

			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			M_ipStabilization.apply ( *M_matrixStabilization, u_starRepeated, false );
			M_matrixStabilization->globalAssemble();
			M_resetStabilization = false;
			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
 		}
		else if(M_oseenData->stabilizationType() == "SUPG")
		{

			M_Displayer.leaderPrint ( "  F-  Updating the SUPG stabilization terms ... " );
			chrono.start();

			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_supgStabilization->applySUPG_Matrix_semi_implicit(M_matrixStabilizationET, u_star, alpha);
			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();

			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_supgStabilization->applySUPG_RHS_semi_implicit(M_rhsStabilization, u_star, *M_velocityRhs);

			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

		}
        else if(M_oseenData->stabilizationType() == "VMSLES")
		{
            
			M_Displayer.leaderPrint ( "  F-  Updating the VMSLES stabilization terms ... " );
			chrono.start();
            
			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_VMSLESStabilization->applyVMSLES_Matrix_semi_implicit(M_matrixStabilizationET,
                                                                    u_star,
                                                                    alpha,
                                                                    *M_velocityPreviousTimestep,
                                                                    *M_pressurePreviousTimestep);
            
			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();
            
			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_VMSLESStabilization->applyVMSLES_RHS_semi_implicit(M_rhsStabilization,
                                                                 u_star,
                                                                 *M_velocityRhs,
                                                                 alpha,
                                                                 *M_velocityPreviousTimestep,
                                                                 *M_pressurePreviousTimestep);
            
			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
            
		}
	}
	else
	{
		if (M_oseenData->stabilizationType() == "IP")
		{
			M_Displayer.leaderPrint ( "  F-  Updating the IP stabilization terms ... " );
			chrono.start();

			if ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) )
			{
				vector_Type u_starUnique ( u_star, Unique );
				M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
				M_ipStabilization.apply ( *M_matrixStabilization, u_starUnique, false );
				M_resetStabilization = false;
				chrono.stop();
				M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
			}
			else
			{
				M_Displayer.leaderPrint ( "reusing\n" );
			}
		}
		else if(M_oseenData->stabilizationType() == "SUPG")
		{

			M_Displayer.leaderPrint ( "  F-  Updating the SUPG stabilization terms ... " );
			chrono.start();

			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_supgStabilization->applySUPG_Matrix_semi_implicit(M_matrixStabilizationET, u_star, alpha);

			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();

			// comment because if steady simulation this is not needed
			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_supgStabilization->applySUPG_RHS_semi_implicit(M_rhsStabilization, u_star, *M_velocityRhs);

			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

		}
        else if(M_oseenData->stabilizationType() == "VMSLES")
		{
            
			M_Displayer.leaderPrint ( "  F-  Updating the VMSLES stabilization terms ... " );
			chrono.start();
            
			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_VMSLESStabilization->applyVMSLES_Matrix_semi_implicit(M_matrixStabilizationET,
                                                                    u_star,
                                                                    alpha,
                                                                    *M_velocityPreviousTimestep,
                                                                    *M_pressurePreviousTimestep);
            
			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();
            
			// comment because if steady simulation this is not needed
			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_VMSLESStabilization->applyVMSLES_RHS_semi_implicit(M_rhsStabilization,
                                                                 u_star,
                                                                 *M_velocityRhs,
                                                                 alpha,
                                                                 *M_velocityPreviousTimestep,
                                                                 *M_pressurePreviousTimestep);
            
			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
            
		}
	}
}


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::updateStabilization ( matrix_Type& matrixFull )
{
    if ( M_stabilization )
    {
    	M_matrixStabilization->globalAssemble();
        matrixFull += *M_matrixStabilization;
    }

}

template <typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::updateSourceTerm ( source_Type const& source )
{
    vector_Type rhs ( M_velocityFESpace.map() + M_pressureFESpace.map() + M_fluxMap );

    VectorElemental M_elvec (M_velocityFESpace->fe().nbFEDof(), nDimensions);
    UInt nc = nDimensions;

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= M_velocityFESpace->mesh()->numVolumes(); ++i )
    {

        M_velocityFESpace->fe().updateFirstDerivQuadPt ( M_velocityFESpace->mesh()->volumeList ( i ) );

        M_elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
        {
            compute_vec ( source, M_elvec, M_velocityFESpace->fe(),  M_oseenData->dataTime()->time(), ic ); // compute local vector
            assembleVector ( *rhs, M_elvec, M_velocityFESpace->fe(), M_velocityFESpace->dof(), ic, ic * M_velocityFESpace->getDim() ); // assemble local vector into global one
        }
    }
    M_rightHandSideNoBC *= 0;
    M_rightHandSideNoBC += rhs;
}




template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::iterate ( bcHandler_Type& bcHandler )
{

    LifeChrono chrono;

    // matrix and vector assembling communication
    M_Displayer.leaderPrint ( "  F-  Updating the boundary conditions ...     " );

    // HERE I SHOULD APPLY THE STABILIZATION ON THE RHS

    chrono.start();

    M_matrixNoBC->globalAssemble();
    
    matrixPtr_Type matrixFull ( new matrix_Type ( M_velocityFESpace.map() + M_pressureFESpace.map() + M_fluxMap ) );

    updateStabilization ( *matrixFull );
    getFluidMatrix ( *matrixFull );

    vector_Type rightHandSideFull ( *M_rightHandSideNoBC );

    if(M_stabilization)
    {
        if(M_oseenData->stabilizationType() == "SUPG" || M_oseenData->stabilizationType() == "VMSLES" )
        {
            rightHandSideFull += *M_rhsStabilization;
        }
    }

    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );

    // boundary conditions update
    M_Displayer.leaderPrint ("  F-  Applying boundary conditions ...         ");

    chrono.start();
    
    applyBoundaryConditions ( *matrixFull, rightHandSideFull, bcHandler );

    matrixFull->globalAssemble();
    
    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

    // solving the system
    M_linearSolver->setMatrix ( *matrixFull );

    boost::shared_ptr<MatrixEpetra<Real> > staticCast = boost::static_pointer_cast<MatrixEpetra<Real> > (matrixFull);

    Int numIter = M_linearSolver->solveSystem ( rightHandSideFull, *M_solution, staticCast );

    // if the preconditioner has been rese the stab terms are to be updated
    if ( numIter < 0 || numIter > M_iterReuseStabilization )
    {
        resetStabilization();
    }

    if (M_stabilization && M_oseenData->stabilizationType() == "VMSLES" )
    {
        *M_velocityPreviousTimestep *= 0;
        *M_pressurePreviousTimestep *= 0;
        
        M_velocityPreviousTimestep->subset(*M_solution, M_velocityFESpace.map(), 0, 0);
        M_pressurePreviousTimestep->subset(*M_solution, M_pressureFESpace.map(), 0, 0);
    }
    
    *M_residual  = *M_rightHandSideNoBC;
    *M_residual -= (*M_matrixNoBC) * (*M_solution);

    //M_residual.spy("residual");
} // iterate()

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::h1normVelocity(Real& uh1error )
{
	boost::shared_ptr<uExactFunctor> uExactFct( new uExactFunctor );
	boost::shared_ptr<vector_Type> uExactVec(new vector_Type(M_velocityFESpace.map(),Unique));
	M_velocityFESpace.interpolate(RossEthierSteinmanUnsteadyDec::uexact, *uExactVec, 0.0);

	vector_Type uComputed ( M_fespaceUETA->map() , Unique );
	uComputed.subset( *M_solution );

	Real errorH1SquaredLocal( 0.0 );
	Real errorH1Squared( 0.0 );

	boost::shared_ptr<gradUExactFunctor> gradUExactFct( new gradUExactFunctor );

	{
		using namespace ExpressionAssembly;
		integrate (
				elements (M_fespaceUETA->mesh() ), // Mesh
				M_velocityFESpace.qr(), // QR
				dot ( ( eval ( gradUExactFct, value(M_oseenData->dataTime()->time()), X) + (-1) * grad( M_fespaceUETA , uComputed ) ) ,
					  ( eval ( gradUExactFct, value(M_oseenData->dataTime()->time()), X) + (-1) * grad( M_fespaceUETA , uComputed ) ) )
		)
		>> errorH1SquaredLocal;

	}

	M_fespaceUETA->mesh()->comm()->Barrier();
	M_fespaceUETA->mesh()->comm()->SumAll (&errorH1SquaredLocal, &errorH1Squared, 1);
	uh1error= std::sqrt(errorH1Squared);
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::reduceSolution ( Vector& velocityVector, Vector& pressureVector )
{
    vector_Type solution ( *M_solution, 0 );

    if ( false /*S_verbose*/ )
    {
        for ( UInt iDof = 0; iDof < M_velocityFESpace.fieldDim() * dimVelocity(); ++iDof )
        {
            velocityVector[ iDof ] = solution[ iDof ];
        }

        for ( UInt iDof = 0; iDof < dimPressure(); ++iDof )
        {
            pressureVector[ iDof ] = solution[ iDof + M_velocityFESpace.fieldDim() * dimVelocity() ];
        }
    }

}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::reduceResidual ( Vector& residualVector )
{
    vector_Type residual ( *M_residual, 0 );

    if ( false /*S_verbose*/ )
    {
        for ( UInt iDof = 0; iDof < M_velocityFESpace.fieldDim() * dimVelocity(); ++iDof )
        {
            residualVector[ iDof ] = residual[ iDof ];
        }

    }
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::getFluidMatrix ( matrix_Type& matrixFull )
{
    M_matrixNoBC->globalAssemble();
    matrixFull += *M_matrixNoBC;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::postProcessingSetArea()
{
    M_postProcessing->set_area();
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::postProcessingSetNormal()
{
    M_postProcessing->set_normal();
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::postProcessingSetPhi()
{
    M_postProcessing->set_phi();
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::flux ( const markerID_Type& flag )
{
    return flux ( flag, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::flux ( const markerID_Type& flag,
                                          const vector_Type& solution )
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    return M_postProcessing->flux ( velocity, flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::kineticNormalStress ( const markerID_Type& flag )
{
    return kineticNormalStress ( flag, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::kineticNormalStress ( const markerID_Type& flag,
                                                         const vector_Type& solution )
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    return M_postProcessing->kineticNormalStress ( velocity, M_oseenData->density(), flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::area ( const markerID_Type& flag )
{
    return M_postProcessing->measure ( flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Vector
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::normal ( const markerID_Type& flag )
{
    return M_postProcessing->normal ( flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Vector
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::geometricCenter ( const markerID_Type& flag )
{
    return M_postProcessing->geometricCenter ( flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::pressure ( const markerID_Type& flag )
{
    return pressure ( flag, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::pressure (const markerID_Type& flag,
                                             const vector_Type& solution)
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type pressure ( this->M_pressureFESpace.map(), Repeated );
    pressure.subset ( velocityAndPressure,
                      this->M_velocityFESpace.dim() *this->M_velocityFESpace.fieldDim() );

    // third argument is 1, to use the pressure finite element space (see PostProcessingBoundary docs)
    return M_postProcessing->average ( pressure, flag, 1 ) [0];
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::meanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return meanNormalStress ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::meanNormalStress (const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution )
{
    if ( bcHandler.findBCWithFlag ( flag ).type() == Flux )
    {
        return -lagrangeMultiplier ( flag, bcHandler, solution );
    }
    else
    {
#ifdef HAVE_LIFEV_DEBUG
        M_Displayer.leaderPrint ( " !!! WARNING - OseenSolver::meanNormalStress( flag, bcHandler, solution) is returning an approximation \n" );
#endif
        return -pressure ( flag, solution ); // TODO: This is an approximation of the mean normal stress as the pressure.
        // A proper method should be coded in the PostprocessingBoundary class
        // to compute the exact mean normal stress.
    }
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::meanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return meanTotalNormalStress ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::meanTotalNormalStress (const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution )
{
    return meanNormalStress ( flag, bcHandler, solution ) - kineticNormalStress ( flag, solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::lagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return lagrangeMultiplier ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::lagrangeMultiplier ( const markerID_Type&  flag,
                                                        bcHandler_Type& bcHandler,
                                                        const vector_Type& solution )
{
    // Create a list of Flux bcName_Type ??
    std::vector< bcName_Type > fluxBCVector = bcHandler.findAllBCWithType ( Flux );
    bcName_Type fluxbcName_Type = bcHandler.findBCWithFlag ( flag ).name();

    // Create a Repeated vector for looking to the lambda
    vector_Type velocityPressureLambda ( solution, Repeated );

    // Find the index associated to the correct Lagrange multiplier
    for ( UInt lmIndex = 0; lmIndex < static_cast <UInt> ( fluxBCVector.size() ); ++lmIndex )
        if ( fluxbcName_Type.compare ( fluxBCVector[ lmIndex ] ) == 0 )
            return velocityPressureLambda[M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof()
                                          + M_pressureFESpace.dof().numTotalDof() + lmIndex];

    // If lmIndex has not been found a warning message is printed
    M_Displayer.leaderPrint (  "!!! Warning - Lagrange multiplier for Flux BC not found!\n" );
    return 0;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
VectorSmall<2> 
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::computeDrag ( const markerID_Type&  flag,
                                                                               bcHandler_Type& bcHandlerDrag,
                                                                               bcHandler_Type& bcHandlerLift,
                                                                               const Real& velocityInfty,
                                                                               const Real& Area)
{
    // 1) flag che dice in quale componente calcolare la resistenza (0 in x, 1 in y, 2 in z)
    // 2) creare vettore pieno di zeri lungo quanto il vettore soluzione -> vettore soluzione u organizzato come ( ux | uy | uz | p)
    // 3) mettere degli 1 nei nodi corrispondenti al bordo nel blocco del vettore, ovvero:
    //      se la resistenza é in direzione x: vettore da creare ( 1suBody | 0 | 0 | 0 )
    // 4) prodotto scalare del vettore del residuo con questo vettore definito al punto 3)
    
    
    vector_Type onesOnBodyDrag(M_localMap, Unique);
    onesOnBodyDrag *= 0;
    
    vector_Type onesOnBodyLift(M_localMap, Unique);
    onesOnBodyLift *= 0;
    
    bcManageRhs ( onesOnBodyDrag, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHandlerDrag, M_velocityFESpace.feBd(), 1., 0.);
    bcManageRhs ( onesOnBodyLift, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHandlerLift, M_velocityFESpace.feBd(), 1., 0.);
    
    Real drag (0.0);
    Real lift (0.0);
    
    drag = M_residual->dot(onesOnBodyDrag);
    lift = M_residual->dot(onesOnBodyLift);
    
    drag /= (0.5*M_oseenData->density()*velocityInfty*velocityInfty*Area);
    lift /= (0.5*M_oseenData->density()*velocityInfty*velocityInfty*Area);
    
    M_Displayer.leaderPrint ( "  F-  Value of the drag:          ", drag );
    M_Displayer.leaderPrint ( "  F-  Value of the lift:          ", lift );
                             
    VectorSmall<2> Coefficients;
    Coefficients[0] = drag;
    Coefficients[1] = lift;
    
    return Coefficients;
}
    
template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::removeMean ( vector_Type& x )
{

    LifeChrono chrono;
    chrono.start();

    const UInt numVelocityComponent ( M_velocityFESpace.fieldDim() );
    const UInt velocityTotalDof ( M_velocityFESpace.dof().numTotalDof() );


    if ( M_pressureMatrixMass.get() == 0 )
    {
        M_pressureMatrixMass.reset ( new matrix_Type ( M_localMap ) );
    }

    for ( UInt iElement = 0; iElement < M_velocityFESpace.mesh()->numElements(); iElement++ )
    {
        chrono.start();
        // just to provide the id number in the assem_mat_mixed
        M_pressureFESpace.fe().update ( M_pressureFESpace.mesh()->element ( iElement ) );

        M_elementMatrixPreconditioner.zero();
        // mass
        chrono.start();
        mass ( 1, M_elementMatrixPreconditioner, M_pressureFESpace.fe(), 0, 0, M_velocityFESpace.fieldDim() );
        chrono.stop();

        chrono.start();
        assembleMatrix ( *M_pressureMatrixMass,
                         M_elementMatrixPreconditioner,
                         M_pressureFESpace.fe(),
                         M_pressureFESpace.fe(),
                         M_pressureFESpace.dof(),
                         M_pressureFESpace.dof(),
                         numVelocityComponent,
                         numVelocityComponent,
                         numVelocityComponent * velocityTotalDof,
                         numVelocityComponent * velocityTotalDof );
        chrono.stop();
    }

    M_pressureMatrixMass->globalAssemble();

    vector_Type ones ( *M_solution );
    ones = 1.0;

    Real mean;
    mean = ones * ( M_pressureMatrixMass * x );
    x += ( -mean );

    ASSERT ( std::fabs ( ones * ( M_pressureMatrixMass * x ) ) < 1e-9 , "after removeMean the mean pressure should be zero!");

    return mean;

} // removeMean()

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::applyBoundaryConditions ( matrix_Type&       matrix,
                                                             vector_Type&       rightHandSide,
                                                             bcHandler_Type& bcHandler )
{
    // M_rightHandSideFull = M_rightHandSideNoBC;

    // BC manage for the velocity
    if ( !bcHandler.bcUpdateDone() || M_recomputeMatrix )
    {
        bcHandler.bcUpdate ( *M_velocityFESpace.mesh(),
                             M_velocityFESpace.feBd(),
                             M_velocityFESpace.dof() );
    }

    // ignoring non-local entries, Otherwise they are summed up lately
    //vector_Type rightHandSideFull( rightHandSide, Repeated, Zero );
    // ignoring non-local entries, Otherwise they are summed up lately
    vector_Type rightHandSideFull ( rightHandSide, Unique );

    bcManage ( matrix, rightHandSideFull,
               *M_velocityFESpace.mesh(),
               M_velocityFESpace.dof(),
               bcHandler,
               M_velocityFESpace.feBd(),
               1.,
               M_oseenData->dataTime()->time() );

    rightHandSide = rightHandSideFull;

    if ( bcHandler.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize ( M_velocityFESpace.fieldDim() *dimVelocity(),
                             M_diagonalize,
                             rightHandSide,
                             0. );
    }

} // applyBoundaryCondition

// ===================================================
// Set Methods
// ===================================================


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::setupPostProc( )
{
    M_postProcessing.reset ( new PostProcessingBoundary<mesh_Type> ( M_velocityFESpace.mesh(),
                                                                     &M_velocityFESpace.feBd(),
                                                                     &M_velocityFESpace.dof(),
                                                                     &M_pressureFESpace.feBd(),
                                                                     &M_pressureFESpace.dof(),
                                                                     M_localMap ) );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::setTolMaxIteration ( const Real& tolerance, const Int& maxIteration )
{
    M_linearSolver->setTolerance ( tolerance );
    M_linearSolver->setMaxNumIterations ( maxIteration );
}

} //end namespace LifeV
