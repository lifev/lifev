#include <lifev/fsi_blocks/solver/FSIcouplingCE.hpp>

namespace LifeV
{

FSIcouplingCE::FSIcouplingCE ( const commPtr_Type& communicator ) :
M_comm ( communicator )
{
}

FSIcouplingCE::~FSIcouplingCE ( )
{
}

void
FSIcouplingCE::setUp ( const Real& timeStep, const Real& interfaceDofs, const Real& coefficientFirstDerivative,
					   const mapPtr_Type& interfaceMap, const FESpacePtr_Type& fluidVelocityFESpace,
					   const FESpacePtr_Type& structureDisplacementFESpace, const vectorPtr_Type& numerationInterface)
{
	M_timeStep = timeStep;
	M_interface = interfaceDofs; // Number of interface dofs per component
	M_fluidVelocityFESpace = fluidVelocityFESpace;
	M_structureDisplacementFESpace = structureDisplacementFESpace;
	M_interfaceMap = interfaceMap;
	M_numerationInterface = numerationInterface;
	M_coefficientFirstDerivative = coefficientFirstDerivative;
}

void
FSIcouplingCE::buildBlocks ( std::map<ID, ID> const& localDofMap )
{
	Real value = -1.0;
	std::map<ID, ID>::const_iterator ITrow;

	M_lambdaToFluidMomentum.reset ( new MatrixEpetra<Real> ( M_fluidVelocityFESpace->map() ) );
	for (UInt dim = 0; dim < 3; ++dim)
	{
		for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
		{
			if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second) ) >= 0 )
			{
				M_lambdaToFluidMomentum->addToCoefficient( ITrow->first + dim * M_fluidVelocityFESpace->dof().numTotalDof(),
														   (int) (*M_numerationInterface)[ITrow->second] + dim * M_interface,
														   value );
			}
		}
	}
	M_lambdaToFluidMomentum->globalAssemble(M_interfaceMap, M_fluidVelocityFESpace->mapPtr());

	// ----------------------------------------------------------------------------------------------------------------------------

	value = 1.0;
	M_lambdaToStructureMomentum.reset ( new MatrixEpetra<Real> ( M_structureDisplacementFESpace->map() ) );
	for (UInt dim = 0; dim < 3; ++dim)
	{
		for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
		{
			if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second) ) >= 0 )
			{
				M_lambdaToStructureMomentum->addToCoefficient( ITrow->second + dim * M_structureDisplacementFESpace->dof().numTotalDof(),
													 	 	   (int) (*M_numerationInterface)[ITrow->second] + dim * M_interface,
													 	 	   value );
			}
		}
	}
	M_lambdaToStructureMomentum->globalAssemble(M_interfaceMap, M_structureDisplacementFESpace->mapPtr());

	// ----------------------------------------------------------------------------------------------------------------------------

	value = 1.0;
	M_fluidVelocityToLambda.reset ( new MatrixEpetra<Real> ( *M_interfaceMap ) );
	for (UInt dim = 0; dim < 3; ++dim)
	{
		for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
		{
			if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second) ) >= 0 )
			{
				M_fluidVelocityToLambda->addToCoefficient( (int) (*M_numerationInterface)[ITrow->second] + dim * M_interface,
						                                   ITrow->first + dim * M_fluidVelocityFESpace->dof().numTotalDof(),
						                                   value );
			}
		}
	}
	M_fluidVelocityToLambda->globalAssemble(M_fluidVelocityFESpace->mapPtr(), M_interfaceMap);

	// ----------------------------------------------------------------------------------------------------------------------------

	value = -1.0;
	value /= M_timeStep;
	value *= M_coefficientFirstDerivative;
	M_structureDisplacementToLambda.reset ( new MatrixEpetra<Real> ( *M_interfaceMap ) );
	for (UInt dim = 0; dim < 3; ++dim)
	{
		for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
		{
			if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second) ) >= 0 )
			{
				M_structureDisplacementToLambda->addToCoefficient( (int) (*M_numerationInterface)[ITrow->second] + dim * M_interface,
																   ITrow->second + dim * M_structureDisplacementFESpace->dof().numTotalDof(),
																   value );
			}
		}
	}
	M_structureDisplacementToLambda->globalAssemble(M_structureDisplacementFESpace->mapPtr(), M_interfaceMap);

	// ----------------------------------------------------------------------------------------------------------------------------

	value = -1.0;
	M_structureDisplacementToFluidDisplacement.reset( new MatrixEpetra<Real> ( M_fluidVelocityFESpace->map() ) );
	for (UInt dim = 0; dim < 3; ++dim)
	{
		for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
		{
			if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second) ) >= 0 )
			{
				M_structureDisplacementToFluidDisplacement->addToCoefficient( ITrow->first + dim * M_fluidVelocityFESpace->dof().numTotalDof(),
																			  ITrow->second + dim * M_structureDisplacementFESpace->dof().numTotalDof(),
																			  value );
			}
		}
	}
	M_structureDisplacementToFluidDisplacement->globalAssemble(M_structureDisplacementFESpace->mapPtr(), M_fluidVelocityFESpace->mapPtr());
}



}
