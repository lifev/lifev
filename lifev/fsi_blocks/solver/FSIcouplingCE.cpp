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
FSIcouplingCE::setUp ( const Real& timeStep, const Real& interfaceDofs, const Real& beta, const Real& gamma,
					   const mapPtr_Type& interfaceMap, const FESpacePtr_Type& fluidVelocityFESpace,
					   const FESpacePtr_Type& structureDisplacementFESpace, const vectorPtr_Type& numerationInterface)
{
	M_timeStep = timeStep;
	M_interface = interfaceDofs; // Number of interface dofs per component
	M_fluidVelocityFESpace = fluidVelocityFESpace;
	M_structureDisplacementFESpace = structureDisplacementFESpace;
	M_interfaceMap = interfaceMap;
	M_numerationInterface = numerationInterface;
	M_beta = beta;
	M_gamma = gamma;
}

void
FSIcouplingCE::setUp ( const Real& timeStep, const Real& interfaceDofs, const Real& coefficientBDF,
	     	 	 	   const mapPtr_Type& interfaceMap, const FESpacePtr_Type& fluidVelocityFESpace,
	     	 	 	   const FESpacePtr_Type& structureDisplacementFESpace, const vectorPtr_Type& numerationInterface )
{
	M_timeStep = timeStep;
	M_interface = interfaceDofs; // Number of interface dofs per component
	M_fluidVelocityFESpace = fluidVelocityFESpace;
	M_structureDisplacementFESpace = structureDisplacementFESpace;
	M_interfaceMap = interfaceMap;
	M_numerationInterface = numerationInterface;
	M_coefficientBDF = coefficientBDF;
}

void
FSIcouplingCE::buildBlocks ( std::map<ID, ID> const& localDofMap, const bool& lambda_num_structure, bool useBDF )
{
	// if lambda_num_structure = true, lambda follows the numeration of the structure
	// if lambda_num_structure = false, lambda follows the numeration of the fluid

	Real value = 1.0;
	std::map<ID, ID>::const_iterator ITrow;

	M_lambdaToFluidMomentum.reset ( new MatrixEpetra<Real> ( M_fluidVelocityFESpace->map() ) );

	if ( lambda_num_structure )
	{
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
	}
	else
	{
		for (UInt dim = 0; dim < 3; ++dim)
		{
			for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
			{
				if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->first) ) >= 0 )
				{
					M_lambdaToFluidMomentum->addToCoefficient( ITrow->first + dim * M_fluidVelocityFESpace->dof().numTotalDof(),
															   (int) (*M_numerationInterface)[ITrow->first] + dim * M_interface,
															   value );
				}
			}
		}
	}

	M_lambdaToFluidMomentum->globalAssemble(M_interfaceMap, M_fluidVelocityFESpace->mapPtr());

	// ----------------------------------------------------------------------------------------------------------------------------

	value = -1.0;
	M_lambdaToStructureMomentum.reset ( new MatrixEpetra<Real> ( M_structureDisplacementFESpace->map() ) );

	if ( lambda_num_structure )
	{
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
	}
	else
	{
		for (UInt dim = 0; dim < 3; ++dim)
		{
			for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
			{
				if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->first) ) >= 0 )
				{
					M_lambdaToStructureMomentum->addToCoefficient( ITrow->second + dim * M_structureDisplacementFESpace->dof().numTotalDof(),
																   (int) (*M_numerationInterface)[ITrow->first] + dim * M_interface,
																   value );
				}
			}
		}
	}

	M_lambdaToStructureMomentum->globalAssemble(M_interfaceMap, M_structureDisplacementFESpace->mapPtr());

	// ----------------------------------------------------------------------------------------------------------------------------

	value = 1.0;
	M_fluidVelocityToLambda.reset ( new MatrixEpetra<Real> ( *M_interfaceMap ) );

	if ( lambda_num_structure )
	{
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
	}
	else
	{
		for (UInt dim = 0; dim < 3; ++dim)
		{
			for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
			{
				if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->first) ) >= 0 )
				{
					M_fluidVelocityToLambda->addToCoefficient( (int) (*M_numerationInterface)[ITrow->first] + dim * M_interface,
															   ITrow->first + dim * M_fluidVelocityFESpace->dof().numTotalDof(),
															   value );
				}
			}
		}
	}

	M_fluidVelocityToLambda->globalAssemble(M_fluidVelocityFESpace->mapPtr(), M_interfaceMap);

	// ----------------------------------------------------------------------------------------------------------------------------

	if ( useBDF )
	{
		value = ( (-1.0) / ( M_timeStep ) * M_coefficientBDF );
	}
	else
	{
		value = ( (-1.0 * M_gamma) / ( M_timeStep * M_beta ) );
	}

	M_structureDisplacementToLambda.reset ( new MatrixEpetra<Real> ( *M_interfaceMap ) );

	if ( lambda_num_structure )
	{
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
	}
	else
	{
		for (UInt dim = 0; dim < 3; ++dim)
		{
			for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
			{
				if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->first) ) >= 0 )
				{
					M_structureDisplacementToLambda->addToCoefficient( (int) (*M_numerationInterface)[ITrow->first] + dim * M_interface,
																	   ITrow->second + dim * M_structureDisplacementFESpace->dof().numTotalDof(),
																	   value );
				}
			}
		}
	}

	M_structureDisplacementToLambda->globalAssemble(M_structureDisplacementFESpace->mapPtr(), M_interfaceMap);

	// ----------------------------------------------------------------------------------------------------------------------------

	value = -1.0;
	M_structureDisplacementToFluidDisplacement.reset( new MatrixEpetra<Real> ( M_fluidVelocityFESpace->map() ) );

	if ( lambda_num_structure )
	{
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
	}
	else
	{
		for (UInt dim = 0; dim < 3; ++dim)
		{
			for (ITrow = localDofMap.begin(); ITrow != localDofMap.end(); ++ITrow)
			{
				if ( M_numerationInterface->map().map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->first) ) >= 0 )
				{
					M_structureDisplacementToFluidDisplacement->addToCoefficient( ITrow->first + dim * M_fluidVelocityFESpace->dof().numTotalDof(),
																				  ITrow->second + dim * M_structureDisplacementFESpace->dof().numTotalDof(),
																				  value );
				}
			}
		}
	}

	M_structureDisplacementToFluidDisplacement->globalAssemble(M_structureDisplacementFESpace->mapPtr(), M_fluidVelocityFESpace->mapPtr());
}



}
