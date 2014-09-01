//@HEADER
/*!
    @file
    @brief Solver class for hyperbolic scalar equations multirate extension.
	@version 1.0

    @date 27-12-2013

    @author Luca Oldani <luca.oldani@mail.polimi.it>
    
    @contributor

    @mantainer Luca Oldani <luca.oldani@mail.polimi.it>
*/

#ifndef _MULTIRATESOLVER_H_
#define _MULTIRATESOLVER_H_ 1

#include <lifev/core/solver/HyperbolicSolver.hpp>
#include <lifev/core/solver/HyperbolicData.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/core/solver/HyperbolicData.hpp>
#include <lifev/core/fem/HyperbolicFlux.hpp>
#include <lifev/core/fem/HyperbolicNumericalFlux.hpp>

// LifeV namespace.
namespace LifeV
{

const Real ConstLevel=4;

//MultirateSolver Class
template < typename MeshType, typename NumericalFluxType , typename PhysicalFluxType >
class MultirateSolver: 
public HyperbolicSolver<MeshType,NumericalFluxType,PhysicalFluxType> 
{

public:

	//! @name Public Types
    	//@{

    	//! Typedef for mesh template.
    	typedef typename HyperbolicSolver<MeshType,NumericalFluxType,PhysicalFluxType>::mesh_Type mesh_Type;

	typedef boost::shared_ptr < mesh_Type > meshPtr_Type;
	
	//! Typedef for partitioner template
	typedef MeshPartitioner < mesh_Type > partitioner_Type;

	//! Shared pointer to partitioner
	typedef boost::shared_ptr < partitioner_Type > partitionerPtr_Type;	

	//! Typedef for numerical flux template.
    	typedef typename HyperbolicSolver<MeshType,NumericalFluxType,PhysicalFluxType>::userNumericalFlux_Type userNumericalFlux_Type;

    	//! Typedef for numerical flux template.
    	typedef typename HyperbolicSolver<MeshType,NumericalFluxType,PhysicalFluxType>::userPhysicalFlux_Type userPhysicalFlux_Type;

    	//! Self typedef.
    	typedef MultirateSolver < mesh_Type, userNumericalFlux_Type, userPhysicalFlux_Type > multirateSolver_Type;

	//! Map type.
    	typedef typename HyperbolicSolver<MeshType,NumericalFluxType,PhysicalFluxType>::map_Type map_Type;

	//! Scalar field.
    	//typedef typename HyperbolicSolver<MeshType,NumericalFluxType,PhysicalFluxType>::scalarField_Type scalarField_Type;

	//! Scalar field.
    	typedef FEScalarField < mesh_Type, map_Type > scalarField_Type;
	
	//! Shared pointer to a scalar field.
    	typedef boost::shared_ptr < scalarField_Type > scalarFieldPtr_Type;

	//! Distributed vector.
    	typedef VectorEpetra vector_Type;

    	//! Shared pointer to a distributed vector.
    	typedef boost::shared_ptr < vector_Type > vectorPtr_Type;
	
	//! Time advance scheme.
    	typedef TimeAdvanceBDF < vector_Type > timeAdvance_Type;

    	//! Shared pointer to a time advance scheme.
    	typedef boost::shared_ptr < timeAdvance_Type > timeAdvancePtr_Type;
	
	//! Shared pointer to vector of UInt	
	typedef boost::shared_ptr <std::vector < UInt > > stdVectorPtr_Type;	

	//! Typedef for graph type
	typedef std::vector <std::vector < UInt > > graph_Type;

	//! Shared pointer to a graph type	
	typedef boost::shared_ptr < graph_Type > graphPtr_Type;

	//! Vector of pointer to Mesh Type to contein al subMesh	
	//typedef std::vector<meshPtr_Type> partMesh_Type;
	
	//! Shared Pointer to partMesh_Type
	//typedef boost::shared_ptr<partMesh_Type> partMeshPtr_Type;
	//@}


	//! @name Set Methos
    	//@{

	//! Set the CFL field.
    	/*!
      	@param set pointer to CFL Field.
    	*/
    	void setTimeStepField ( const scalarFieldPtr_Type& CFLField, const MapEpetra& ghostMap)
    	{
        	M_timeStepField = CFLField;

		// Set the localTimeStepField.
        	M_localTimeStepField.reset ( new vector_Type ( ghostMap, Repeated ) );

    	}

	/*!
      	@param set pointer to Region Field.
    	*/
    	void setRegionField ( const scalarFieldPtr_Type& RegionField, const MapEpetra& ghostMap)
    	{
        	M_regionField = RegionField;

		M_regionFieldOld.reset ( new scalarField_Type ( this->M_fESpace, ghostMap ) ); 
	
		M_localRegionField.reset ( new vector_Type ( ghostMap, Repeated ) );
    	}
	
	//@ param set slowTimeStep
	void setSlowTimeStep (const Real& slowTimeStep)
	{
		M_slowTimeStep = slowTimeStep;
	}

	//@ param set slowTimeStep
	void setFastTimeStep (const Real& fastTimeStep)
	{
		M_fastTimeStep = fastTimeStep;
	}

	//@ param set setOnlyFast
	void setIsSlow (bool& ifSlow)
	{
		M_ifSlow=ifSlow;
	}	
	
	//@ param set setOnlyFast
	void setIsSlow (bool ifSlow)
	{
		M_ifSlow=ifSlow;
	}


	//@}

	
	//! @name Get Methos
    	//@{

	//! Get the TimeStep field.
      	// @param return scalarFiledPtr to TimeStepField
    	scalarFieldPtr_Type& getTimeStepField ()  const
    	{
        	// Get the pointer to TimeStepField.
        	return M_timeStepField;
    	}

	
      	//! @param return scalarFiledPtr to regionField
    	scalarFieldPtr_Type& getRegionField ()  const
    	{
        	// Get the pointer to regionField
        	return M_regionField;
    	}

	//! @param return scalarFiledPtr to regionFieldOld
    	scalarFieldPtr_Type& getRegionFieldOld ()  const
    	{
        	// Get the pointer to regionField
        	return M_regionFieldOld;
    	}

	//! @param return the number of level set in class	
	UInt getLevel() const
	{
		return M_level;
	}

	//! @param return true if there is some elements in slow region
	bool isSlow() const
	{
		return M_ifSlow;
	}
	
	//@}

	//! @name Constructors & Destructor
    	//@{

    	//! Empty constructor.
    	MultirateSolver ();

    	//! Virtual destructor.
    	virtual ~MultirateSolver () {};

	//@}

	//! @name Methods
    	//@{
	
	//! Compute 2 region	
	void computeRegionAccordingTimeStep (int level=ConstLevel);

	//! Compute the global time step according to CFL condition, or return the user defined time step.
    	Real computeTimeStep () const;

	//! solve one slow time step
	void solveOneSlowTimeStep ();

	//! solve one fast time step
	void solveOneFastTimeStep ();

	//! setup Time Steps
	void setupTimeSteps ()
	{
		M_fastTimeStep = this->M_data->dataTimePtr()->timeStep();		
		M_slowTimeStep = M_level*M_fastTimeStep; 
	};
	//@}

protected:
	
	//! @name Protected Methods
    	//@{

	//! Compute the time step according to CFL condition
    	Real computeTimeStepAccordingCFL () const;

	//! Expand fast region of one cell
	void expandFastRegion ();
	
	//! Update M_elementDomains from M_regionField
	void updateDomains();	
	//@}

	// Data of the problem.
    	//! @name Data of the problem
    	//@{ 

	//! @name # of fast iteration in a slow iteration
	Real M_level;

	//! @name TimeStep of slow region
	Real M_slowTimeStep;

	//! @name TimeStep of fast region
	Real M_fastTimeStep;

	//! @name Field of TimeStep condition
	scalarFieldPtr_Type M_timeStepField;

	//! @name Field of TimeStep condition
	vectorPtr_Type M_localTimeStepField;

	//! @name Field of multiregion 1 correspond to fast region and 0 correspond to slow region
	scalarFieldPtr_Type M_regionField;

	//! @name Old Field of multiregion 1 correspond to fast region and 0 correspond to slow region
	scalarFieldPtr_Type M_regionFieldOld;

	//! @name local Field of region subdivision
	vectorPtr_Type M_localRegionField;

	//! @name Flag M_ifSlow is true if there is a slow region
	bool M_ifSlow;
	
	//@}

};//end of MultirateSolver Class

// ===================================================
// Constructor
// ===================================================

template < typename MeshType, typename NumericalFluxType, typename PhysicalFluxType >
MultirateSolver < MeshType, NumericalFluxType, PhysicalFluxType >::
MultirateSolver ()
{
	timeAdvancePtr_Type M_timeAdvance(new timeAdvance_Type );
	//graphPtr_Type M_elementDomains(new graph_Type );
	
} //constructor


// ===================================================
// Public Methods
// ===================================================

template < typename MeshType, typename NumericalFluxType, typename PhysicalFluxType >
void
MultirateSolver < MeshType, NumericalFluxType, PhysicalFluxType >::
computeRegionAccordingTimeStep (int level)
{
	M_level=level;
	M_ifSlow=false;

	stdVectorPtr_Type M_elementSlowDomain(new std::vector <UInt > );
	stdVectorPtr_Type M_elementFastDomain(new std::vector <UInt > );
	
	// Total number of elements in the mesh.
    	const UInt meshNumberOfElements = this->M_fESpace->mesh()->numElements();
	
	M_localRegionField->zero();

	VectorElemental MyRegion ( this->M_fESpace->refFE().nbDof(), 1 );

	VectorElemental valueCurrent ( this->M_fESpace->refFE().nbDof(), 1 );

	// Loop on all the elements to set two regions.
    	for ( UInt iElem(0); iElem < meshNumberOfElements; ++iElem )
    	{
		MyRegion.zero();
		
		const typename mesh_Type::element_Type& element = this->M_fESpace->mesh()->element ( iElem );

		this->M_fESpace->fe().update( element, UPDATE_QUAD_NODES | UPDATE_WDET | UPDATE_PHI );
	
		extract_vec(M_timeStepField->getVector(),
				valueCurrent,
				this->M_fESpace->refFE(),
				this->M_fESpace->dof(),
				iElem,0 );
		
		if (valueCurrent[0] >= M_level * (this->M_data->dataTimePtr()->timeStep()))
		{
			//then we can put this element in slow region
			MyRegion[0]=1;
			
			M_ifSlow=true;
			
		}
		else
		{
			
			MyRegion[0]=0;
		}

		assembleVector ( *M_localRegionField ,
       			this->M_fESpace->fe().currentLocalId(),
                       	MyRegion,
                       	this->M_fESpace->refFE().nbDof(),
                       	this->M_fESpace->dof(), 0 );		
	}

	M_localRegionField->globalAssemble();

	M_regionField->setVector(*M_localRegionField);
	
	M_localRegionField->zero();	
	
	if (M_ifSlow==true)
		for (UInt i(0); i<M_level; i++)
			{
				expandFastRegion();
			}	
	
}

// Compute the time step according to global CFL condition.
template < typename MeshType, typename NumericalFluxType, typename PhysicalFluxType >
Real
MultirateSolver < MeshType, NumericalFluxType, PhysicalFluxType >::
computeTimeStep () const
{
    // Check if the CFL has to be computed, otherwise use the user defined one.
    if ( this->M_data->getComputeCFL() )
    {
        // If the CFL has to be computed once, do it the first time step.
        if ( this->M_data->getTimeStepFrozen() == 0 && this->M_data->dataTimePtr()->timeStepNumber() == 0 )
        {
            return computeTimeStepAccordingCFL();
        }

        // Check if the current time step the CFL has to be computed.
        if ( this->M_data->getTimeStepFrozen() != 0 &&
             this->M_data->dataTimePtr()->timeStepNumber() % this->M_data->getTimeStepFrozen() == 0. )
        {
            return computeTimeStepAccordingCFL();
        }
    }

    // The CFL has not to be computed, return the computed or the user defined time step.
    return this->M_data->dataTimePtr()->timeStep();

}// computeTimeStep


// Solve one fast time step of the hyperbolic problem.
template < typename MeshType, typename NumericalFluxType, typename PhysicalFluxType >
void
MultirateSolver < MeshType, NumericalFluxType, PhysicalFluxType >::
solveOneSlowTimeStep ()
{
	
	this->M_displayer->leaderPrint ("Solve Slow TimeStep\n");	


    // Check if the initial solution is setted or not.
    ASSERT ( this->M_initialSolutionFct.get(), "HyperbolicSolver : initial solution not set." );

    // Check if the solution is setted or not.
    ASSERT ( this->M_solutionField.get(), "HyperbolicSolver : solution not set." );

    // Clean the vector of fluxes.
    this->M_globalFlux->zero();

	// # of elements
	const UInt meshNumberOfElements = this->M_fESpace->mesh()->numElements();

    // Local flux for local computations.
    VectorElemental localFlux ( this->M_fESpace->refFE().nbDof(), 1 );

    VectorElemental regionCurrent ( this->M_fESpace->refFE().nbDof(), 1 );

	// pointer tovector slowElemenDomain where are stored element that belong to slow region
	//boost::shared_ptr <std::vector < UInt > > fastElementDomain (&(*M_elementDomains)[1]);
	
    // Loop on all the elements to perform the fluxes.
    //for ( UInt i(0); i < fastElementDomain->size(); ++i )
	for ( UInt iElem(0); iElem < meshNumberOfElements ; ++iElem)    
	{

        // Clean the local flux.
        localFlux.zero();

        // Element of current ID.
        const typename mesh_Type::element_Type& element = this->M_fESpace->mesh()->element ( iElem );

		extract_vec(M_regionField->getVector(),
				regionCurrent,
				this->M_fESpace->refFE(),
				this->M_fESpace->dof(),
				iElem,0 );

	if (regionCurrent[0]==1){
        	// Update the property of the current element.
		this->M_fESpace->fe().update( element, UPDATE_QUAD_NODES | UPDATE_WDET | UPDATE_PHI );
	
        	// Reconstruct step of the current element.
        	this->localReconstruct ( iElem );
	
        	// Evolve step of the current element.
        	this->localEvolve ( iElem, localFlux );
	
        	// Put the total flux of the current element in the global vector of fluxes.
        	assembleVector ( *this->M_globalFlux,
        	                 this->M_fESpace->fe().currentLocalId(),
        	                 localFlux,
        	                 this->M_fESpace->refFE().nbDof(),
        	                 this->M_fESpace->dof(), 0 );
	
        	// Average step of the current element.
        	this->localAverage ( iElem );
		}
	}

    // Assemble the global hybrid vector.
    this->M_globalFlux->globalAssemble ();

    // Update the value of the solution.
    this->M_solutionField->setVector( this->M_solutionFieldOld->getVector() -
                                M_slowTimeStep * (*this->M_globalFlux) );

    	// alternative: instead of modifying M_globalFlux.map, we can make a local copy with the correct map
    	// this is needed since M_uOld.map != M_globalFlux.map
	//    vector_Type fluxCopy ( M_uOld->map() );
	//    fluxCopy = *M_globalFlux;
	//    *M_u = *M_uOld - M_data.dataTime()->timeStep() * fluxCopy;

    // Clean the vector of fluxes.
    this->M_globalFlux->zero();

    // Update the solution at previous time step.
    this->M_solutionFieldOld->setVector(  this->M_solutionField->getVector() );

} // solveOneSlowStep

// Solve one fast time step of the hyperbolic problem.
template < typename MeshType, typename NumericalFluxType, typename PhysicalFluxType >
void
MultirateSolver < MeshType, NumericalFluxType, PhysicalFluxType >::
solveOneFastTimeStep ()
{
	this->M_displayer->leaderPrint ("Solve One Fast Time Step\n");	


    // Check if the initial solution is setted or not.
    ASSERT ( this->M_initialSolutionFct.get(), "HyperbolicSolver : initial solution not set." );

    // Check if the solution is setted or not.
    ASSERT ( this->M_solutionField.get(), "HyperbolicSolver : solution not set." );

    // Clean the vector of fluxes.
    this->M_globalFlux->zero();

	// # of elements
	const UInt meshNumberOfElements = this->M_fESpace->mesh()->numElements();

    // Local flux for local computations.
    VectorElemental localFlux ( this->M_fESpace->refFE().nbDof(), 1 );

	VectorElemental valueCurrent ( this->M_fESpace->refFE().nbDof(), 1 );

    // Loop on all the elements to perform the fluxes.
    //for ( UInt i(0); i < fastElementDomain->size(); ++i )
	for ( UInt iElem(0); iElem < meshNumberOfElements ; ++iElem)    
	{

        // Clean the local flux.
        localFlux.zero();

        // Element of current ID.
        const typename mesh_Type::element_Type& element = this->M_fESpace->mesh()->element ( iElem );

		extract_vec(M_regionField->getVector(),
				valueCurrent,
				this->M_fESpace->refFE(),
				this->M_fESpace->dof(),
				iElem,0 );

	if (valueCurrent[0]==0){

        	// Update the property of the current element.
		this->M_fESpace->fe().update( element, UPDATE_QUAD_NODES | UPDATE_WDET | UPDATE_PHI );

        	// Reconstruct step of the current element.
        	this->localReconstruct ( iElem );

        	// Evolve step of the current element.
        	this->localEvolve ( iElem, localFlux );

	        // Put the total flux of the current element in the global vector of fluxes.
	        assembleVector ( *this->M_globalFlux,
	                         this->M_fESpace->fe().currentLocalId(),
	                         localFlux,
	                         this->M_fESpace->refFE().nbDof(),
	                         this->M_fESpace->dof(), 0 );

	        // Average step of the current element.
	        this->localAverage ( iElem );
		}
	}

    // Assemble the global hybrid vector.
    this->M_globalFlux->globalAssemble ();

    // Update the value of the solution.
    this->M_solutionField->setVector( this->M_solutionFieldOld->getVector() -
                                M_fastTimeStep * (*this->M_globalFlux) );

    	// alternative: instead of modifying M_globalFlux.map, we can make a local copy with the correct map
    	// this is needed since M_uOld.map != M_globalFlux.map
	//    vector_Type fluxCopy ( M_uOld->map() );
	//    fluxCopy = *M_globalFlux;
	//    *M_u = *M_uOld - M_data.dataTime()->timeStep() * fluxCopy;

    // Clean the vector of fluxes.
    this->M_globalFlux->zero();

    // Update the solution at previous time step.
    this->M_solutionFieldOld->setVector(  this->M_solutionField->getVector() );

} // solveOneFastStep

// ===================================================
// Protected Methods
// ===================================================

// Compute the CFL condition
template < typename MeshType, typename NumericalFluxType, typename PhysicalFluxType >
Real
MultirateSolver < MeshType, NumericalFluxType, PhysicalFluxType >::
computeTimeStepAccordingCFL () const
{
    // Total number of elements in the mesh.
    const UInt meshNumberOfElements = this->M_fESpace->mesh()->numElements();

    // Number of local faces for the geometrical element.
    const UInt elementNumberFacets = mesh_Type::element_Type::S_numLocalFacets;

    // The local value for the CFL condition, without the time step.
    Real localCFL ( 0. ), localCFLOld ( localCFL - 1. ), localCFLregion ( 0. );
	
    //clean my vector of timestepfield
    M_localTimeStepField->zero();
	
    // Loop on all the elements to perform the fluxes.
    for ( UInt iElem(0); iElem < meshNumberOfElements; ++iElem )
    {
	// Element of current ID.
        const typename mesh_Type::element_Type& element = this->M_fESpace->mesh()->element ( iElem );
	
        // Update the property of the current element.
        this->M_fESpace->fe().update ( element, UPDATE_QUAD_NODES | UPDATE_WDET | UPDATE_PHI );

        // Volumetric measure of the current element.
        const Real K = this->M_fESpace->fe().measure();

        // Loop on the faces of the element iElem and compute the local contribution.
        for ( ID iFacet(0); iFacet < elementNumberFacets; ++iFacet )
        {
            // Id mapping.
            const UInt facetId ( this->M_fESpace->mesh()->localFacetId( iElem, iFacet ) );

            // Element of current ID.
            const typename mesh_Type::facet_Type& elementFacet = this->M_fESpace->mesh()->facet ( facetId );

            // Take one element facing on the current facet.
            UInt elemCurrent = this->M_fESpace->mesh()->faceElement( elementFacet, 0 );

            // Take the other element facing on the facet.
            UInt elemNeighbor = this->M_fESpace->mesh()->faceElement( elementFacet, 1 );

            // Assign to the elemCurrent iElem and the other ID to elemNeighor.
            if ( elemCurrent != iElem )
            {
                std::swap ( elemCurrent, elemNeighbor );
            }

            // Update the measure and quadrature points of current facet.
            this->M_fESpace->feBd().update (elementFacet, UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES );
            
		// Solution in the left element.
            VectorElemental valueCurrent  ( this->M_fESpace->refFE().nbDof(), 1 );

            // Solution in the right element.
            VectorElemental valueNeighbor ( this->M_fESpace->refFE().nbDof(), 1 );

            // Extract the solution in the current element, now is the elemCurrent.
            extract_vec ( this->M_solutionFieldOld->getVector(), valueCurrent, this->M_fESpace->refFE(),
                          this->M_fESpace->dof(), elemCurrent , 0 );

            // Flag type of the current face.
            const flag_Type facetFlag = elementFacet.flag();

            if ( !Flag::testOneSet ( facetFlag, EntityFlags::PHYSICAL_BOUNDARY | EntityFlags::SUBDOMAIN_INTERFACE ) )
            {
                /* The current element is not a boundary element.
                   Extract the solution in the right element. */
                extract_vec ( this->M_solutionFieldOld->getVector(), valueNeighbor, this->M_fESpace->refFE(),
                              this->M_fESpace->dof(), elemNeighbor , 0 );
            }
            else if ( Flag::testOneSet ( facetFlag, EntityFlags::SUBDOMAIN_INTERFACE ) )
            {
                valueNeighbor[ 0 ] = this->M_solutionFieldOld->getVector()[ elemNeighbor ];
            }
            else
            {
                // The current element is a boundary face.
                valueNeighbor = this->inflow ( valueCurrent, iElem, iFacet, facetId );

                // The normal is oriented outward
                elemNeighbor = elemCurrent;
            }

            // Area of the current face, for the 1d case e=1
            const Real e ( ( MeshType::S_geoDimensions == 1 )? 1.: this->M_fESpace->feBd().measure() );

            // Number of quadrature points.
            const UInt boundaryQuadPoints = this->M_fESpace->feBd().nbQuadPt ();

            // Normal vector.
            const Vector3D normal = this->M_fESpace->mesh()->element( iElem ).normal( iFacet );

            // Loop on all the quadrature points.
            for ( UInt ig(0); ig < boundaryQuadPoints; ++ig )
            {
                // Current quadrature point.
                Vector3D quadPoint;

                for ( UInt icoor(0); icoor < MeshType::S_geoDimensions; ++icoor )
                {
                    quadPoint ( icoor ) = this->M_fESpace->feBd().quadPt ( ig, icoor );
                }

                 // Compute the mass coefficient.
                const Real massValue = this->M_massFct->eval ( iElem, quadPoint );

                // Compute the local CFL without the time step.
                localCFL = e / ( K * massValue ) *
                           this->M_numericalFlux.physicalFlux().maxAbsFluxPrimeDotNormal ( iElem, quadPoint,
                                                                                         this->M_data->dataTimePtr()->time(),
                                                                                         valueCurrent[0], valueNeighbor[0],
                                                                                         normal );

                // Check the integrity of the CFL.
                ASSERT ( localCFL >= 0., "Wrong value in compute CFL " );

                // Select the maximum between the old CFL condition and the new CFL condition.
                if ( localCFL > localCFLOld  )
                {
                    localCFLOld = localCFL;
                }
                else
                {
                    localCFL = localCFLOld;
                }

            }

        }

	VectorElemental MyCFL ( M_timeStepField->getFESpacePtr()->refFE().nbDof(), 1 );

	MyCFL[0]=(this->M_data->getCFLRelaxParameter()/localCFLOld);
	
	assembleVector ( *M_localTimeStepField ,
                         this->M_fESpace->fe().currentLocalId(),
                         MyCFL,
                         this->M_fESpace->refFE().nbDof(),
                         this->M_fESpace->dof(), 0 );

	if ( localCFLregion < localCFLOld  )
        {
               	localCFLregion = localCFLOld;
        }

	localCFLOld = -1.;

   }

M_localTimeStepField->globalAssemble();

M_timeStepField->setVector( *M_localTimeStepField );
	
// Compute the time step according to CLF for the current process.
Real timeStepLocal[] = { this->M_data->getCFLRelaxParameter() / localCFLregion };
Real timeStepGlobal[] = { 0. };

// Compute the minimum of the computed time step for all the processes.
this->M_displayer->comm()->MinAll( timeStepLocal, timeStepGlobal, 1 );

this->M_displayer->leaderPrint ("TIME STEP GLOBAL", *timeStepGlobal, "\n");	

// Return the computed value.
return *timeStepGlobal;
	
} // computeTimeStepAccordingCFL

//Expand Fast Region
template < typename MeshType, typename NumericalFluxType, typename PhysicalFluxType >
void
MultirateSolver < MeshType, NumericalFluxType, PhysicalFluxType >::
expandFastRegion (){

	// Total number of elements in the mesh.
    	const UInt meshNumberOfElements = this->M_fESpace->mesh()->numElements();
		
	// Number of local faces for the geometrical element.
    	const UInt elementNumberFacets = mesh_Type::element_Type::S_numLocalFacets;

	// Construct a ghost copy of M_regionField called M_regionFieldOld
	M_regionFieldOld->setVector(  M_regionField->getVector() );

	// Local flux for local computations.
    	VectorElemental OldRegionCurrent ( this->M_fESpace->refFE().nbDof(), 1 );

	VectorElemental OldRegionNeighbor ( this->M_fESpace->refFE().nbDof(), 1 );

	VectorElemental RegionCurrent ( this->M_fESpace->refFE().nbDof(), 1 );

	// Loop on all element
    	for ( UInt iElem(0); iElem < meshNumberOfElements; ++iElem )
    	{
		extract_vec ( M_regionField->getVector(), 
				OldRegionCurrent, 
				this->M_fESpace->refFE(),
                          	this->M_fESpace->dof(), 
				iElem , 0 );

		// Element of current ID.
	        const typename mesh_Type::element_Type& element = this->M_fESpace->mesh()->element ( iElem );
	
	        // Update the property of the current element.
	        this->M_fESpace->fe().update ( element, UPDATE_QUAD_NODES | UPDATE_WDET | UPDATE_PHI );
		
		if (OldRegionCurrent[0]==1){

			UInt iFacet(0);        	
			RegionCurrent[0]=1;

			while ( iFacet < elementNumberFacets && RegionCurrent[0] == 1)
			{
				// Id mapping.
        	    		const UInt facetId ( this->M_fESpace->mesh()->localFacetId( iElem, iFacet ) );
				
				// Element of current ID.
        	    		typename mesh_Type::facet_Type& elementFacet = this->M_fESpace->mesh()->facet ( facetId );

				// Take one element facing on the current facet.
        	    		UInt elemCurrent = elementFacet.firstAdjacentElementIdentity();

				// Take the other element facing on the current facet.
        	    		UInt elemNeighbor = elementFacet.secondAdjacentElementIdentity();

				// Assign to the elemCurrent iElem and the other ID to elemNeighor.
        	    		if ( elemCurrent != iElem )
        	    		{
        	        		std::swap ( elemCurrent, elemNeighbor );
        	    		}
		
				this->M_fESpace->feBd().update(elementFacet, UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES);

				// Flag type of the current face.
        	    		const flag_Type facetFlag = elementFacet.flag();

				// Control that current face is not a boundary face
				if ( !Flag::testOneSet ( facetFlag, EntityFlags::PHYSICAL_BOUNDARY | EntityFlags::SUBDOMAIN_INTERFACE ) )
        	    		{        	    								
				
		        		//The current element is not a boundary element.
        	           		//	Extract the regionField in the right element. 
        	        		extract_vec ( M_regionField->getVector(), 
							OldRegionNeighbor, 
							this->M_fESpace->refFE(),
        	                      			this->M_fESpace->dof(), 
							elemNeighbor , 0 );

					elemNeighbor=this->M_fESpace->mesh()->element (elemNeighbor).id();        	    		

				}
        	    		else if ( Flag::testOneSet ( facetFlag, EntityFlags::SUBDOMAIN_INTERFACE ) )
        	    		{
					OldRegionNeighbor[0]=M_regionFieldOld->getVector()[elemNeighbor];
				}

				if ( OldRegionNeighbor[0]== 0 )
				{
					
					if (elemCurrent==iElem && elemNeighbor<=meshNumberOfElements*2)					
					{
						RegionCurrent[0] = 0;
					}
				}
				else
				{
						RegionCurrent[0] = 1;
				}
				
				++iFacet;		
				
			}
		
		}
		else
		{
			RegionCurrent[0]=0;
		}

		assembleVector ( *M_localRegionField ,
       			this->M_fESpace->fe().currentLocalId(),
                       	RegionCurrent,
                       	this->M_fESpace->refFE().nbDof(),
                       	this->M_fESpace->dof(), 0 );

	}

	M_localRegionField->globalAssemble();

	//update real regions
	M_regionField->setVector(*M_localRegionField);

	M_localRegionField->zero();

}//expand Region

}//end of LifeV namespace

#endif /*_MULTIRATESOLVER_H_ */
