//@HEADER

#include <boost/numeric/ublas/matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesPreconditionerOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aPCDOperator.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesOperator.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/interpolation/RBFInterpolation.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

#ifndef _FSIAPPLYOPERATORNONCONFORMING_H_
#define _FSIAPPLYOPERATORNONCONFORMING_H_

namespace LifeV
{
namespace Operators
{
class FSIApplyOperatorNonConforming: public LinearOperatorAlgebra
{
public:
    //! @name Public Types
    //@{

    typedef  LinearOperatorAlgebra                     super;
    typedef  Epetra_CrsMatrix                          matrix_Type;
    typedef  boost::shared_ptr<matrix_Type>            matrixPtr_Type;
    typedef  MatrixEpetra<Real>                        matrixEpetra_Type;
    typedef  boost::shared_ptr<matrixEpetra_Type>      matrixEpetraPtr_Type;
    typedef  Epetra_Vector                             lumpedMatrix_Type;
    typedef  boost::shared_ptr<lumpedMatrix_Type>      lumpedMatrixPtr_Type;
    typedef  super::comm_Type                          comm_Type;
    typedef  super::commPtr_Type                       commPtr_Type;
    typedef  boost::shared_ptr<Teuchos::ParameterList> parameterListPtr_Type;
    typedef  MapEpetra                                 mapEpetra_Type;
    typedef  boost::shared_ptr<mapEpetra_Type>         mapEpetraPtr_Type;
    typedef  VectorEpetra                              VectorEpetra_Type;
    typedef  boost::shared_ptr<VectorEpetra_Type>      VectorEpetraPtr_Type;
    typedef  boost::shared_ptr<RBFInterpolation<RegionMesh<LinearTetra> > > interpolationPtr_Type;

	typedef LifeV::Preconditioner                  basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type>       basePrecPtr_Type;
	typedef LifeV::PreconditionerIfpack            prec_Type;
	typedef boost::shared_ptr<prec_Type>           precPtr_Type;
	typedef Teuchos::RCP< Teuchos::ParameterList > parameterListRCP_Type;

	typedef boost::shared_ptr<GetPot> datafilePtr_Type;

    //@}

    //! @name Constructors
    //@{
    //! Empty constructor
    FSIApplyOperatorNonConforming();
    //@}
    virtual ~FSIApplyOperatorNonConforming();

    //! @name SetUp
    //@{
    //! SetUp
    /*!
     */

    //! @name Set Methods
    //@{

    //! set the communicator
    void setComm ( const commPtr_Type & comm ) { M_comm = comm; };
    //! \warning Transpose of this operator is not supported
    int SetUseTranspose(bool UseTranspose){M_useTranspose = UseTranspose; return 0;}
    //! set the domain map
    void setDomainMap(const boost::shared_ptr<BlockEpetra_Map> & domainMap){M_operatorDomainMap = domainMap;}
    //! set the range map
    void setRangeMap(const boost::shared_ptr<BlockEpetra_Map> & rangeMap){M_operatorRangeMap = rangeMap;}
    //@}


    //! @name
    //@{
    //! Returns the Application of the Jacobian applied to \c X.
    int Apply(const vector_Type &X, vector_Type &Y) const;
    //! \warning No method \c Apply defined for this operator. It return an error code.
    int ApplyInverse(const vector_Type &/*X*/, vector_Type &/*Y*/) const {return -1;};
    //! \warning Infinity norm not defined for this operator
    double NormInf() const {return -1.0;}
    //@}

    // @name Attribute access functions
    //@{
    //! Return a character string describing the operator
    const char * Label() const {return M_label.c_str();}
    //! Return the current UseTranspose setting \warning Not Supported Yet.
    bool UseTranspose() const {return M_useTranspose;}
    //! Return false.
    bool HasNormInf() const {return false;}
    //! return a reference to the Epetra_Comm communicator associated with this operator
    const comm_Type & Comm() const {return *M_comm;}
    //! Returns the Epetra_Map object associated with the domain of this operator
    const map_Type & OperatorDomainMap() const {return *(M_monolithicMap->map(Unique));}
    //! Returns the Epetra_Map object associated with the range of this operator
    const map_Type & OperatorRangeMap() const {return *(M_monolithicMap->map(Unique));}
    //@}

    // @name Set
    //@{

    //! Set the monolithic map
    void setMonolithicMap(const mapEpetraPtr_Type& monolithicMap);

    //! Set the map of each component of the residual
    void setMaps( const mapEpetraPtr_Type& fluid_velocity_map,
    			  const mapEpetraPtr_Type& fluid_pressure_map,
    			  const mapEpetraPtr_Type& structure_displacement_map,
    			  const mapEpetraPtr_Type& lagrange_multipliers_map,
    			  const mapEpetraPtr_Type& ALE_map,
				  const mapEpetraPtr_Type& structure_interface_map);

    //! Set the shape derivatives
    void setUseShapeDerivatives( bool use ) { M_useShapeDerivatives = use; };

    //! Set the shape derivatives
    void setShapeDerivativesBlocks( const matrixEpetraPtr_Type & ShapeVelocity,
    								const matrixEpetraPtr_Type & ShapePressure);

    //! Copy the pointer of the interpolation objects
    void setInterpolants( interpolationPtr_Type fluidToStructure,
    		              interpolationPtr_Type structureToFluid);

    //! Set the blocks of the fluid Jacobian when stabilization is used
    void setFluidBlocks (   const matrixEpetraPtr_Type &  block00,
							const matrixEpetraPtr_Type &  block01,
							const matrixEpetraPtr_Type &  block10,
							const matrixEpetraPtr_Type &  block11);

    //! Set the blocks of the fluid Jacobian when stabilization is not used
    void setFluidBlocks ( 	const matrixEpetraPtr_Type &  block00,
							const matrixEpetraPtr_Type &  block01,
							const matrixEpetraPtr_Type &  block10);

    //! Set the block of the Jacobian of the structure
    void setStructureBlock ( const matrixEpetraPtr_Type &   structure ) { M_S = structure;};

    //! Set the block of the Jacobian of the ALE
    void setALEBlock ( const matrixEpetraPtr_Type & ale ) { M_G = ale;};

    //! Set the datafile needed by the solver of the interface mass
    void setDatafile( const GetPot& dataFile) { M_datafile = dataFile;};

    //! Set the timestep
    void setTimeStep( Real timeStep ) { M_timeStep = timeStep;};

    //! Set the blocks of the fluid Jacobian when stabilization is used
    void setInterfaceMassMatrices (  const matrixEpetraPtr_Type &  fluid_interface_mass,
    								 const matrixEpetraPtr_Type &  structure_interface_mass);

    //@}

    // @name Others
    //@{

    //! Apply the inverse of the fluid interface mass to a vector
    //  defined on the fluid interface
    void applyInverseInterfaceFluidMass(const VectorEpetraPtr_Type&X, VectorEpetraPtr_Type&Y) const;

    //@}

    //! Show information about the class
    void showMe();

    // @name Set the parameters
    //@{

    //@}

private:

    // @name Set the parameters
    //@{


    //@}

    //! Create the domain and the range maps
    void setMaps();

    boost::shared_ptr<BlockEpetra_Map> M_operatorDomainMap;
    //! Range Map
    boost::shared_ptr<BlockEpetra_Map> M_operatorRangeMap;

    //! Communicator
    commPtr_Type M_comm;

    bool M_useTranspose;

    // @name Maps
    //@{

    mapEpetraPtr_Type M_u_map;
    mapEpetraPtr_Type M_p_map;
    mapEpetraPtr_Type M_ds_map;
    mapEpetraPtr_Type M_lambda_map;
    mapEpetraPtr_Type M_ale_map;
    mapEpetraPtr_Type M_structure_interface_map;

    //@}

    // @name Offsets
    //@{

    UInt M_fluidVelocity;
    UInt M_fluid;
    UInt M_structure;
    UInt M_lambda;

    //@}

    // @name Parameter list of each block
    //@{

    //! Fluid blocks

    matrixEpetraPtr_Type M_F_00;
    matrixEpetraPtr_Type M_F_01;
    matrixEpetraPtr_Type M_F_10;
    matrixEpetraPtr_Type M_F_11;

    //! Interface masses

    matrixEpetraPtr_Type M_fluid_interface_mass;
    matrixEpetraPtr_Type M_structure_interface_mass;

    //! Shape derivatives

    matrixEpetraPtr_Type M_shapeVelocity;
    matrixEpetraPtr_Type M_shapePressure;

    //! Structure block

    matrixEpetraPtr_Type M_S;

    //! ALE block

    matrixEpetraPtr_Type M_G;

    //! Interpolants

    interpolationPtr_Type M_FluidToStructureInterpolant;
    interpolationPtr_Type M_StructureToFluidInterpolant;

    //@}

    mapEpetraPtr_Type M_monolithicMap;

    //! datafile
    GetPot M_datafile;

    //! Label
    const std::string M_label;

    //! Vectors needed for the Apply - input vectors associated to each part of the residual
    boost::shared_ptr<VectorEpetra_Type > M_X_velocity;
    boost::shared_ptr<VectorEpetra_Type > M_X_pressure;
    boost::shared_ptr<VectorEpetra_Type > M_X_displacement;
    boost::shared_ptr<VectorEpetra_Type > M_X_lambda;
    boost::shared_ptr<VectorEpetra_Type > M_X_geometry;

    //! Vectors needed for the Apply - output vectors associated to the application of the
    // Jacobian to each part of the residual
    boost::shared_ptr<VectorEpetra_Type > M_Y_velocity;
    boost::shared_ptr<VectorEpetra_Type > M_Y_pressure;
    boost::shared_ptr<VectorEpetra_Type > M_Y_displacement;
    boost::shared_ptr<VectorEpetra_Type > M_Y_lambda;
    boost::shared_ptr<VectorEpetra_Type > M_Y_geometry;

    //! If using the stabilization for the fluid
    bool M_useStabilization;

    //! If using the shape derivatives for the Jacobian
    bool M_useShapeDerivatives;

    Real M_timeStep;

};

} /* end namespace Operators */
} //end namespace
#endif
