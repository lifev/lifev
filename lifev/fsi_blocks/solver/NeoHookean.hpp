#ifndef NEOHOOKEAN_H
#define NEOHOOKEAN_H 1

/*
 *  author: DAVIDE FORTI, davide.forti@epfl.ch
 *  neohookean structure
 *
 */

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/fem/BCManage.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/fsi_blocks/solver/AssemblyElementalStructure.hpp>

#include <lifev/fsi_blocks/solver/ExpressionDefinitions.hpp>

namespace LifeV
{

class NeoHookean
{

	typedef Epetra_Comm comm_Type;

	typedef boost::shared_ptr< comm_Type > commPtr_Type;

    typedef VectorEpetra vector_Type;

    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;

    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef boost::shared_ptr<BCHandler> bcPtr_Type;

    typedef RegionMesh<LinearTetra> mesh_Type;

    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

    typedef MapEpetra map_Type;

    typedef ETFESpace< mesh_Type, map_Type, 3, 3 > ETFESpace_displacement;

    typedef ExpressionDefinitions::deformationGradient_Type tensorF_Type;

    typedef ExpressionDefinitions::determinantTensorF_Type determinantF_Type;

    typedef ExpressionDefinitions::rightCauchyGreenTensor_Type tensorC_Type;

    typedef ExpressionDefinitions::minusTransposedTensor_Type minusT_Type;

    typedef ExpressionDefinitions::traceTensor_Type traceTensor_Type;

public:

    // empty constructor
    NeoHookean( const commPtr_Type& communicator );

    // empty destructor
    ~NeoHookean();

    void setCoefficients ( const Real density, const Real young, const Real poisson);

    void setup( const meshPtr_Type& mesh, const std::string dOrder );

    void evaluate_residual( const vectorPtr_Type& solution, const Real& coefficient, const vectorPtr_Type& csi, vectorPtr_Type& residual );

    void update_jacobian(const vectorPtr_Type& solution, const Real& coefficient, matrixPtr_Type& jacobian );

    const boost::shared_ptr<FESpace<mesh_Type, map_Type> >& fespace ( ) const { return M_displacementFESpace; };

    const boost::shared_ptr<ETFESpace_displacement >& et_fespace ( ) const { return M_displacementFESpace_ETA; };

private:

    // communicator
    commPtr_Type M_comm;

    // Coefficients for the structure
    Real M_density;

    // Lame coefficients computed on M_young and M_poisson
    Real M_bulk;
    Real M_mu;
    Real M_offset;

    matrixSmall_Type M_identity;

    // order FE space displacement
    std::string M_dOrder;

    // FE space
    boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_displacementFESpace;

    // ET FE Space
    boost::shared_ptr<ETFESpace_displacement > M_displacementFESpace_ETA;


    // jacobian matrix
    matrixPtr_Type M_jacobian;

};

} // end namespace LifeV

#endif
