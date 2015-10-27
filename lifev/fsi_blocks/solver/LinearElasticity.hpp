#ifndef LINEARELASTICITY_H
#define LINEARELASTICITY_H 1

/*
 *  author: DAVIDE FORTI, davide.forti@epfl.ch
 *  Lightweighted class to Handle the time advancing scheme (based on Newmark).
 *
 */

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/fem/BCManage.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/eta/expression/Integrate.hpp>

namespace LifeV
{

class LinearElasticity
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

public:

    // empty constructor
    LinearElasticity( const commPtr_Type& communicator );

    // empty destructor
    ~LinearElasticity();

    void setCoefficients ( const Real density, const Real young, const Real poisson);

    void setup( const meshPtr_Type& mesh, const std::string dOrder );

    void assemble_matrices ( const Real timestep, const Real beta, bcPtr_Type & bc );

    matrixPtr_Type const& mass_matrix_no_bc ( ) const { return M_mass_no_bc; };

    matrixPtr_Type const& stiffness_matrix_no_bc ( ) const { return M_stiffness_no_bc; };

    matrixPtr_Type const& jacobian ( ) const { return M_jacobian; };

    const boost::shared_ptr<FESpace<mesh_Type, map_Type> >& fespace ( ) const { return M_displacementFESpace; };

    const boost::shared_ptr<ETFESpace_displacement >& et_fespace ( ) const { return M_displacementFESpace_ETA; };

private:

    // communicator
    commPtr_Type M_comm;

    // Coefficients for the structure
    Real M_density;
    Real M_young;
    Real M_poisson;

    // Lame coefficients computed on M_young and M_poisson
    Real M_lambda;
    Real M_mu;

    // order FE space displacement
    std::string M_dOrder;

    // FE space
    boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_displacementFESpace;

    // ET FE Space
    boost::shared_ptr<ETFESpace_displacement > M_displacementFESpace_ETA;

    // Matrices without bc applied
    matrixPtr_Type M_mass_no_bc;
    matrixPtr_Type M_stiffness_no_bc;

    // jacobian matrix
    matrixPtr_Type M_jacobian;

};

} // end namespace LifeV

#endif
