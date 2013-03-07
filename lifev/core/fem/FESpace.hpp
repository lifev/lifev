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
    @brief This files contains the description and the implementation of the FESpace class.

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @date 00-06-2007

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @contributor Mauro Perego <mperego@fsu.edu>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */




#ifndef _FESPACE_H_
#define _FESPACE_H_

#include <iomanip>
#include <sstream>
#include <utility>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/fem/CurrentBoundaryFE.hpp>
#include <lifev/core/fem/DOF.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>


namespace LifeV
{


//! FESpace - Short description here please!
/*!
 *  @author Gilles Fourestey
 *
 *  Class representing the FE space, i.e. the reference FE and the geometric mapping.
 *
 */

template <typename MeshType, typename MapType>
class FESpace
{

public:

    //! @name Public Types
    //@{

    typedef boost::function < Real ( Real const&, Real const&, Real const&,
                                     Real const&, ID const& ) > function_Type;
    typedef MeshType                                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                     meshPtr_Type;
    typedef MapType                                          map_Type;
    typedef boost::shared_ptr<map_Type>                      mapPtr_Type;
    typedef typename map_Type::comm_ptrtype                  commPtr_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    /*!
      \param data_file GetPot data file
      \param refFE_u reference FE for the velocity
      \param refFE_p reference FE for the pressure
      \param Qr_u element quadrature rule for the velocity
      \param bdQr_u surface quadrature rule for the velocity
      \param Qr_p element quadrature rule for the pressure
      \param bdQr_p surface quadrature rule for the pressure
      \param BCh_fluid boundary conditions for the fluid
      \param ord_bdf order of the bdf time advancing scheme and incremental pressure approach (default: Backward Euler)
    */

    LIFEV_DEPRECATED ( FESpace ( MeshPartitioner<MeshType>& mesh,
                                 const ReferenceFE&         refFE,
                                 const QuadratureRule&      Qr,
                                 const QuadratureRule&      bdQr,
                                 const Int                  fDim,
                                 const commPtr_Type&        commptr
                               ) );

    LIFEV_DEPRECATED ( FESpace ( MeshPartitioner<MeshType>& mesh,
                                 const std::string&         space,
                                 const Int                  fDim,
                                 const commPtr_Type&        commptr
                               ) );

    FESpace ( meshPtr_Type            mesh,
              const ReferenceFE&      refFE,
              const QuadratureRule&   Qr,
              const QuadratureRule&   bdQr,
              const Int               fDim,
              const commPtr_Type&     commptr
            );

    FESpace ( meshPtr_Type           mesh,
              const std::string&     space,
              const Int              fDim,
              const commPtr_Type&    commptr
            );

    //! Do nothing destructor
    virtual ~FESpace() {}

    //@}


    //! @name Methods
    //@{

    //! Interpolate a given velocity function nodally onto a velocity vector
    template<typename vector_type>
    void interpolate ( const function_Type& fct, vector_type& vect, const Real time = 0.);

    template<typename vector_type>
    void interpolateBC ( BCHandler& BCh, vector_type& vect, const Real time);

    //! Interpolation method for FEFunctions
    /*!
      @param fEFunction Pointer to an FEFunction
      @param vector Interpolated function
      @param time Time in the interpolation
    */
    template < typename FEFunctionType, typename vector_Type >
    void interpolate ( const FEFunctionType* fEFunction, vector_Type& vector, const Real time = 0. );

    //! calculate L2 velocity error for given exact velocity function
    //! \param pexact the exact velocity as a function
    //! \param time the time
    //! \param relError Real* to store the relative error in

    // this computes vec = \int fct phi_i
    template<typename vector_type>
    void l2ScalarProduct ( const function_Type& fct,
                           vector_type&         vec,
                           const Real           t );

    template<typename vector_type>
    Real l20Error ( const function_Type& fexact,
                    const vector_type&   vec,
                    const Real           time,
                    Real*                relError = 0 );

    template<typename vector_type>
    Real l2Error ( const function_Type& fexact,
                   const vector_type&   vec,
                   const Real           time,
                   Real*                relError = 0 );

    template<typename function, typename vector_type>
    Real h1Error ( const function&     fexact,
                   const vector_type&  vec,
                   const Real          time,
                   Real*               relError = 0 );


    template<typename vector_type>
    Real l2Norm ( const vector_type& vec );

    template<typename function>
    Real l2NormFunction ( const function& f, const Real time = 0);

    template<typename vector_type>
    Real h1Norm ( const vector_type& vec );


    //! Method to computes the L2 error when using a weight function
    /*!
      The scope of this method is to compute \f$ \left( \int w (u_{exact} - u_h) \right)^{1/2} \f$.
      The usual L2 error norm can be retrieved by using \f$ w=1 \f$
     */
    template<typename vector_type>
    Real l2ErrorWeighted ( const function_Type& exactSolution,
                           const vector_type&   solution,
                           const function_Type& weight,
                           const Real           time);


    //! This method computes the interpolate value of a given FE function in a given point.
    /*!
     * The user of this function has to provide the element and the vector of the DOF values.
     * Given these informations and a point P, the function compute the value:
     * value = sum_{dof i} v(i) phi_i(P)
     *
     * Note that the point P can be outside of the element considered (this is NOT checked).
     *
     * Warning: this method has been only checked in 3D.
     *
     *  @param elementID  The ID of the element considered. The ID is the local ID
     * (not the global one, see in the MeshEntity class).
     *  @param solutionVector  The vector containing the values in all nodes (not only of
     * the considered element).
     *  @param pt  The point where to perform the interpolation. Note that pt must allow
     * for STL-type accessor [] and must have the size() method (typically a std::vector<Real>).
     *  @param component  The component for which the interpolation has to be performed
     * (usefull only for the vectorial FE, set to 0 if the FE is scalar).
     *  @return  The interpolated value in the point.
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename point_type, typename vector_type>
    Real feInterpolateValue (const ID& elementID, const vector_type& solutionVector,
                             const point_type& pt, const UInt& component = 0) const;


    //! This method computes the interpolate value of a given FE function in a given point.
    /*!
     * This method is the same as feInterpolateValue, but for the definition of the solution
     * vector. Here, the solution is given only for the degrees of freedom of the given
     * element. The parameter solutionVector is typically a std::vector<Real> containing
     * the values in the dofs.
     *
     * It is not possible to specify the component, the user has to provide directly the
     * values for the wanted component in the solutionVector.
     *
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename point_type, typename vector_type>
    Real feInterpolateValueLocal (const ID& elementID, const vector_type& solutionVector,
                                  const point_type& pt ) const;


    //! This method computes the interpolated gradient of a given FE function in a given point.
    /*!
     * This method is the same as feInterpolateValue, but it computes the gradient insted of
     * the value. Therefor, the desired element of the gradient has to be expressed using the
     * parameter gradientElement.
     *
     * For example, if gradientElement=i and component=j, then the results corresponds to the
     * value of partial u_j / partial x_i.
     *
     * Warning: This method has not been tested deeply, buggy behaviour is possible.
     *
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename point_type, typename vector_type>
    Real feInterpolateGradient (const ID& elementID, const vector_type& solutionVector, const point_type& pt,
                                const UInt& gradientElement, const UInt& component = 0 ) const;


    //! This method computes the interpolated gradient of a given FE function in a given point.
    /*!
     * This method is the same as feInterpolateGradient, but it works with local values only,
     * just as feInterpolateValueLocal.
     *
     * Warning: This method has not been tested deeply, buggy behaviour is possible.
     *
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template<typename point_type, typename vector_type>
    Real feInterpolateGradientLocal (const ID& elementID, const vector_type& solutionVector, const point_type& pt,
                                     const UInt& gradientElement ) const;


    //! This method enables to pass a solution from one FESpace to the present one.
    /*!
     * This method interpolates the values of a FE solution (given by the vector and the FESpace)
     * in the dofs of the present FESpace: we note \f$ v_i \f$ the component of the vector given in
     * argument and \f$ phi_i \f$ the basis functions of the FESpace given in argument. This defines a
     * function on the domain given by: \f$ f(x) = \sum_i v_i \phi_i(x) \f$.
     *
     * We search then the function \f$ g(x) \f$ belonging to the present FESpace such that \f$ g(x_j) = f(x_j) \f$
     * for all \f$ x_j \f$ dofs of the present FESpace. This function returns the vector  \f$ w_j \f$ such that
     * \f$ g(x) = sum_j w_j psi_j(x) \f$ where \f$ psi_j \f$ are the basis functions of the present space.
     *
     * Warning: It is NOT true in general that \f$ w_j = f(x_j) \f$ (it is for lagrangian FEs). For
     * example, if we want to pass the solution from the P_1 to the P_1 with bubbles, then
     * the value associated to the bubble function will be always 0 even if the value of f
     * is not 0 at the location of the dof of the bubble!
     *
     *  @param originalSpace  The space where the solution is defined originally.
     *  @param originalVector The vector of the solution in the original space.
     *  @param outputMapType  The map type (default: Unique) of the the returned vector.
     *  @return The vector in the current FESpace having map type outputMapType.
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */

    template <typename vector_type>
    vector_type feToFEInterpolate (const FESpace<mesh_Type, map_Type>& originalSpace,
                                   const vector_type& originalVector, const MapEpetraType& outputMapType = Unique) const;

    //! This method reconstruct a gradient of a solution in the present FE space.
    /*!
      The goal of this method is to build an approximation of the gradient of a given
      FE function in this FESpace. Typically, when one use P1 elements for approximating
      the solution of a given problem, the gradient is only piecewise constant. However, one
      could need continuous gradient. The solutions to this problem is either to use specific
      finite elements (like Hermite FE) or rely on a recovery procedure for the gradient.

      This method implements a recovery procedure that performs a local average with weights
      corresponding to the areas of the elements:
      \f[ Gr(P) = \frac{\sum_{P \in T} |T| \nabla u(P)}{\sum_{P \in T} |T|} \f]
      See Zienkiewicz and Zhu (1987) for more details.

      @Note Results might be very wrong if you are not using lagrangian FE for tetrahedra
     */
    template <typename vector_type>
    vector_type gradientRecovery (const vector_type& solution, const UInt& component) const;

    //! Reconstruction of the laplacian using gradientRecovery procedures.
    /*!
      This method simply uses the FESpace::gradientRecovery method several times so
      that one can get a continuous approximation of the laplacian of the given solution.

      @Note Results might be very wrong if you are not using lagrangian FE for tetrahedra
     */
    template <typename vector_type>
    vector_type laplacianRecovery (const vector_type& solution) const;

    //! Return the polynomial degree of the finite element used
    UInt polynomialDegree() const;

    //@}

    //! @name Set Methods
    //@{

    //! Method to set replace the quadrule
    /*!
      @param Qr The new quadrule to be used in the FESpace
     */
    void setQuadRule (const QuadratureRule& Qr);

    //@}


    //! @name Get Methods
    //@{

    //    const mesh_type& mesh()       {return M_mesh;}
    //    const mesh_type& mesh() const {return *M_mesh;}
    //    mesh_type& mesh() {return *M_mesh;}

    const meshPtr_Type  mesh()  const
    {
        return M_mesh;
    }
    meshPtr_Type        mesh()
    {
        return M_mesh;
    }

    //! Returns map
    const map_Type&     map()   const
    {
        return *M_map;
    }
    map_Type&           map()
    {
        return *M_map;
    }

    const mapPtr_Type& mapPtr() const
    {
        return M_map;
    }

    //! Returns the velocity dof
    const DOF&          dof()   const
    {
        return *M_dof;
    }
    DOF&                dof()
    {
        return *M_dof;
    }
    const boost::shared_ptr<DOF>& dofPtr() const
    {
        return M_dof;
    }

    //! Returns the current FE
    const CurrentFE&    fe()    const
    {
        return *M_fe;
    }
    CurrentFE&          fe()
    {
        return *M_fe;
    }

    //! Returns the current boundary FE
    CurrentBoundaryFE&      feBd()
    {
        return *M_feBd;
    }

    //! Returns the res FE
    const ReferenceFE&      refFE() const
    {
        return *M_refFE;
    }

    //! Returns the volumic quadratic rule
    const QuadratureRule&       qr()    const
    {
        return *M_Qr;
    }

    //! Returns the surfasic quadratic rule
    const QuadratureRule&       bdQr()  const
    {
        return *M_bdQr;
    }

    //! Returns FE space dimension
    const UInt&         dim()      const
    {
        return M_dim;
    }
    const UInt&         fieldDim() const
    {
        return M_fieldDim;
    }

    //@}



private:

    //! @name Private Methods
    //@{

    //! copy constructor
    FESpace ( const FESpace& fespace );

    //! Creates the map for interprocessor communication
    void  createMap (const commPtr_Type& commptr);

    //! Resets boundary data if necessary
    void resetBoundaryFE();

    //! Set space
    inline void setSpace ( const std::string& space, UInt dimension );


    //! This is a generic function called by feToFEInterpolate method.
    //! It allows to interpolate vectors between any two continuous and scalar FE spaces. It is used when other specialized functions are not provided
    /*!
     *  @param originalSpace  The space where the solution is defined originally.
     *  @param originalVector The vector of the solution in the original space (must be a Repeated vector).
     *  @return The Repeated vector in the current FESpace.
     *  @author Mauro Perego
     *  @date 27-06-2011
     */
    template<typename vector_type>
    vector_type interpolateGeneric ( const FESpace<mesh_Type, map_Type>& OriginalSpace,
                                     const vector_type& OriginalVector) const;

    //! This is a specialized function called by feToFEInterpolate method.
    //! It allows to interpolate P1bubble, P2 vectors into P1 vectors, P1 into P1bubble vectors and Q2 vectors into Q1 vectors.
    /*!
     *  @param originalSpace  The space where the solution is defined originally.
     *  @param originalVector The vector of the solution in the original space.
     *  @return The Unique vector in the current FESpace.
     *
     *  @author Mauro Perego
     *  @date 27-06-2011
     */

    template<typename vector_type>
    vector_type linearInterpolate (const FESpace<mesh_Type, map_Type>& originalSpace,
                                   const vector_type& originalVector) const;

    //! This is a specialized function called by feToFEInterpolate method.
    //! It allows to interpolate P1, P1bubble vectors into P2 vectors, and Q1 vectors into Q2 vectors.
    /*!
     *  @param originalSpace  The space where the solution is defined originally.
     *  @param originalVector The vector of the solution in the original space (must be a Repeated vector).
     *  @return The Repeated vector in the current FESpace.
     *  @author Samuel Quinodoz
     *  @date 13-01-2010
     */
    template<typename vector_type>
    vector_type P2Interpolate (const FESpace<mesh_Type, map_Type>& original_space,
                               const vector_type& original_vector) const;


    //! This is a specialized function called by FESpace::feToFEInterpolate method for RT0 to P0 interpolation
    /*!
     *  @author Alessio Fumagalli
     *  @date 13-01-2010
     */
    template<typename vector_type>
    vector_type RT0ToP0Interpolate (const FESpace<mesh_Type, map_Type>& original_space,
                                    const vector_type& original_vector) const;
    //@}




    //! Set space Map (useful for switch syntax with strings)
    enum spaceData {P1, P1_HIGH, P1Bubble, P2, P2_HIGH};
    std::map<std::string, spaceData>        M_spaceMap;

    //! reference to the mesh
    meshPtr_Type                            M_mesh;

    //! Reference FE for the velocity
    const ReferenceFE*                      M_refFE;

    //! Quadrature rule for volumic elementary computations
    const QuadratureRule*                   M_Qr;

    //! Quadrature rule for surface elementary computations
    const QuadratureRule*                   M_bdQr;

    //! dimension of the field variable ( scalar/vector field)
    UInt                                    M_fieldDim;

    //! A shared pointer to the DOF object
    boost::shared_ptr<DOF>                  M_dof;

    //! The number of total dofs
    UInt                                    M_dim;

    //! Current FE
    boost::shared_ptr<CurrentFE>            M_fe;
    boost::shared_ptr<CurrentBoundaryFE>    M_feBd;

    //! Map
    mapPtr_Type                             M_map;

};

// ===================================================
// Constructors & Destructor
// ===================================================

template <typename MeshType, typename MapType>
FESpace<MeshType, MapType>::
FESpace ( MeshPartitioner<MeshType>&  mesh,
          const ReferenceFE&          refFE,
          const QuadratureRule&       Qr,
          const QuadratureRule&       bdQr,
          const Int                   fDim,
          const commPtr_Type&         commptr
        ) :
    M_mesh          ( mesh.meshPartition() ),
    M_refFE         ( &refFE ),
    M_Qr            ( &Qr ),
    M_bdQr          ( &bdQr ),
    M_fieldDim      ( fDim ),
    M_dof           ( new DOF ( *M_mesh, *M_refFE ) ),
    M_dim           ( M_dof->numTotalDof() ),
    M_fe            ( new CurrentFE  ( *M_refFE, getGeometricMap ( *M_mesh ), *M_Qr ) ),
    M_feBd          ( ),
    M_map           ( new map_Type() )
{
    resetBoundaryFE();
    createMap (commptr);
}

template <typename MeshType, typename MapType>
FESpace<MeshType, MapType>::
FESpace ( MeshPartitioner<MeshType>& mesh,
          const std::string&         space,
          const Int                  fDim,
          const commPtr_Type&        commptr
        ) :
    M_mesh          ( mesh.meshPartition() ),
    M_fieldDim      ( fDim ),
    M_dof           ( ),
    M_fe            ( ),
    M_feBd          ( ),
    M_map           ( new map_Type() )
{
    // Set spaceMap
    M_spaceMap["P1"]        = P1;
    M_spaceMap["P1_HIGH"]   = P1_HIGH;
    M_spaceMap["P1Bubble"]  = P1Bubble;
    M_spaceMap["P2"]        = P2;
    M_spaceMap["P2_HIGH"]   = P2_HIGH;

    // Set space
    setSpace ( space, mesh_Type::S_geoDimensions );

    // Set other quantities
    M_dof.reset ( new DOF ( *M_mesh, *M_refFE ) );
    M_dim = M_dof->numTotalDof();
    M_fe.reset ( new CurrentFE ( *M_refFE, getGeometricMap ( *M_mesh ), *M_Qr ) );
    resetBoundaryFE();
    createMap (commptr);
}

template <typename MeshType, typename MapType>
FESpace<MeshType, MapType>::
FESpace ( meshPtr_Type           mesh,
          const ReferenceFE&     refFE,
          const QuadratureRule&  Qr,
          const QuadratureRule&  bdQr,
          const Int              fDim,
          const commPtr_Type&    commptr
        ) :
    M_mesh          ( mesh ),
    M_refFE         ( &refFE ),
    M_Qr            ( &Qr ),
    M_bdQr          ( &bdQr ),
    M_fieldDim      ( fDim ),
    M_dof           ( new DOF ( *M_mesh, *M_refFE ) ),
    M_dim           ( M_dof->numTotalDof() ),
    M_fe            ( new CurrentFE ( *M_refFE, getGeometricMap ( *M_mesh ), *M_Qr ) ),
    M_feBd          ( ),
    M_map           ( new map_Type() )
{
    createMap (commptr);
    resetBoundaryFE();
}

template <typename MeshType, typename MapType>
FESpace<MeshType, MapType>::
FESpace ( meshPtr_Type         mesh,
          const std::string&   space,
          const Int            fDim,
          const commPtr_Type&  commptr
        ) :
    M_mesh          ( mesh ),
    M_fieldDim      ( fDim ),
    M_dof           ( ),
    M_fe            ( ),
    M_feBd          ( ),
    M_map           ( new map_Type() )
{
    // Set spaceMap
    M_spaceMap["P1"]        = P1;
    M_spaceMap["P1_HIGH"]   = P1_HIGH;
    M_spaceMap["P1Bubble"]  = P1Bubble;
    M_spaceMap["P2"]        = P2;
    M_spaceMap["P2_HIGH"]   = P2_HIGH;

    // Set space
    setSpace ( space, mesh_Type::S_geoDimensions );

    // Set other quantities
    M_dof.reset ( new DOF ( *M_mesh, *M_refFE ) );
    M_dim = M_dof->numTotalDof();
    M_fe.reset ( new CurrentFE ( *M_refFE, getGeometricMap ( *M_mesh ), *M_Qr ) );
    resetBoundaryFE();
    createMap (commptr);
}

// ===================================================
// Methods
// ===================================================


template <typename MeshType, typename MapType>
template<typename vector_type>
void
FESpace<MeshType, MapType>::interpolate ( const function_Type& fct,
                                          vector_type&    vect,
                                          const Real      time)
{
    // First, we build a "quadrature" that consists in the nodes (0 weight)
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_refFE->shape() ), M_refFE->shape() );
    interpQuad.setPoints (M_refFE->refCoor(), std::vector<Real> (M_refFE->nbDof(), 0) );

    // Then, we define a currentFE with nodes on the reference nodes
    CurrentFE interpCFE (*M_refFE, getGeometricMap (*M_mesh ), interpQuad);

    // Some constants
    UInt totalNumberElements (M_mesh->numElements() );
    UInt numberLocalDof (M_dof->numLocalDof() );

    // Storage for the values
    std::vector<Real> nodalValues (numberLocalDof, 0);
    std::vector<Real> FEValues (numberLocalDof, 0);

    // Do the loop over the cells
    for (UInt iterElement (0); iterElement < totalNumberElements; ++iterElement)
    {
        // We update the CurrentFE so that we get the coordinates of the nodes
        interpCFE.update (M_mesh->element (iterElement), UPDATE_QUAD_NODES);

        // Loop over the dimension of the field
        for (UInt iDim (0); iDim < M_fieldDim; ++iDim)
        {
            // Loop over the degrees of freedom (= quadrature nodes)
            for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
            {
                // Store the nodal value
                nodalValues[iterDof] =  fct (time, interpCFE.quadNode (iterDof, 0), interpCFE.quadNode (iterDof, 1)
                                             , interpCFE.quadNode (iterDof, 2), iDim);
            }

            // Transform the nodal values in FE values
            FEValues = M_refFE->nodalToFEValues (nodalValues);

            // Then on the dimension of the FESpace (scalar field vs vectorial field)
            for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
            {
                // Find the ID of the considered DOF
                ID globalDofID (M_dof->localToGlobalMap (iterElement, iterDof) + iDim * M_dim);

                // Compute the value of the function and set it
                vect.setCoefficient (globalDofID, FEValues[iterDof]);

            }
        }
    }
}

template < typename MeshType, typename MapType>
template < typename FEFunctionType, typename vector_Type >
void FESpace<MeshType, MapType>::
interpolate ( const FEFunctionType* fEFunction, vector_Type& vector, const Real time )
{

    // First, we build a "quadrature" that consists in the nodes (0 weight)
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape ( shapeDimension (M_refFE->shape() ), M_refFE->shape() );
    interpQuad.setPoints ( M_refFE->refCoor(), std::vector<Real> (M_refFE->nbDof(), 0) );

    // Then, we define a currentFE with nodes on the reference nodes
    CurrentFE interpCFE ( *M_refFE, getGeometricMap ( *M_mesh ), interpQuad );

    // Some constants
    const UInt totalNumberElements ( M_mesh->numElements() );
    const UInt numberLocalDof ( M_dof->numLocalDof() );

    // Storage for the values
    typename FEFunctionType::point_Type point;
    std::vector<Real> nodalValues (numberLocalDof, 0);
    std::vector<Real> FEValues (numberLocalDof, 0);

    // Do the loop over the cells
    for (UInt iterElement (0); iterElement < totalNumberElements; ++iterElement)
    {
        // We update the CurrentFE so that we get the coordinates of the nodes
        interpCFE.update ( M_mesh->element (iterElement), UPDATE_QUAD_NODES );

        // Loop over the dimension of the field
        for (UInt iDim (0); iDim < M_fieldDim; ++iDim)
        {
            // Loop over the degrees of freedom (= quadrature nodes)
            for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
            {
                point [0] = interpCFE.quadNode ( iterDof, 0 );
                point [1] = interpCFE.quadNode ( iterDof, 1 );
                point [2] = interpCFE.quadNode ( iterDof, 2 );
                // Store the nodal value
                nodalValues[iterDof] =  fEFunction->eval ( iterElement, point, time );
            }

            // Transform the nodal values in FE values
            FEValues = M_refFE->nodalToFEValues ( nodalValues );

            // Then on the dimension of the FESpace (scalar field vs vectorial field)
            for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
            {
                // Find the ID of the considered DOF
                const ID globalDofID ( M_dof->localToGlobalMap ( iterElement, iterDof ) + iDim * M_dim );

                // Compute the value of the function and set it
                vector.setCoefficient ( globalDofID, FEValues[iterDof] );

            }
        }
    }

} // interpolate

template <typename MeshType, typename MapType>
template<typename vector_type>
void
FESpace<MeshType, MapType>::interpolateBC ( BCHandler& BCh,
                                            vector_type&    vect,
                                            const Real      time)
{
    //
    // ID   nbComp   = M_fieldDim; // Number of components of the mesh velocity
    UInt totalDof = dof().numTotalDof();

    //
    ID idDof;

    if ( !BCh.bcUpdateDone() )
    {
        BCh.bcUpdate ( *mesh(), feBd(), dof() );
    }


    for ( UInt ibc = 0; ibc < BCh.size(); ++ibc )
    {
        if (BCh[ ibc ].type() == Essential)
        {
            // Number of components involved in this boundary condition
            UInt nComp = BCh[ibc].numberOfComponents();

            for ( ID i = 0; i < BCh[ibc].list_size(); ++i )
            {
                // Coordinates of the node where we impose the value
                Real x = static_cast< const BCIdentifierEssential* > ( BCh[ibc][ i ] ) ->x();
                Real y = static_cast< const BCIdentifierEssential* > ( BCh[ibc][ i ] ) ->y();
                Real z = static_cast< const BCIdentifierEssential* > ( BCh[ibc][ i ] ) ->z();

                for ( ID j = 0; j < nComp; ++j )
                {
                    // Global Dof
                    idDof = BCh[ibc][ i ] ->id() + BCh[ibc].component ( j ) * totalDof;
                    Real val = BCh[ibc] ( time, x, y, z, BCh[ibc].component ( j ) );

                    vect.setCoefficient (idDof, val);
                }
            }
        }
    }



}



// this computes vec = \int fct phi_i
template <typename MeshType, typename MapType>
template<typename vector_type>
void
FESpace<MeshType, MapType>::l2ScalarProduct ( const function_Type& fct, vector_type& vec, const Real t)
{

    for ( UInt iVol = 0; iVol < this->mesh()->numElements(); iVol++ )
    {
        this->fe().update ( this->mesh()->element ( iVol ), UPDATE_QUAD_NODES | UPDATE_PHI | UPDATE_WDET );

        Real f, x, y, z;

        UInt i, inod, iQuadPt, ic;
        UInt eleID = this->fe().currentLocalId();
        Real u_ig;

        for ( iQuadPt = 0; iQuadPt < this->fe().nbQuadPt(); iQuadPt++ )
        {
            x = this->fe().quadNode (iQuadPt, 0);
            y = this->fe().quadNode (iQuadPt, 1);
            z = this->fe().quadNode (iQuadPt, 2);
            for ( ic = 0; ic < M_fieldDim; ic++ )
            {
                f = fct ( t, x, y, z, ic );
                u_ig = 0.;
                for ( i = 0; i < this->fe().nbFEDof(); i++ )
                {
                    inod = this->dof().localToGlobalMap ( eleID, i ) + ic * dim();
                    u_ig = f * this->fe().phi ( i, iQuadPt );
                    vec.sumIntoGlobalValues (inod, u_ig * this->fe().weightDet ( iQuadPt ) );
                }

            }
        }
    }

    vec.globalAssemble();
}



template <typename MeshType, typename MapType>
template<typename vector_type>
Real
FESpace<MeshType, MapType>::l20Error ( const function_Type& fexact,
                                       const vector_type& vec,
                                       const Real time,
                                       Real* relError )
{
    Real normU     = 0.;
    Real meanU     = 0.;
    Real mass      = 0.;
    Real sumExact2 = 0.;
    Real sumExact1 = 0.;

    for ( UInt iVol = 0; iVol < this->mesh()->numElements(); iVol++ )
    {
        this->fe().update ( this->mesh()->element ( iVol ), UPDATE_QUAD_NODES | UPDATE_PHI | UPDATE_WDET );

        normU += elementaryDifferenceL2NormSquare ( vec, fexact, this->fe(), this->dof(), time, M_fieldDim );

        meanU += elementaryDifferenceIntegral ( vec, fexact, this->fe(), this->dof(), time, M_fieldDim );
        mass += this->fe().measure();
        if (relError)
        {
            sumExact2 += elementaryFctL2NormSquare ( fexact, this->fe(), time, M_fieldDim );
            sumExact1 += elementaryFctIntegral ( fexact, this->fe(), time, M_fieldDim );
        }
    }


    Real sendbuff[5] = {mass, meanU, normU, sumExact1, sumExact2};
    Real recvbuff[5];

    this->map().comm().SumAll (&sendbuff[0], &recvbuff[0], 5);


    mass       = recvbuff[0];
    meanU      = recvbuff[1];
    normU      = recvbuff[2];
    sumExact1 = recvbuff[3];
    sumExact2 = recvbuff[4];

    Real absError = std::sqrt ( normU - meanU * meanU / mass );

    if (relError)
    {
        Real normExact = std::sqrt ( sumExact2 - sumExact1 * sumExact1 / mass );
        *relError = absError / normExact;
    }
    return absError;
}



template <typename MeshType, typename MapType>
template<typename vector_type>
Real
FESpace<MeshType, MapType>::l2Error ( const function_Type&    fexact,
                                      const vector_type& vec,
                                      const Real       time,
                                      Real*            relError )
{
    Real normU       = 0.;
    Real sumExact    = 0.;

    for ( UInt iVol  = 0; iVol < this->mesh()->numElements(); iVol++ )
    {
        //this->fe().updateFirstDeriv( this->mesh()->element( iVol ) );

        // CurrentFE newFE(this->fe().refFE(),this->fe().geoMap(),quadRuleTetra64pt);
        this->fe().update (this->mesh()->element ( iVol ),  UPDATE_QUAD_NODES | UPDATE_PHI | UPDATE_WDET);

        normU += elementaryDifferenceL2NormSquare ( vec, fexact,
                                                    this->fe(),
                                                    //newFE,
                                                    this->dof(),
                                                    time,
                                                    M_fieldDim );
        if (relError)
        {
            sumExact += elementaryFctL2NormSquare ( fexact,
                                                    this->fe(),
                                                    time,
                                                    M_fieldDim );
        }
    }


    Real sendbuff[2] = {normU, sumExact};
    Real recvbuff[2];

    this->map().comm().SumAll (&sendbuff[0], &recvbuff[0], 2);

    normU    = recvbuff[0];
    sumExact = recvbuff[1];

    if (relError)
    {
        *relError = std::sqrt ( normU / sumExact );
    }

    return std::sqrt ( normU );
}

template <typename MeshType, typename MapType>
template<typename function>
Real
FESpace<MeshType, MapType>::l2NormFunction ( const function& f, const Real time)
{
    //
    //ID nbComp = M_fieldDim; // Number of components of the mesh velocity
    //
    Real sumExact = 0.;
    //
    for ( UInt ielem = 0; ielem < this->mesh()->numElements(); ielem++ )
    {
        this->fe().update ( this->mesh()->element ( ielem ),  UPDATE_QUAD_NODES | UPDATE_PHI | UPDATE_WDET  );

        sumExact += elementaryFctL2NormSquare ( f, this->fe(), time, M_fieldDim );
    }

    Real sendbuff[1] = {sumExact};
    Real recvbuff[1];

    this->map().comm().SumAll (&sendbuff[0], &recvbuff[0], 1);

    sumExact    = recvbuff[0];

    return std::sqrt ( sumExact );

}

template<typename MeshType, typename MapType>
template<typename vector_type>
Real
FESpace<MeshType, MapType>:: l2ErrorWeighted (const function_Type&    exactSolution,
                                              const vector_type& solution,
                                              const function_Type&    weight,
                                              const Real         time)
{
    // Check that the vector is repeated (needed!)
    if (solution.mapType() == Unique)
    {
        return l2ErrorWeighted (exactSolution, VectorEpetra (solution, Repeated), weight, time);
    }

    Real sumOfSquare (0.0);

    // Compute the integral on this processor

    for (UInt iVol (0); iVol < this->mesh()->numElements(); ++iVol)
    {
        this->fe().update (this->mesh()->element (iVol), UPDATE_QUAD_NODES | UPDATE_PHI | UPDATE_WDET);

        for (UInt iQuad (0); iQuad < this->fe().nbQuadPt(); ++iQuad)
        {
            Real solutionInQuadNode (0.0);

            for (UInt iDim (0); iDim < M_fieldDim; ++iDim)
            {
                for (UInt iDof (0); iDof < this->fe().nbFEDof(); ++iDof)
                {
                    UInt dofID (this->dof().localToGlobalMap (iVol, iDof) + iDim * this->dof().numTotalDof() );
                    solutionInQuadNode += this->fe().phi (iDof, iQuad) * solution[dofID];
                }

                Real x (this->fe().quadNode (iQuad, 0) );
                Real y (this->fe().quadNode (iQuad, 1) );
                Real z (this->fe().quadNode (iQuad, 2) );
                Real weightInQuadNode (weight (time, x, y, z, iDim) );
                Real exactInQuadNode (exactSolution (time, x, y, z, iDim) );

                sumOfSquare += weightInQuadNode
                               * (exactInQuadNode - solutionInQuadNode)
                               * (exactInQuadNode - solutionInQuadNode)
                               * this->fe().wDetJacobian (iQuad);
            }
        }
    }

    // Communicate the results

    Real sendbuff (sumOfSquare);
    Real recvbuff;

    this->map().comm().SumAll (&sendbuff, &recvbuff, 1);

    sumOfSquare = recvbuff;

    return std::sqrt (sumOfSquare);
}




template <typename MeshType, typename MapType>
template<typename function, typename vector_type>
Real
FESpace<MeshType, MapType>::h1Error ( const function&    fexact,
                                      const vector_type& vec,
                                      const Real       time,
                                      Real*            relError )
{
    Real normU       = 0.;
    Real sumExact    = 0.;

    for ( UInt iVol  = 0; iVol < this->mesh()->numElements(); iVol++ )
    {
        this->fe().updateFirstDerivQuadPt ( this->mesh()->element ( iVol ) );

        normU += elementaryDifferenceH1NormSquare ( vec, fexact,
                                                    this->fe(),
                                                    //p1,
                                                    this->dof(),
                                                    time,
                                                    M_fieldDim );
        if (relError)
        {
            sumExact += elementaryFctH1NormSquare ( fexact,
                                                    this->fe(),
                                                    //p1,
                                                    time,
                                                    M_fieldDim );
        }
    }


    Real sendbuff[2] = {normU, sumExact};
    Real recvbuff[2];

    this->map().comm().SumAll (&sendbuff[0], &recvbuff[0], 2);


    //    int me = this->map().Comm().MyPID();

    //     if (me == 0)
    //     {
    normU    = recvbuff[0];
    sumExact = recvbuff[1];

    if (relError)
    {
        *relError = std::sqrt ( normU / sumExact );
    }
    //     }

    return std::sqrt ( normU );
}

// Warning! When using this function in parallel, the vector vec HAS TO BE REPEATED!
template <typename MeshType, typename MapType>
template<typename vector_type>
Real
FESpace<MeshType, MapType>::l2Norm ( const vector_type& vec)
{
    //
    ID nbComp = M_fieldDim; // Number of components of the mesh velocity
    //
    Real norm = 0.;
    //
    for ( UInt ielem = 0; ielem < this->mesh()->numElements(); ielem++ )
    {
        //UInt elem = M_FESpace.mesh()->element( ielem ).id();
        this->fe().update ( this->mesh()->element ( ielem ), UPDATE_QUAD_NODES | UPDATE_PHI | UPDATE_WDET );
        //
        norm += elementaryL2NormSquare ( vec, this->fe(), this->dof(), nbComp );
    }

    Real sendbuff[1] = {norm};
    Real recvbuff[1];

    this->map().comm().SumAll (&sendbuff[0], &recvbuff[0], 1);

    norm    = recvbuff[0];

    return std::sqrt ( norm );

}


template <typename MeshType, typename MapType>
template<typename vector_type>
Real
FESpace<MeshType, MapType>::h1Norm (const vector_type& vec)
{
    //
    ID nbComp = M_fieldDim; // Number of components of the mesh velocity
    //
    Real norm = 0.;
    //
    for ( UInt ielem = 0; ielem < this->mesh()->numElements(); ielem++ )
    {
        //UInt elem = M_FESpace.mesh()->element( ielem ).id();
        this->fe().updateFirstDerivQuadPt ( this->mesh()->element ( ielem ) );
        //
        //      for ( UInt j = 0; j < nbComp; ++j)
        norm += elementaryH1NormSquare ( vec, this->fe(), this->dof(), nbComp );
    }

    Real sendbuff[1] = {norm};
    Real recvbuff[1];

    this->map().comm().SumAll (&sendbuff[0], &recvbuff[0], 1);

    norm    = recvbuff[0];

    return std::sqrt ( norm );

}



template<typename MeshType, typename MapType>
template<typename point_type, typename vector_type>
Real
FESpace<MeshType, MapType>::
feInterpolateValue (const ID& elementID, const vector_type& solutionVector, const point_type& pt, const UInt& component ) const
{
    // The vector has to be repeated, so if it is not, we make is repeated and call this function again.
    if (solutionVector.mapType() != Repeated )
    {
        vector_type repeatedSolutionVector (solutionVector, Repeated);
        return feInterpolateValue (elementID, repeatedSolutionVector, pt, component);
    }

    // Make sur everything is up to date
    M_fe->update ( M_mesh->element ( elementID ), UPDATE_PHI);

    // Map the point back to the ref FE
    Real hat_x (0);
    Real hat_y (0);
    Real hat_z (0);
    Real x (0);
    Real y (0);
    Real z (0);

    if (pt.size() >= 1)
    {
        x = pt[0];
    }
    if (pt.size() >= 2)
    {
        y = pt[1];
    }
    if (pt.size() >= 3)
    {
        z = pt[2];
    }

    M_fe->coorBackMap (x, y, z, hat_x, hat_y, hat_z);

    // Store the number of local DoF
    UInt nDof (dof().numLocalDof() );
    UInt totalDof (dof().numTotalDof() );

    // Initialization
    Real value (0);

    // Loop over the DoFs
    for (UInt iter_dof (0); iter_dof < nDof ; ++iter_dof)
    {
        // The global ID of the selected dof
        // This should be changed in case of a VECTORIAL FE

        ID globalDofID (component * totalDof + dof().localToGlobalMap (elementID, iter_dof) ); // iter_dof -> dofID

        // Make the accumulation
        //        std::cout << M_refFE->phi(iter_dof, hat_x, hat_y, hat_z) << " " << iter_dof << " " << hat_x << " " << hat_y << " " << hat_z << std::endl;
        value += solutionVector[globalDofID] * M_refFE->phi (iter_dof, hat_x, hat_y, hat_z);
    }

    return value;
}


template<typename MeshType, typename MapType>
template<typename point_type, typename vector_type>
Real
FESpace<MeshType, MapType>::
feInterpolateValueLocal (const ID& elementID, const vector_type& solutionVector, const point_type& pt ) const
{
    // Make sur everything is up to date
    M_fe->update ( M_mesh->element ( elementID ), UPDATE_PHI);

    // Map the point back to the ref FE
    Real hat_x (0);
    Real hat_y (0);
    Real hat_z (0);

    Real x (0);
    Real y (0);
    Real z (0);

    if (pt.size() >= 1)
    {
        x = pt[0];
    }
    if (pt.size() >= 2)
    {
        y = pt[1];
    }
    if (pt.size() >= 3)
    {
        z = pt[2];
    }

    M_fe->coorBackMap (x, y, z, hat_x, hat_y, hat_z);

    // Store the number of local DoF
    UInt nDof (dof().numLocalDof() );

    // Initialization
    Real value (0);

    // Loop over the DoFs
    for (UInt iter_dof (0); iter_dof < nDof ; ++iter_dof)
    {
        // Make the accumulation
        value += solutionVector[iter_dof] * M_refFE->phi (iter_dof, hat_x, hat_y, hat_z);
    }

    return value;
}


template<typename MeshType, typename MapType>
template<typename point_type, typename vector_type>
Real
FESpace<MeshType, MapType>::
feInterpolateGradient (const ID& elementID, const vector_type& solutionVector, const point_type& pt,
                       const UInt& gradientElement, const UInt& component ) const
{
    // The vector has to be repeated, so if it is not, we make is repeated and call this function again.
    if (solutionVector.getMaptype() != Repeated )
    {
        vector_type repeatedSolutionVector (solutionVector, Repeated);
        return feInterpolateGradient (elementID, repeatedSolutionVector, pt, gradientElement, component);
    }


    // Make sur everything is up to date
    M_fe->updateFirstDeriv ( M_mesh->element ( elementID ) );

    // Map the point back to the ref FE
    Real hat_x (0);
    Real hat_y (0);
    Real hat_z (0);
    Real x (0);
    Real y (0);
    Real z (0);

    if (pt.size() >= 1)
    {
        x = pt[0];
    }
    if (pt.size() >= 2)
    {
        y = pt[1];
    }
    if (pt.size() >= 3)
    {
        z = pt[2];
    }

    M_fe->coorBackMap (x, y, z, hat_x, hat_y, hat_z);

    // Store the number of local DoF
    UInt nDof (dof().numLocalDof() );
    UInt totalDof (dof().numTotalDof() );

    // Initialization
    Real grad (0);

    // Loop over the DoFs
    std::vector<Real> invJac (3, 0);
    invJac[0] = M_fe->pointInverseJacobian (hat_x, hat_y, hat_z, gradientElement, 0);
    invJac[1] = M_fe->pointInverseJacobian (hat_x, hat_y, hat_z, gradientElement, 1);
    invJac[2] = M_fe->pointInverseJacobian (hat_x, hat_y, hat_z, gradientElement, 2);

    for (UInt iter_dof (0); iter_dof < nDof ; ++iter_dof)
    {
        // The global ID of the selected dof
        ID globalDofID (totalDof * component + dof().localToGlobalMap (elementID, iter_dof) );

        for (UInt iter_dim (0); iter_dim < 3; ++iter_dim)
        {
            grad += solutionVector (globalDofID) * M_refFE->dPhi (iter_dof, iter_dim, hat_x, hat_y, hat_z)
                    * invJac[iter_dim];
        }
    }

    return grad;
}


template<typename MeshType, typename MapType>
template<typename point_type, typename vector_type>
Real
FESpace<MeshType, MapType>::
feInterpolateGradientLocal (const ID& elementID, const vector_type& solutionVector, const point_type& pt,
                            const UInt& gradientElement ) const
{
    // Make sur everything is up to date
    M_fe->updateFirstDeriv ( M_mesh->element ( elementID ) );

    // Map the point back to the ref FE
    Real hat_x (0);
    Real hat_y (0);
    Real hat_z (0);
    Real x (0);
    Real y (0);
    Real z (0);

    if (pt.size() >= 1)
    {
        x = pt[0];
    }
    if (pt.size() >= 2)
    {
        y = pt[1];
    }
    if (pt.size() >= 3)
    {
        z = pt[2];
    }

    M_fe->coorBackMap (x, y, z, hat_x, hat_y, hat_z);

    // Store the number of local DoF
    UInt nDof (dof().numLocalDof() );

    // Initialization
    Real grad (0);

    // Loop over the DoFs
    for (UInt iter_dof (0); iter_dof < nDof ; ++iter_dof)
    {

        for (UInt iter_dim (0); iter_dim < 3; ++iter_dim)
        {
            grad += solutionVector[iter_dof] * M_refFE->dPhi (iter_dof, iter_dim, hat_x, hat_y, hat_z)
                    * M_fe->pointInverseJacobian (hat_x, hat_y, hat_z, gradientElement, iter_dim);
        }
    }

    return grad;
}


template<typename MeshType, typename MapType>
template<typename vector_type>
vector_type
FESpace<MeshType, MapType>::
feToFEInterpolate (const FESpace<mesh_Type, map_Type>& OriginalSpace,
                   const vector_type& OriginalVector, const MapEpetraType& outputMapType) const
{
    // This method just check that everything is all right and then call the
    // appropriate method __To__Interpolate if needed.

    // First, check that the interpolation is possible
    ASSERT (fieldDim() == OriginalSpace.fieldDim() ||
            OriginalSpace.refFE().type() == FE_RT0_TETRA_3D ||
            OriginalSpace.refFE().type() == FE_RT0_TRIA_2D, "Incompatible field dimension for interpolation");
    ASSERT (refFE().shape() == OriginalSpace.refFE().shape() , "Incompatible element shape for interpolation");


    boost::shared_ptr<vector_type> InterpolatedVectorPtr;

    // If the spaces are the same, just return the original vector
    if (refFE().type() == OriginalSpace.refFE().type() )
    {
        //return OriginalVector;
        InterpolatedVectorPtr.reset ( new vector_type (OriginalVector) );
    }

    // Distinguish the other cases
    else if ( ( (refFE().type() == FE_P1_3D) && ( (OriginalSpace.refFE().type() == FE_P1bubble_3D) || (OriginalSpace.refFE().type() == FE_P2_3D) ) ) ||
              ( (refFE().type() == FE_P1_2D) && ( (OriginalSpace.refFE().type() == FE_P1bubble_2D) || (OriginalSpace.refFE().type() == FE_P2_3D) ) ) ||
              ( (refFE().type() == FE_P1_2D) && ( OriginalSpace.refFE().type() == FE_P2_2D ) ) ||
              ( (refFE().type() == FE_Q1_3D) && ( OriginalSpace.refFE().type() == FE_Q2_3D ) ) ||
              ( (refFE().type() == FE_Q1_2D) && ( OriginalSpace.refFE().type() == FE_Q2_2D ) ) ||
              ( (refFE().type() == FE_P1bubble_3D) && (OriginalSpace.refFE().type() == FE_P1_3D) ) ||
              ( (refFE().type() == FE_P1bubble_2D) && (OriginalSpace.refFE().type() == FE_P1_2D) ) )
    {
        InterpolatedVectorPtr.reset (new vector_type ( linearInterpolate (OriginalSpace, OriginalVector) ) );
    }
    else
    {
        // The following methods work only if the map is repeated.
        // if the original vector is repeated, make it repeated and recall this method.
        if ( OriginalVector.mapType() == Unique )
        {
            const vector_type OriginalRepeated (OriginalVector, Repeated);
            return feToFEInterpolate (OriginalSpace, OriginalRepeated);
        }

        if ( ( (refFE().type() == FE_P2_3D) && ( (OriginalSpace.refFE().type() == FE_P1bubble_3D) || (OriginalSpace.refFE().type() == FE_P1_3D) ) ) ||
                ( (refFE().type() == FE_P2_2D) && ( (OriginalSpace.refFE().type() == FE_P1bubble_2D) || (OriginalSpace.refFE().type() == FE_P1_2D) ) ) ||
                ( (refFE().type() == FE_P2_2D) && ( OriginalSpace.refFE().type() == FE_P1_2D ) ) )
        {
            InterpolatedVectorPtr.reset (new vector_type ( P2Interpolate (OriginalSpace, OriginalVector) ) );
        }

        else if ( ( OriginalSpace.refFE().type() == FE_RT0_TETRA_3D) || (refFE().type() == FE_RT0_TETRA_3D) ||
                  ( OriginalSpace.refFE().type() == FE_RT0_TRIA_2D) || (refFE().type() == FE_RT0_TRIA_2D) )
        {
            if (refFE().type() == FE_P0_3D || refFE().type() == FE_P0_2D)
            {
                InterpolatedVectorPtr.reset (new vector_type ( RT0ToP0Interpolate (OriginalSpace, OriginalVector) ) );
            }
            else
            {
                ERROR_MSG (" The interpolation with this host space has not been yet implemented. Please, add it!");
            }
        }
        else
        {
            InterpolatedVectorPtr.reset (new vector_type ( interpolateGeneric ( OriginalSpace, OriginalVector) ) );
        }
    }

    if (InterpolatedVectorPtr->mapType() == outputMapType)
    {
        return *InterpolatedVectorPtr;
    }
    else
    {
        // Here we do need to use the combine mode "Insert": the default combine mode
        // is "Add". In fact, when we convert from a Repeated to a unique Map, as we pass several times on the same DoF, the values would be added, what would
        // be wrong. Instead, we keep only one value (which anyway should be always the same for continuous
        // finite element).
        return vector_type (*InterpolatedVectorPtr, outputMapType, Insert);
    }
}


template<typename MeshType, typename MapType>
template<typename vector_type>
vector_type
FESpace<MeshType, MapType>::
gradientRecovery (const vector_type& solution, const UInt& dxi) const
{
    if (solution.mapType() != Repeated)
    {
        return gradientRecovery (vector_type (solution, Repeated), dxi);
    }

    Real refElemArea (0); //area of reference element

    //compute the area of reference element
    for (UInt iq = 0; iq < M_Qr->nbQuadPt(); iq++)
    {
        refElemArea += M_Qr->weight (iq);
    }

    // Define a special quadrature rule for the interpolation
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_refFE->shape() ), M_refFE->shape() );
    Real wQuad (refElemArea / M_refFE->nbDof() );

    for (UInt i (0); i < M_refFE->nbDof(); ++i) //nbRefCoor
    {
        interpQuad.addPoint (QuadraturePoint (M_refFE->xi (i), M_refFE->eta (i), M_refFE->zeta (i), wQuad) );
    }

    // Initialization of the vectors
    vector_type patchArea (solution, Repeated);
    patchArea *= 0.0;

    vector_type gradientSum (solution, Repeated);
    gradientSum *= 0.0;


    UInt totalNumberElements (M_mesh->numElements() );
    UInt numberLocalDof (M_dof->numLocalDof() );

    CurrentFE interpCFE (*M_refFE, getGeometricMap (*M_mesh ), interpQuad);

    // Loop over the cells
    for (UInt iterElement (0); iterElement < totalNumberElements; ++iterElement)
    {
        interpCFE.update (mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );

        for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
        {
            for (UInt iDim (0); iDim < M_fieldDim; ++iDim)
            {
                ID globalDofID (dof().localToGlobalMap (iterElement, iterDof) + iDim * dof().numTotalDof() );

                patchArea[globalDofID] += interpCFE.measure();
                for (UInt iterDofGrad (0); iterDofGrad < numberLocalDof; ++iterDofGrad)
                {
                    ID globalDofIDGrad (dof().localToGlobalMap (iterElement, iterDofGrad) + iDim * dof().numTotalDof() );
                    gradientSum[globalDofID] += interpCFE.measure() * solution[globalDofIDGrad] * interpCFE.dphi (iterDofGrad, dxi, iterDof);
                }
            }
        }
    }

    // Assembly
    return vector_type (gradientSum, Unique, Add) / vector_type (patchArea, Unique, Add);
}

template<typename MeshType, typename MapType>
template<typename vector_type>
vector_type
FESpace<MeshType, MapType>::
laplacianRecovery (const vector_type& solution) const
{
    if (solution.getMaptype() != Repeated)
    {
        return laplacianRecovery (vector_type (solution, Repeated) );
    }

    vector_type laplacian (solution);
    laplacian *= 0.0;

    for (UInt iterDim (0); iterDim < 3; ++iterDim)
    {
        laplacian += gradientRecovery (gradientRecovery (solution, iterDim), iterDim);
    }
    return laplacian;
}



// ===================================================
// Set Methods
// ===================================================

// Set FE space (default standard parameters)
template <typename MeshType, typename MapType>
inline void
FESpace<MeshType, MapType>::setSpace ( const std::string& space, UInt dimension )
{

    if (dimension == 2)
    {
        switch ( M_spaceMap[space] )
        {
            case P1 :
                M_refFE = &feTriaP1;
                M_Qr    = &quadRuleTria3pt;
                M_bdQr  = &quadRuleSeg2pt;
                break;

            case P1_HIGH :
                M_refFE = &feTriaP1;
                M_Qr    = &quadRuleTria6pt;
                M_bdQr  = &quadRuleSeg3pt;
                break;

            case P1Bubble :
                M_refFE = &feTriaP1bubble;
                M_Qr    = &quadRuleTria6pt;
                M_bdQr  = &quadRuleSeg2pt;
                break;

            case P2 :
                M_refFE = &feTriaP2;
                M_Qr    = &quadRuleTria6pt;
                M_bdQr  = &quadRuleSeg3pt;
                break;

            case P2_HIGH :
                M_refFE = &feTriaP2;
                M_Qr    = &quadRuleTria7pt;
                M_bdQr  = &quadRuleSeg3pt;
                break;

            default :
                std::cout << "!!! WARNING: Space " << space << "not implemented in FESpace class !!!" << std::endl;
        }
    }
    else
    {
        switch ( M_spaceMap[space] )
        {
            case P1 :
                M_refFE = &feTetraP1;
                M_Qr    = &quadRuleTetra4pt;
                M_bdQr  = &quadRuleTria3pt;
                break;

            case P1_HIGH :
                M_refFE = &feTetraP1;
                M_Qr    = &quadRuleTetra15pt;
                M_bdQr  = &quadRuleTria4pt;
                break;

            case P1Bubble :
                M_refFE = &feTetraP1bubble;
                M_Qr    = &quadRuleTetra64pt;
                M_bdQr  = &quadRuleTria3pt;
                break;

            case P2 :
                M_refFE = &feTetraP2;
                M_Qr    = &quadRuleTetra15pt;
                M_bdQr  = &quadRuleTria4pt;
                break;

            case P2_HIGH :
                M_refFE = &feTetraP2;
                M_Qr    = &quadRuleTetra64pt;
                M_bdQr  = &quadRuleTria4pt;
                break;

            default :
                std::cout << "!!! WARNING: Space " << space << "not implemented in FESpace class !!!" << std::endl;
        }
    }
}

template<typename MeshType, typename MapType>
void
FESpace<MeshType, MapType>::
setQuadRule (const QuadratureRule& Qr)
{
    M_Qr = &Qr;
    M_fe.reset ( new CurrentFE ( *M_refFE, getGeometricMap ( *M_mesh ), *M_Qr ) );
}


// ===================================================
// Private Methods
// ===================================================
/*
template<typename MeshType, typename MapType>
void
FESpace<MeshType,MapType>::
createMap(const commPtr_Type& commptr)
{
// Build Map
MapType map( *M_refFE, *M_mesh, commptr );
// If more than one field is present the map is
// duplicated by offsetting the DOFs
for ( UInt ii = 0; ii < M_fieldDim; ++ii )
    *M_map += map;
}
*/
template<typename MeshType, typename MapType>
void
FESpace<MeshType, MapType>::
createMap (const commPtr_Type& commptr)
{
    // Against dummies
    ASSERT_PRE (this->M_dof->numTotalDof() > 0, " Cannot create FeSpace with no degrees of freedom");
    // get globalElements list from DOF
    std::vector<Int> myGlobalElements ( this->M_dof->globalElements ( *this->M_mesh ) );
    // Create the map
    MapType map ( -1, myGlobalElements.size(), &myGlobalElements[0], commptr );
    // Store the map. If more than one field is present the map is
    // duplicated by offsetting the DOFs
    for ( UInt ii = 0; ii < M_fieldDim; ++ii )
    {
        *M_map += map;
    }
}


template<typename MeshType, typename MapType>
void
FESpace<MeshType, MapType>::
resetBoundaryFE()
{
    if (M_refFE->hasBoundaryFE() )
    {
        M_feBd.reset (new CurrentBoundaryFE ( M_refFE->boundaryFE(), getGeometricMap ( *M_mesh ).boundaryMap(), *M_bdQr ) );
    }
}


template<typename MeshType, typename MapType>
template<typename vector_type>
vector_type
FESpace<MeshType, MapType>::
interpolateGeneric ( const FESpace<mesh_Type, map_Type>& OriginalSpace,
                     const vector_type& OriginalVector) const
{

    vector_type Interpolated (map(), Repeated);

    // Initialization of the dimensions of the vectors
    UInt size_comp_in = OriginalSpace.dim();

    // number of local DOF
    UInt numberLocalDof (M_dof->numLocalDof() );

    Real xi, yi, zi;

    //vector containing basis function evaluated in the nodes of the current FE element
    std::vector<Real> phi (OriginalSpace.dof().numLocalDof() *numberLocalDof, 0.);

    //vector containing the nodal values of the current FE element
    std::vector<std::vector<Real> > nodalValues (M_fieldDim, std::vector<Real> (numberLocalDof, 0.) );

    //vector containing basis function weights on the current FE element
    std::vector<Real> FEValues (numberLocalDof, 0.);

    //Computing basis function values in the nodes of the current FE element (same values for every element)
    for ( UInt iterDof = 0 ; iterDof < numberLocalDof ; ++iterDof )
    {
        xi = M_refFE->xi (iterDof);
        yi = M_refFE->eta (iterDof);
        zi = M_refFE->zeta (iterDof);
        for (UInt iterDofOrig (0) ; iterDofOrig < OriginalSpace.dof().numLocalDof(); ++iterDofOrig)
        {
            phi[iterDofOrig * numberLocalDof + iterDof] += OriginalSpace.refFE().phi (iterDofOrig, xi, yi, zi);
        }
    }

    //Loop over the mesh elements
    for ( UInt iElem (0); iElem < M_mesh->numElements(); ++iElem )
    {
        //computation the nodal values on each element
        for ( UInt iterDof = 0 ; iterDof < numberLocalDof ; ++iterDof )
            for (UInt iterDofOrig (0) ; iterDofOrig < OriginalSpace.dof().numLocalDof(); ++iterDofOrig)
            {
                size_t index = OriginalSpace.dof().localToGlobalMap ( iElem, iterDofOrig );
                for ( UInt iDim = 0; iDim < M_fieldDim; ++iDim )
                {
                    nodalValues[iDim][iterDof] += OriginalVector[iDim * size_comp_in + index] * phi[iterDofOrig * numberLocalDof + iterDof];
                }
            }

        //loop over field dimension
        for ( UInt iDim = 0; iDim < M_fieldDim; ++iDim )
        {
            FEValues = M_refFE->nodalToFEValues (nodalValues[iDim]); //needed for non-Lagrangian elements
            for ( UInt iterDof (0) ; iterDof < numberLocalDof; ++iterDof )
            {
                ID globalDofId (M_dof->localToGlobalMap ( iElem, iterDof ) + iDim * M_dim);

                // compute the value of the function and set it
                Interpolated.setCoefficient (globalDofId, FEValues[iterDof]);
                nodalValues[iDim][iterDof] = 0;
            }
        }
    }

    return Interpolated;
}

// This method returns a vector with Unique map.
//P2, P1bubble -> P1
//P1 -> P1bubble
//Q2  -> Q1
template<typename MeshType, typename MapType>
template<typename vector_type>
vector_type
FESpace<MeshType, MapType>::
linearInterpolate (const FESpace<mesh_Type, map_Type>& OriginalSpace,
                   const vector_type& OriginalVector) const
{
    // Create a zero vector.
    vector_type Interpolated (map(), Unique);
    UInt totalDofsOriginal (OriginalSpace.dof().numTotalDof() );
    UInt totalDofsPresent (dof().numTotalDof() );
    UInt myLinearDofs;

    if (totalDofsPresent > totalDofsOriginal) //P1 -> P1bubble case
    {
        myLinearDofs = OriginalSpace.map().map (Unique)->NumMyElements() / fieldDim();
    }
    else                               //other cases
    {
        myLinearDofs = map().map (Unique)->NumMyElements() / fieldDim();
    }

    //we exploit the fact that the first DOFs of P1b, P2 vectors are the same of P1 ones.
    for (UInt j = 0; j <  myLinearDofs; j++)
    {
        UInt ig1 = map().map (Unique)->MyGlobalElements() [j];
        UInt ig2 = OriginalSpace.map().map (Unique)->MyGlobalElements() [j];
        for (UInt iComponent (0); iComponent < fieldDim(); ++iComponent)
        {
            Interpolated[ig1 + iComponent * totalDofsPresent] = OriginalVector[ig2 + iComponent * totalDofsOriginal];
        }
    }

    return Interpolated;
}


//P1 -> P2
//P1bubble -> P2
template<typename MeshType, typename MapType>
template<typename vector_type>
vector_type
FESpace<MeshType, MapType>::
P2Interpolate (const FESpace<mesh_Type, map_Type>& OriginalSpace,
               const vector_type& OriginalVector) const
{
    // Create a zero vector.
    vector_type Interpolated (map(), Repeated);
    Interpolated *= 0.0;

    // Some constants (avoid recomputing them each time used)
    UInt FieldDim (fieldDim() );

    UInt numElements (mesh()->numElements() );

    UInt totalDofsOriginal (OriginalSpace.dof().numTotalDof() );
    UInt totalDofsPresent (dof().numTotalDof() );
    UInt localDofsPresent (dof().numLocalDof() );

    // A vector to store the values to set in the dofs
    std::vector<Real> DofValues (localDofsPresent, 0.0);

    // Loop over the elements to get the values
    for ( ID iElem = 0; iElem < numElements ; ++iElem )
    {
        UInt elemId (mesh()->element (iElem).localId() );

        // In the file /lifefem/ReferenceElement.hpp, we can see that the dofs for P1
        // are the first ones of the P2 dofs. We have then to recompute the values
        // in the "face" dofs of the P2 element. Moreover, the bubble is zero on the
        // faces of the element. It gives then no contribution.

        for (UInt iComponent (0); iComponent < FieldDim; ++iComponent)
        {
            // Get the values in the vertices
            for (UInt iP1dof (0); iP1dof < 4; ++iP1dof)
            {
                ID globalDofID_original (iComponent * totalDofsOriginal + OriginalSpace.dof().localToGlobalMap (elemId, iP1dof) );

                DofValues[iP1dof]  =  OriginalVector[globalDofID_original];
            }

            // Compute the values in the faces
            DofValues[4] = 0.5 * (DofValues[0] + DofValues[1]);
            DofValues[5] = 0.5 * (DofValues[1] + DofValues[2]);
            DofValues[6] = 0.5 * (DofValues[0] + DofValues[2]);
            DofValues[7] = 0.5 * (DofValues[0] + DofValues[3]);
            DofValues[8] = 0.5 * (DofValues[1] + DofValues[3]);
            DofValues[9] = 0.5 * (DofValues[2] + DofValues[3]);

            // Now set them
            for (UInt iP2dof (0); iP2dof < localDofsPresent; ++iP2dof)
            {
                ID globalDofID_present (iComponent * totalDofsPresent + dof().localToGlobalMap (elemId, iP2dof) );

                Interpolated[globalDofID_present] = DofValues[iP2dof];
            }

        }
    }

    return Interpolated;
}




template<typename MeshType, typename MapType>
template<typename vector_type>
vector_type
FESpace<MeshType, MapType>::
RT0ToP0Interpolate (const FESpace<mesh_Type, map_Type>& OriginalSpace,
                    const vector_type& OriginalVector) const
{
    // Create a zero vector.
    vector_type Interpolated (map(), Repeated);
    Interpolated *= 0.0;

    // Field dimension of the problem.
    const UInt FieldDim (fieldDim() );

    // Total number of mesh elements.
    const UInt numElements (mesh()->numElements() );

    // Total d.o.f. in the present FESpace.
    const UInt totalDofsPresent (dof().numTotalDof() );

    // Loop over the elements to get the values. To compute the value we use the Piola transformation.
    for ( ID iElem (0); iElem < numElements ; ++iElem )
    {

        // Map between local and global mesh.
        const UInt elemId (mesh()->element (iElem).localId() );

        // Current and reference barycenter, Jacobian of the map.
        Vector3D barCurrentFE, barRefFE, Jac;

        // Update the current element of the P0 vector space.
        M_fe->update (this->mesh()->element (elemId), UPDATE_QUAD_NODES | UPDATE_PHI | UPDATE_WDET);

        // Store the number of local DoF
        const UInt nDof (OriginalSpace.dof().numLocalDof() );

        // Get the coordinate of the barycenter of the current element of ID iElem.
        M_fe->barycenter ( barCurrentFE[0], barCurrentFE[1], barCurrentFE[2] );

        // Update the barycenter to the reference finite element.
        M_fe->coorBackMap ( barCurrentFE[0], barCurrentFE[1], barCurrentFE[2],
                            barRefFE[0], barRefFE[1], barRefFE[2] );

        // Compute the determinant of the Jacobian in the barycenter of the reference finite element.
        const Real det = M_fe->pointDetJacobian (barRefFE[0], barRefFE[1], barRefFE[2]);

        // Loop on the vector component of P0 element.
        for (UInt iComponent (0); iComponent < FieldDim; ++iComponent)
        {
            Real value (0);

            // The Jacobian evaluated in the three points of the barycenter of each component.
            for ( UInt localComponent (0); localComponent < FieldDim; ++localComponent )
            {
                Jac[localComponent] = M_fe->pointJacobian ( barRefFE[0], barRefFE[1], barRefFE[2],
                                                            iComponent, localComponent);
            }

            const UInt iGlobalFacePresent ( iComponent * totalDofsPresent + dof().localToGlobalMap (elemId, 0) );

            // At each face loop on all the d.o.f.
            for (UInt iter_dof (0); iter_dof < nDof ; ++iter_dof)
            {
                // Map between local d.o.f. and global d.o.f.
                const ID globalDofID ( OriginalSpace.dof().localToGlobalMap ( elemId, iter_dof) );

                // Find the correct position in the final vector.
                const UInt iGlobalFacet ( mesh()->localFacetId ( elemId, iter_dof ) );

                // Select if the current facet is coherent or not with the orientation. If yes use +, if not use -.
                if ( mesh()->facet ( iGlobalFacet ).firstAdjacentElementIdentity() != iElem )
                {
                    // Loop on each component of the selected finite element.
                    for (UInt jComponent (0); jComponent < FieldDim; ++jComponent)
                    {
                        value -= OriginalVector[globalDofID] * Jac[jComponent] *
                                 OriginalSpace.refFE().phi (iter_dof, jComponent, barRefFE[0], barRefFE[1], barRefFE[2]);
                    }
                }
                else
                {
                    // Loop on each component of the selected finite element.
                    for (UInt jComponent (0); jComponent < FieldDim; ++jComponent)
                    {
                        value += OriginalVector[globalDofID] * Jac[jComponent] *
                                 OriginalSpace.refFE().phi (iter_dof, jComponent, barRefFE[0], barRefFE[1], barRefFE[2]);
                    }
                }
            }

            // Insert the final value.
            Interpolated[iGlobalFacePresent] = value / det;
        }
    }

    return Interpolated;
}


template<typename MeshType, typename MapType>
UInt
FESpace<MeshType, MapType>::
polynomialDegree() const
{
    switch (M_refFE->type() )
    {
        case FE_P0_0D:
        case FE_P0_2D:
        case FE_Q0_2D:
        case FE_P0_3D:
        case FE_Q0_3D:
            return 0;
            break;

        case FE_P1_1D:
        case FE_P1_2D:
        case FE_Q1_2D:
        case FE_P1_3D:
        case FE_Q1_3D:
            return 1;
            break;

        case FE_P2_1D:
        case FE_P2_2D:
        case FE_Q2_2D:
        case FE_P2_3D:
        case FE_Q2_3D:
            return 2;
            break;

        case FE_P1bubble_2D:
            return 3;
            break;

        case FE_P1bubble_3D:
        case FE_P2tilde_3D:
            return 4;
            break;

        default:
            std::cerr << " FESpace: No polynomial degre for this type of finite element " << std::endl;
            std::cerr << " FESpace: " << M_refFE->name() << std::endl;
            abort();
    };

    return 0;
}

} // end of the namespace
#endif
