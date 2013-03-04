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
    @brief File containing the IP stabilization for the ADR problem.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 06-10-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    A more detailed description of the file (if necessary)
 */

#ifndef ADRASSEMBLERIP_H
#define ADRASSEMBLERIP_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/scoped_ptr.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>

namespace LifeV
{

//! ADRAssemblerIP - This class is used to add IP stabilization to an Advection-Diffusion-Reaction problem.
/*!

  <b> Scope </b>

  This is class has been designed to add the IP stabilization for the ADR problem (see in the ADRAssembler class
  for more detail about that problem). IP stands for Interior Penalty. Indeed, the idea of this stabililization
  is to penalize the solutions that would exhibit oscillations. To quantify the smoothness of the solution, one
  can rely on the computation of the jump of gradient across the faces:

  \f[ \sum_K \sum_{F \subset K} \int_F [\![ \nabla u ]\!]  [\![ \nabla v ]\!]  \f]

  where K represents an element and F one of its faces. The IP stabilization consists then in adding
  the new bilinear form

  \f[ \sum_K \sum_{F \subset K} \int_F k [\![ \nabla u ]\!]  [\![ \nabla v ]\!]  \f]

  in the formulation of the ADR problem to force the solution to be smooth. The coefficient \f$ k \f$
  is used to control the strength of the penalization as well as the good convergence behaviour of the
  stabilized scheme.

  In this class, two choices are possible for \f$ k \f$

  <ul>
  <li> Usually, the choice is \f$ k = \gamma h_F^2 \f$ to ensure the near-optimal convergence of the scheme, \f$ \gamma \f$
  being a fixed parameter.

  <li> For advection dominated flows, it might be easier to use the formula \f$ k = \gamma |\beta \cdot n_F| h_F^2 \f$.
  The choice for a suitable constant \f$ \gamma \f$ does not depend then on the magnitude of the field \f$ \beta \f$ and
  the penalization is somehow better balanced in the domain.
  </ul>

  <b> Stencil problems </b>

  One of the issues of the IP stabilization is the enlargement of the stencil with respect to the standard
  Galerkin approximation. Indeed, the IP stabilization couples that would not be coupled otherwise: if we
  consider a given face F, all the degrees of freedom of the two adjacent elements are coupled together, including
  degrees of freedom that do not belong to a common element.

  From the computational point of view, this can sometimes be not acceptable: the enlarged stencil increases the
  size of the matrix and (often) the size of the preconditioner (and so the time to compute it). For example, if
  the LU factorization of matrix is computed, the factors will have much more non-zero entries because the new
  entries in the matrix are usually "far" from the diagonal.

  To eleviate this problem in some situation, this class provides the possibility to assemble the stabilization
  in two different matrices:
  <ul>
  <li> In one of the matrices (matrixGalerkin), the terms concerning basis functions that belong to the same element are added
  <li> In the other matrice (matrixExtended), the terms concerning basis functions that do not belong to the same element
  are added.
  </ul>

  This allows to use the first matrix in the system and try to explicit the term corresponding to the second matrix
  (putting it in the right hand side). For example, when dealing with parabolic ADR problems,
  one can use the solution of the previous time step (or an extrapolation of the solution)
  to explicit the stabilization by adding in the right hand side the
  product of the extended matrix with the approximated solution.

  <b> Use </b>

  There are two main kind of method: setup methods and stabilization assembly methods.

  The setup method is used to simply build the few internal data structures that are needed for the assembly.
  One should always start by setting up the FESpaces needed.

  Once this is done, the assembly can be performed. There are here 4 methods, distinguished by 2 options.

  The first option is whether the scaling of the stabilization should depend on the
  advection field. If one wants to have a scaling \f$ |\beta \cdot n| \f$ in front of the
  stabilization, then \beta has to appear in the arguments of the function.

  The second option (more advanced) is whether all the terms of the stabilization have to
  be added to the same matrix. If one want to handle the stencil problem, then two
  distinct matrices have to be passed to the methods.

    @author Samuel Quinodoz

    @remark: no choice for the quadrature here, because of currentBDFE.
 */

template< typename mesh_type, typename matrix_type, typename vector_type>
class ADRAssemblerIP
{
public:

    //! @name Public Types
    //@{

    typedef MapEpetra                             map_type;

    typedef FESpace<mesh_type, map_type> fespace_type;
    typedef boost::shared_ptr<fespace_type>              fespace_ptrType;

    typedef boost::shared_ptr<matrix_type>               matrix_ptrType;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    ADRAssemblerIP();

    //! Destructor
    ~ADRAssemblerIP() {};

    //@}



    //! @name Methods
    //@{

    //! Setup method for the FESpaces
    void setup ( const fespace_ptrType& fespace, const fespace_ptrType& betaFEspace);

    //! This method adds the IP stabilization terms into two matrices depending on the stencil.
    /*!
      This method is the main method of this class. It adds the IP stabilization terms into
      two matrices, matrixGalerkin and matrixExtended. In the matrix Galerkin, the terms
      for the terms of the Galerkin stencil are added. In the matrix Extended, the terms
      for the eventually extra stencil are added. Beta is the transport field (given for the
      betaFESpace), coef is a coefficient that is put in front of all the terms.
     */
    void addIPStabilizationStencil (const matrix_ptrType& matrixGalerkin,
                                    const matrix_ptrType& matrixExtended,
                                    const vector_type& beta,
                                    const Real& coef);

    //! This method adds the IP stabilization on the two matrices, depening on the stencil.
    /*!
      However, it does not take into account the advection field to balance the stabilization.
     */
    void addIPStabilizationStencil (const matrix_ptrType& matrixGalerkin,
                                    const matrix_ptrType& matrixExtended,
                                    const Real& coef);

    //! This method builds the IP stabilization on a unique matrix (both stencils)
    /*!
      This method takes also into account the transport field for wieghting the
      stabilization.
     */
    void addIPStabilization (const matrix_ptrType& matrix,
                             const vector_type& beta,
                             const Real& coef)
    {
        addIPStabilizationStencil (matrix, matrix, beta, coef);
    }

    //! Simplest method to add IP stabilization on the matrix.
    /*!
      Terms for both stencils are added on the same matrix and the coefficient
      acts uniformly on the domain (no scaling with the transport field)
     */
    void addIPStabilization (const matrix_ptrType& matrix,
                             const Real& coef)
    {
        addIPStabilizationStencil (matrix, matrix, coef);
    }

    //@}

private:

    typedef MatrixElemental                               localMatrix_type;
    typedef boost::scoped_ptr<localMatrix_type>          localMatrix_ptrType;

    typedef CurrentFE                             currentFE_type;
    typedef boost::scoped_ptr<currentFE_type>            currentFE_ptrType;

    typedef CurrentBoundaryFE                           currentBdFE_type;
    typedef boost::scoped_ptr<currentBdFE_type>          currentBdFE_ptrType;


    // Finite element space for the unknown
    fespace_ptrType M_fespace;

    // Finite element space for the advection
    fespace_ptrType M_betaFESpace;

    currentBdFE_ptrType M_IPFaceCFE;

    currentFE_ptrType M_IPQuad1CFE;
    currentFE_ptrType M_IPQuad2CFE;

    currentFE_ptrType M_IP1CFE;
    currentFE_ptrType M_IP2CFE;

    currentFE_ptrType M_IPBetaCFE;

    // local matrix for the Galerkin stencil entries
    localMatrix_ptrType M_localIPGalerkin_11;
    localMatrix_ptrType M_localIPGalerkin_22;

    // local matrix for the Extended stencil entries
    localMatrix_ptrType M_localIPExtended_12;
    localMatrix_ptrType M_localIPExtended_21;
};

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename mesh_type, typename matrix_type, typename vector_type>
ADRAssemblerIP< mesh_type, matrix_type, vector_type>::
ADRAssemblerIP() :

    M_fespace(),
    M_betaFESpace(),

    M_IPFaceCFE(),

    M_IPQuad1CFE(),
    M_IPQuad2CFE(),

    M_IP1CFE(),
    M_IP2CFE(),

    M_IPBetaCFE(),

    M_localIPGalerkin_11(),
    M_localIPGalerkin_22(),
    M_localIPExtended_12(),
    M_localIPExtended_21()
{}

// ===================================================
// Methods
// ===================================================

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssemblerIP< mesh_type, matrix_type, vector_type>::
setup ( const fespace_ptrType& fespace, const fespace_ptrType& betaFESpace )
{
    M_fespace = fespace;
    M_betaFESpace = betaFESpace;

    M_IPFaceCFE.reset (new CurrentBoundaryFE (M_fespace->feBd().refFE, M_fespace->feBd().geoMap, M_fespace->feBd().qr) );

    // For the two next CurrentFEs, the quadrature plays no role
    M_IPQuad1CFE.reset (new CurrentFE (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_IPQuad2CFE.reset (new CurrentFE (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );

    // For the three next CurrentFEs, the quadrature will be replaced when computing the
    // IP stabilization
    M_IP1CFE.reset (new CurrentFE (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_IP2CFE.reset (new CurrentFE (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_IPBetaCFE.reset (new CurrentFE (M_betaFESpace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );

    // Local matrices
    M_localIPGalerkin_11.reset (new localMatrix_type (M_fespace->fe().nbFEDof()
                                                      , M_fespace->fieldDim()
                                                      , M_fespace->fieldDim() ) );
    M_localIPGalerkin_22.reset (new localMatrix_type (M_fespace->fe().nbFEDof()
                                                      , M_fespace->fieldDim()
                                                      , M_fespace->fieldDim() ) );
    M_localIPExtended_12.reset (new localMatrix_type (M_fespace->fe().nbFEDof()
                                                      , M_fespace->fieldDim()
                                                      , M_fespace->fieldDim() ) );
    M_localIPExtended_21.reset (new localMatrix_type (M_fespace->fe().nbFEDof()
                                                      , M_fespace->fieldDim()
                                                      , M_fespace->fieldDim() ) );
}

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssemblerIP< mesh_type, matrix_type, vector_type>::
addIPStabilizationStencil (const matrix_ptrType& matrixGalerkin,
                           const matrix_ptrType& matrixExtended,
                           const vector_type& beta,
                           const Real& coef)
{
    if (beta.mapType() != Repeated)
    {
        addIPStabilizationStencil (matrixGalerkin, matrixExtended, vector_type (beta, Repeated), coef);
        return;
    }

    ASSERT (M_fespace != 0, "No FE space for building the IP stabilization! ");
    ASSERT (M_betaFESpace != 0, "No FE space (beta) for building the IP stabilization! ");

    // Some constants
    const UInt nbBoundaryFaces (M_fespace->mesh()->numBFaces() );
    const UInt nbFaces (M_fespace->mesh()->numFaces() );
    const UInt nbQuadPt (M_fespace->bdQr().nbQuadPt() );
    const UInt nbLocalDof (M_fespace->fe().nbFEDof() );
    const UInt nbLocalBetaDof (M_betaFESpace->fe().nbFEDof() );
    const UInt nbComponents (M_fespace->fieldDim() );
    const UInt betaTotalDof (M_betaFESpace->dof().numTotalDof() );


    // Temporaries
    Real localValue_11 (0.0);
    Real localValue_22 (0.0);
    Real localValue_12 (0.0);
    Real localValue_21 (0.0);
    std::vector<Real> betaN (nbQuadPt, 0.0);
    Real hFace2 (0.0);

    // Here instead of looping over the elements, we loop on the faces.
    for (UInt iFace (nbBoundaryFaces); iFace < nbFaces; ++iFace)
    {
        // Get the adjacent elements ID
        const UInt adjacentElement1 (M_fespace->mesh()->face (iFace).firstAdjacentElementIdentity() );
        const UInt adjacentElement2 (M_fespace->mesh()->face (iFace).secondAdjacentElementIdentity() );

        // Here we check that the face is included in the IP
        // stabilization: it is not a boundary face (we do not
        // stabilize there). We also cannot (not possible) stabilize
        // across the different partitions of the mesh (if they exist).
        // These cases are the excluded.

        if ( Flag::testOneSet ( M_fespace->mesh()->face (iFace).flag(), EntityFlags::SUBDOMAIN_INTERFACE | EntityFlags::PHYSICAL_BOUNDARY ) )
        {
            continue;
        };

        // Now the point is that we need to integrate on the face informations coming
        // from the element (beta is defined on the volume, dphi as well). Therefore,
        // we need to update the currentFEs with a quadrature that lies on the face.
        // First step , we compute this quadrature.

        M_IPFaceCFE->updateMeasNormalQuadPt (M_fespace->mesh()->face (iFace) );
        hFace2 = M_IPFaceCFE->measure();

        // Second step, we take the quadrature back to the reference frame for both
        // adjacent elements

        M_IPQuad1CFE->update ( M_fespace->mesh()->element (adjacentElement1), UPDATE_ONLY_CELL_NODES );

        QuadratureRule faceQR1 ("custom quad 1", TETRA, 3, 0, 0);
        for (UInt iQuad (0); iQuad < nbQuadPt; ++iQuad) // Here we do not use UInt because of KNM, but we should
        {
            Real x (0.0), y (0.0), z (0.0);
            M_IPQuad1CFE->coorBackMap ( M_IPFaceCFE->quadPt (iQuad, 0),
                                        M_IPFaceCFE->quadPt (iQuad, 1),
                                        M_IPFaceCFE->quadPt (iQuad, 2),
                                        x, y, z);
            QuadraturePoint newPoint (x, y, z, M_fespace->bdQr().weight (iQuad) );
            faceQR1.addPoint (newPoint);
        }

        M_IPQuad2CFE->update ( M_fespace->mesh()->element (adjacentElement2), UPDATE_ONLY_CELL_NODES );
        QuadratureRule faceQR2 ("custom quad 2", TETRA, 3, 0, 0);
        for (UInt iQuad (0); iQuad < nbQuadPt; ++iQuad) // Idem here
        {
            Real x (0.0), y (0.0), z (0.0);
            M_IPQuad2CFE->coorBackMap ( M_IPFaceCFE->quadPt (iQuad, 0),
                                        M_IPFaceCFE->quadPt (iQuad, 1),
                                        M_IPFaceCFE->quadPt (iQuad, 2),
                                        x, y, z);
            QuadraturePoint newPoint (x, y, z, M_fespace->bdQr().weight (iQuad) );
            faceQR2.addPoint (newPoint);
        }

        // Third step, we change the quadrature in the CurrentFEs

        M_IP1CFE->setQuadRule (faceQR1);
        M_IP2CFE->setQuadRule (faceQR2);
        M_IPBetaCFE->setQuadRule (faceQR1);

        // Now the CurrentFEs are updated with a quadrature that is
        // actually only on the considered face (iterFace).

        M_IP1CFE->update (M_fespace->mesh()->element (adjacentElement1), UPDATE_DPHI | UPDATE_WDET);
        M_IP2CFE->update (M_fespace->mesh()->element (adjacentElement2), UPDATE_DPHI | UPDATE_WDET);
        M_IPBetaCFE->update (M_fespace->mesh()->element (adjacentElement1), UPDATE_PHI );

        // Before starting the assembly, we compute the values of |beta n|
        // in the quadrature nodes

        for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
        {
            betaN[iQuadPt] = 0.0;
            for (UInt iDof (0); iDof < nbLocalBetaDof; ++iDof)
            {
                for (UInt iDim (0); iDim < 3; ++iDim)
                {
                    betaN[iQuadPt] += beta[M_betaFESpace->dof().localToGlobalMap (adjacentElement1, iDof)
                                           + betaTotalDof * iDim]
                                      * M_IPBetaCFE->phi (iDof, iQuadPt)
                                      * M_IPFaceCFE->normal (iDim, iQuadPt);
                }
            }
            betaN[iQuadPt] = std::fabs (betaN[iQuadPt]);
        }

        // Now we can start the assembly. There are 4 parts, depending on
        // which sides we consider ( side1 with side1, side1 with side2,...)
        // We have also to be carefull about the signs, because we want the jumps
        // Therefore, there is a "-" when considering two different sides.

        // Zero out the local matrices
        M_localIPGalerkin_11->zero();
        M_localIPGalerkin_22->zero();
        M_localIPExtended_12->zero();
        M_localIPExtended_21->zero();

        // Loop on the components
        for (UInt iFieldDim (0); iFieldDim < nbComponents; ++iFieldDim)
        {
            // Extract the views
            localMatrix_type::matrix_view viewIPGalerkin_11 = M_localIPGalerkin_11->block (iFieldDim, iFieldDim);
            localMatrix_type::matrix_view viewIPGalerkin_22 = M_localIPGalerkin_22->block (iFieldDim, iFieldDim);
            localMatrix_type::matrix_view viewIPExtended_12 = M_localIPExtended_12->block (iFieldDim, iFieldDim);
            localMatrix_type::matrix_view viewIPExtended_21 = M_localIPExtended_21->block (iFieldDim, iFieldDim);

            for (UInt iDof (0); iDof < nbLocalDof; ++iDof)
            {
                for (UInt jDof (0); jDof < nbLocalDof; ++jDof)
                {
                    localValue_11 = 0.0;
                    localValue_22 = 0.0;
                    localValue_12 = 0.0;
                    localValue_21 = 0.0;

                    for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                    {
                        for (UInt iDim (0); iDim < 3; ++iDim)
                        {
                            localValue_11 += betaN[iQuadPt]
                                             * M_IP1CFE->dphi (iDof, iDim, iQuadPt)
                                             * M_IP1CFE->dphi (jDof, iDim, iQuadPt)
                                             * M_IP1CFE->wDetJacobian (iQuadPt);

                            localValue_22 += betaN[iQuadPt]
                                             * M_IP2CFE->dphi (iDof, iDim, iQuadPt)
                                             * M_IP2CFE->dphi (jDof, iDim, iQuadPt)
                                             * M_IP2CFE->wDetJacobian (iQuadPt);

                            localValue_12 += betaN[iQuadPt]
                                             * M_IP1CFE->dphi (iDof, iDim, iQuadPt)
                                             * M_IP2CFE->dphi (jDof, iDim, iQuadPt)
                                             * M_IP1CFE->wDetJacobian (iQuadPt);

                            localValue_21 += betaN[iQuadPt]
                                             * M_IP2CFE->dphi (iDof, iDim, iQuadPt)
                                             * M_IP1CFE->dphi (jDof, iDim, iQuadPt)
                                             * M_IP1CFE->wDetJacobian (iQuadPt);
                        }
                    }

                    // Here we put the values in the local matrices
                    // We care for sign (to get jumps) and for the
                    // coefficient here.
                    viewIPGalerkin_11 (iDof, jDof) += coef * hFace2 * localValue_11;
                    viewIPGalerkin_22 (iDof, jDof) += coef * hFace2 * localValue_22;
                    viewIPExtended_12 (iDof, jDof) += -coef * hFace2 * localValue_12;
                    viewIPExtended_21 (iDof, jDof) += -coef * hFace2 * localValue_21;
                }
            }
        }

        // Global Assembly
        // We separate here the assembly of the two
        // contributions for Galerkin and Extended
        // stencil.

        for (UInt iFieldDim (0); iFieldDim < nbComponents; ++iFieldDim)
        {
            assembleMatrix ( *matrixGalerkin,
                             *M_localIPGalerkin_11,
                             *M_IP1CFE,
                             *M_IP1CFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * M_fespace->dof().numTotalDof(), iFieldDim * M_fespace->dof().numTotalDof() );

            assembleMatrix ( *matrixGalerkin,
                             *M_localIPGalerkin_22,
                             *M_IP2CFE,
                             *M_IP2CFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * M_fespace->dof().numTotalDof(), iFieldDim * M_fespace->dof().numTotalDof() );

            assembleMatrix ( *matrixExtended,
                             *M_localIPExtended_12,
                             *M_IP1CFE,
                             *M_IP2CFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * M_fespace->dof().numTotalDof(), iFieldDim * M_fespace->dof().numTotalDof() );

            assembleMatrix ( *matrixExtended,
                             *M_localIPExtended_21,
                             *M_IP2CFE,
                             *M_IP1CFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * M_fespace->dof().numTotalDof(), iFieldDim * M_fespace->dof().numTotalDof() );
        }
    }

}

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssemblerIP< mesh_type, matrix_type, vector_type>::
addIPStabilizationStencil (const matrix_ptrType& matrixGalerkin,
                           const matrix_ptrType& matrixExtended,
                           const Real& coef)
{
    ASSERT (M_fespace != 0, "No FE space for building the IP stabilization! ");

    // Some constants
    const UInt nbBoundaryFaces (M_fespace->mesh()->numBFaces() );
    const UInt nbFaces (M_fespace->mesh()->numFaces() );
    const UInt nbQuadPt (M_fespace->bdQr().nbQuadPt() );
    const UInt nbLocalDof (M_fespace->fe().nbFEDof() );
    const UInt nbComponents (M_fespace->fieldDim() );


    // Temporaries
    Real localValue_11 (0.0);
    Real localValue_22 (0.0);
    Real localValue_12 (0.0);
    Real localValue_21 (0.0);
    Real hFace2 (0.0);

    // Here instead of looping over the elements, we loop on the faces.
    for (UInt iFace (nbBoundaryFaces); iFace < nbFaces; ++iFace)
    {
        // Get the adjacent elements ID
        const UInt adjacentElement1 (M_fespace->mesh()->face (iFace).firstAdjacentElementIdentity() );
        const UInt adjacentElement2 (M_fespace->mesh()->face (iFace).secondAdjacentElementIdentity() );

        // Here we check that the face is included in the IP
        // stabilization: it is not a boundary face (we do not
        // stabilize there). We also cannot (not possible) stabilize
        // across the different partitions of the mesh (if they exist).
        // These cases are the excluded.

        if ( (adjacentElement1 == NotAnId) || (adjacentElement2 == NotAnId) || (adjacentElement1 == adjacentElement2) )
        {
            continue;
        };

        // Now the point is that we need to integrate on the face informations coming
        // from the element (beta is defined on the volume, dphi as well). Therefore,
        // we need to update the currentFEs with a quadrature that lies on the face.
        // First step , we compute this quadrature.

        M_IPFaceCFE->updateMeasNormalQuadPt (M_fespace->mesh()->face (iFace) );
        hFace2 = M_IPFaceCFE->measure();

        // Second step, we take the quadrature back to the reference frame for both
        // adjacent elements

        M_IPQuad1CFE->update ( M_fespace->mesh()->element (adjacentElement1), UPDATE_ONLY_CELL_NODES );

        QuadratureRule faceQR1 ("custom quad 1", TETRA, 3, 0, 0);
        for (int iQuad (0); iQuad < nbQuadPt; ++iQuad) // Here we do not use UInt because of KNM, but we should
        {
            Real x (0.0), y (0.0), z (0.0);
            M_IPQuad1CFE->coorBackMap ( M_IPFaceCFE->quadPt (iQuad, 0),
                                        M_IPFaceCFE->quadPt (iQuad, 1),
                                        M_IPFaceCFE->quadPt (iQuad, 2),
                                        x, y, z);
            QuadraturePoint newPoint (x, y, z, M_fespace->bdQr().weight (iQuad) );
            faceQR1.addPoint (newPoint);
        }

        M_IPQuad2CFE->update ( M_fespace->mesh()->element (adjacentElement2), UPDATE_ONLY_CELL_NODES );
        QuadratureRule faceQR2 ("custom quad 2", TETRA, 3, 0, 0);
        for (int iQuad (0); iQuad < nbQuadPt; ++iQuad) // Idem here
        {
            Real x (0.0), y (0.0), z (0.0);
            M_IPQuad2CFE->coorBackMap ( M_IPFaceCFE->quadPt (iQuad, 0),
                                        M_IPFaceCFE->quadPt (iQuad, 1),
                                        M_IPFaceCFE->quadPt (iQuad, 2),
                                        x, y, z);
            QuadraturePoint newPoint (x, y, z, M_fespace->bdQr().weight (iQuad) );
            faceQR2.addPoint (newPoint);
        }

        // Third step, we change the quadrature in the CurrentFEs

        M_IP1CFE->setQuadRule (faceQR1);
        M_IP2CFE->setQuadRule (faceQR2);

        // Now the CurrentFEs are updated with a quadrature that is
        // actually only on the considered face (iterFace).

        M_IP1CFE->update (M_fespace->mesh()->element (adjacentElement1), UPDATE_DPHI | UPDATE_WDET);
        M_IP2CFE->update (M_fespace->mesh()->element (adjacentElement2), UPDATE_DPHI | UPDATE_WDET);

        // Now we can start the assembly. There are 4 parts, depending on
        // which sides we consider ( side1 with side1, side1 with side2,...)
        // We have also to be carefull about the signs, because we want the jumps
        // Therefore, there is a "-" when considering two different sides.

        // Zero out the local matrices
        M_localIPGalerkin_11->zero();
        M_localIPGalerkin_22->zero();
        M_localIPExtended_12->zero();
        M_localIPExtended_21->zero();

        // Loop on the components
        for (UInt iFieldDim (0); iFieldDim < nbComponents; ++iFieldDim)
        {
            // Extract the views
            localMatrix_type::matrix_view viewIPGalerkin_11 = M_localIPGalerkin_11->block (iFieldDim, iFieldDim);
            localMatrix_type::matrix_view viewIPGalerkin_22 = M_localIPGalerkin_22->block (iFieldDim, iFieldDim);
            localMatrix_type::matrix_view viewIPExtended_12 = M_localIPExtended_12->block (iFieldDim, iFieldDim);
            localMatrix_type::matrix_view viewIPExtended_21 = M_localIPExtended_21->block (iFieldDim, iFieldDim);

            for (UInt iDof (0); iDof < nbLocalDof; ++iDof)
            {
                for (UInt jDof (0); jDof < nbLocalDof; ++jDof)
                {
                    localValue_11 = 0.0;
                    localValue_22 = 0.0;
                    localValue_12 = 0.0;
                    localValue_21 = 0.0;

                    for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                    {
                        for (UInt iDim (0); iDim < 3; ++iDim)
                        {
                            localValue_11 += 1.0
                                             * M_IP1CFE->dphi (iDof, iDim, iQuadPt)
                                             * M_IP1CFE->dphi (jDof, iDim, iQuadPt)
                                             * M_IP1CFE->wDetJacobian (iQuadPt);

                            localValue_22 += 1.0
                                             * M_IP2CFE->dphi (iDof, iDim, iQuadPt)
                                             * M_IP2CFE->dphi (jDof, iDim, iQuadPt)
                                             * M_IP2CFE->wDetJacobian (iQuadPt);

                            localValue_12 += 1.0
                                             * M_IP1CFE->dphi (iDof, iDim, iQuadPt)
                                             * M_IP2CFE->dphi (jDof, iDim, iQuadPt)
                                             * M_IP1CFE->wDetJacobian (iQuadPt);

                            localValue_21 += 1.0
                                             * M_IP2CFE->dphi (iDof, iDim, iQuadPt)
                                             * M_IP1CFE->dphi (jDof, iDim, iQuadPt)
                                             * M_IP1CFE->wDetJacobian (iQuadPt);
                        }
                    }

                    // Here we put the values in the local matrices
                    // We care for sign (to get jumps) and for the
                    // coefficient here.
                    viewIPGalerkin_11 (iDof, jDof) += coef * hFace2 * localValue_11;
                    viewIPGalerkin_22 (iDof, jDof) += coef * hFace2 * localValue_22;
                    viewIPExtended_12 (iDof, jDof) += -coef * hFace2 * localValue_12;
                    viewIPExtended_21 (iDof, jDof) += -coef * hFace2 * localValue_21;
                }
            }
        }

        // Global Assembly
        // We separate here the assembly of the two
        // contributions for Galerkin and Extended
        // stencil.

        for (UInt iFieldDim (0); iFieldDim < nbComponents; ++iFieldDim)
        {
            assembleMatrix ( *matrixGalerkin,
                             *M_localIPGalerkin_11,
                             *M_IP1CFE,
                             *M_IP1CFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * M_fespace->dof().numTotalDof(), iFieldDim * M_fespace->dof().numTotalDof() );

            assembleMatrix ( *matrixGalerkin,
                             *M_localIPGalerkin_22,
                             *M_IP2CFE,
                             *M_IP2CFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * M_fespace->dof().numTotalDof(), iFieldDim * M_fespace->dof().numTotalDof() );

            assembleMatrix ( *matrixExtended,
                             *M_localIPExtended_12,
                             *M_IP1CFE,
                             *M_IP2CFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * M_fespace->dof().numTotalDof(), iFieldDim * M_fespace->dof().numTotalDof() );

            assembleMatrix ( *matrixExtended,
                             *M_localIPExtended_21,
                             *M_IP2CFE,
                             *M_IP1CFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * M_fespace->dof().numTotalDof(), iFieldDim * M_fespace->dof().numTotalDof() );
        }
    }

}

} // Namespace LifeV

#endif /* ADRASSEMBLERIP_H */
