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
    @brief File containing the ADRAssembler class.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 28-09-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */

#ifndef ADRASSEMBLER_H
#define ADRASSEMBLER_H 1


#include <boost/scoped_ptr.hpp>


#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>

namespace LifeV
{

//! ADRAssembler - This class add into given matrices terms corresponding to the space discretization of the ADR problem.
/*!

  <b> Scope </b>

  This class has been designed for assembling the terms arising
  from the space discretization with finite elements of the advection-diffusion-reaction
  problem, i.e. of the PDE

  \f$ - \alpha \Delta u + \beta \cdot \nabla u + \sigma u = f \f$

  where \f$\alpha\f$ and \f$\sigma\f$ are constants, \f$\beta\f$ a vectorial field.

  Time discretization of the parabolic version of this equation is not in the scope of
  this assembler (it is usually not a finite element discretization) and has to be
  treated outside of the class.

  This class is also not supposed to provide any stabilization for the discretization.


  <b> Good Practice </b>

  Here is the way this class is intended to be used. First, one should define the
  assembler. Only the empty constructor is available, so there is no choice.

  \code
  ADRAssembler<mesh_type,matrix_type,vector_type> myAssembler;
  \endcode

  Then, one should setup the assembler by providing the finite element spaces
  required, namely the space for the unknown and the space where the advection
  field \f$\beta\f$ is defined

  \code
  boost::shared_ptr<FESpace< mesh_type, MapEpetra > > uFESpace( new FESpace< mesh_type, MapEpetra >( ... ));
  boost::shared_ptr<FESpace< mesh_type, MapEpetra > > betaFESpace( new FESpace< mesh_type, MapEpetra >( ... ));

  myAssembler.setup(uFESpace,betaFESpace);
  \endcode

  When this setup method is used, the quadrature rules are set as the one
  stored in the FESpace of the unknown. One should then change them if one wants to
  have a more precise or a faster assembly (usually, the default quadrature is too
  precise for diffusion and the advection terms, so one should consider this option
  to make the assembly faster without affecting the results).

  \code
  myAssembler.setQuadRuleForDiffusionMatrix( ... );
  myAssembler.setQuadRuleForAdvectionMatrix( ... );
  \endcode

  Now that everything has been set up, one can assemble the terms:

  \code
  boost::shared_ptr<matrix_type> systemMatrix(new matrix_type( uFESpace->map() ));
  *systemMatrix*=0.0;

  myAssembler.addDiffusion(systemMatrix,4.0);
  myAssembler.addMass(systemMatrix);
  \endcode

  These two last lines have added the required terms into the matrix. Here one should
  be aware that the matrix is NOT finalized after the assembly of any of the terms.
  Therefore, before passing the matrix to a linear solver (and before applying boundary
  conditions), one should close the matrix using

  \code
  systemMatrix->GlobalAssemble();
  \endcode


  <b> Remarks </b>

  <ol>
  <li> Changing the quadratures of the assembler is in many case a good idea. However, try to
  avoid calling the setters for the quadratures when it is not necessary, as these methods change
  some internal containers and are then more expensive than "simple setters".
  <li> In case the FESpace for the unknown is vectorial (and not scalar), the problem is assembled
  for each component of the FESpace.
  <li> The matrices assembled are NOT stored in this class, so calling twice the same assembly method
  will not make it faster. If your code calles repeatedly the assembly procedures, consider storing the
  matrices outside the class.
  </ol>


  <b> Troubleshooting </b>

  <i> Empty FE Space cannot setup the ADR assembler </i> You are passing a nul pointer as FE Space in the setup method.

  <i> Empty beta FE Space cannot setup the ADR assembler </i> You are passing a nul pointer as FE Space for the advection field in the setup method.

  <i> Setting the FE space for the unknown to 0 is not permitted </i> You are passing a nul pointer as FE Space in the setFespace method.

  <i> No FE space for the unknown! Use setFespace before setBetaFespace! </i> You are trying to use the method setBetaFespace before having set the FE space for the unknown.

  <i> No FE space for assembling the mass! </i> You are trying to use the addMass method without having set a FE Space for the unknown. Use the setup or the setFespace methods before calling addMass.

  <i> No FE space for assembling the advection! </i> You are trying to use the addAdvection method without having set a FE Space for the unknown. Use the setup or the setFespace methods before calling addAdvection.

  <i> No FE space for assembling the diffusion! </i> You are trying to use the addDiffusion method without having set a FE Space for the unknown. Use the setup or the setFespace methods before calling addDiffusion.

  <i> No FE space (beta) for assembling the advection! </i> You are trying to add advection without having set the FE space for the advection field. Use the setup of the setBetaFespace methods before calling addAdvection.


  @author Samuel Quinodoz
  @version 1.1
*/

template< typename mesh_type, typename matrix_type, typename vector_type>
class ADRAssembler
{
public:

    //! @name Public Types
    //@{

    typedef MapEpetra                                    map_type;

    typedef FESpace<mesh_type, map_type> fespace_type;
    typedef boost::shared_ptr<fespace_type>              fespace_ptrType;

    typedef boost::shared_ptr<matrix_type>               matrix_ptrType;

    typedef LifeChrono                                       chrono_type;

    // Use the portable syntax of the boost function
    typedef boost::function5<Real, const Real&, const Real&, const Real&, const Real&, const ID&> function_type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    ADRAssembler();

    //! Destructor
    virtual ~ADRAssembler() {}

    //@}



    //! @name Methods
    //@{

    //! Setup for the class (fill the members)
    /*!
      This method can be called either when the class is empty or when is already (partially)
      setup. It changes the internal containers to fit with the FESpace given in argument.
      All the quadrature rules are overriden by the one given in fespace.
      @param fespace The FESpace for the unknown
      @param betaFESpace The FESpace for the advection field
      @remark Nul pointers cannot be passed to this method. Use a more
      specific setter if you want to do that.
     */
    void setup (const fespace_ptrType& fespace, const fespace_ptrType& betaFESpace);

    //! Assembling for the mass
    /*!
      This method adds the mass matrix times the coefficient
      (whose default value is 1.0) to the matrix passed in argument.
      @Remark The matrix is NOT finalized, you have to call globalAssemble
      outside this method when the matrix is finished.
     */
    inline void addMass (matrix_ptrType matrix, const Real& coefficient = 1.0)
    {
        addMass (matrix, coefficient, 0, 0);
    }

    //! Assembling for the mass using offsets
    /*!
      This method adds the mass in the given matrix. Additional arguements are provided
      so that one can chose where to add the mass in the matrix, i.e. somewhere else
      than in the upper left block.
     */
    void addMass (matrix_ptrType matrix, const Real& coefficient, const UInt& offsetLeft, const UInt& offsetUp);

    //! Assembling for the advection
    /*!
      This method adds the advection (with respect to the
      given vector field) matrix to the matrix passed in argument.
      beta represents a finite element function with the beta FE space
      of this assembler.
      @Remark The matrix is NOT finalized, you have to call globalAssemble
      outside this method when the matrix is finished.
     */
    inline void addAdvection (matrix_ptrType matrix, const vector_type& beta)
    {
        addAdvection (matrix, beta, 0, 0);
    }

    //! Assembling for the advection using offsets
    /*!
      This method adds the advection in the given matrix. Additional arguements are provided
      so that one can chose where to add the mass in the matrix, i.e. somewhere else
      than in the upper left block.
     */
    void addAdvection (matrix_ptrType matrix, const vector_type& beta, const UInt& offsetLeft, const UInt& offsetUp);

    //! Assembling for the diffusion
    /*!
      This method adds the diffusion matrix times the coefficient
      (whose default value is 1.0) to the matrix passed in argument.
      @Remark The matrix is NOT finalized, you have to call globalAssemble
      outside this method when the matrix is finished.
     */
    inline void addDiffusion (matrix_ptrType matrix, const Real& coefficient = 1.0)
    {
        addDiffusion (matrix, coefficient, 0, 0);
    }

    //! Assembling for the diffusion using offsets
    /*!
      This method adds the diffusion in the given matrix. Additional arguements are provided
      so that one can chose where to add the mass in the matrix, i.e. somewhere else
      than in the upper left block.
     */
    void addDiffusion (matrix_ptrType matrix, const Real& coefficient, const UInt& offsetLeft, const UInt& offsetUp);

    //! Assembly for the right hand side (mass) with f given in vectorial form.
    /*!
      This method assembles the right hand side for the ADR problem
      where the forcing term is given in the FE space of the unknown.
     */
    void addMassRhs (vector_type& rhs, const vector_type& f);

    //! Assembly for the right hand side (mass) with f given in functional form.
    /*!
      This method assembles the right hand side for the ADR problem
      where f is given as a function of space and time (t,x,y,z,component)
     */
    void addMassRhs (vector_type& rhs, const function_type& f, const Real& t);

    //@}



    //! @name Set Methods
    //@{

    //! Setter for the finite element space used for the unknown (reset the quadratures as well!)
    void setFespace (const fespace_ptrType& fespace);

    //! Setter for the finite element space used for the advection field (reset the quadratures as well!)
    /*!
      Beware that a FE space for the unknown has to be set before calling this method.
     */
    void setBetaFespace (const fespace_ptrType& betaFESpace);

    //! Setter for the quadrature used for the mass matrix
    /*!
      Beware that calling this function might be quite heavy, so avoid using
      it when it is not necessary.
     */
    inline void setQuadRuleForMassMatrix (const QuadratureRule& qr)
    {
        ASSERT (M_massCFE != 0, "No mass currentFE for setting the quadrature rule!");
        M_massCFE->setQuadRule (qr);
    }

    //! Setter for the quadrature used for the diffusion matrix
    /*!
      Beware that calling this function might be quite heavy, so avoid using
      it when it is not necessary.
     */
    inline void setQuadRuleForDiffusionMatrix (const QuadratureRule& qr)
    {
        ASSERT (M_diffCFE != 0, "No diffusion currentFE for setting the quadrature rule!");
        M_diffCFE->setQuadRule (qr);
    }

    //! Setter for the quadrature used for the advection matrix
    /*!
      Beware that calling this function might be quite heavy, so avoid using
      it when it is not necessary.
     */
    inline void setQuadRuleForAdvectionMatrix (const QuadratureRule& qr)
    {
        ASSERT (M_advCFE != 0, "No advection (u) currentFE for setting the quadrature rule!");
        ASSERT (M_advBetaCFE != 0, "No advection (beta) currentFE for setting the quadrature rule!");
        M_advCFE->setQuadRule (qr);
        M_advBetaCFE->setQuadRule (qr);
    }

    //! Setter for the quadrature used for the right hand side
    /*!
      Beware that calling this function might be quite heavy, so avoid using
      it when it is not necessary.
    */
    inline void setQuadRuleForMassRhs (const QuadratureRule& qr)
    {
        ASSERT (M_massRhsCFE != 0, "No Rhs currentFE for setting the quadrature rule!");
        M_massRhsCFE->setQuadRule (qr);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the chrono of the assembly of the mass term
    chrono_type& massAssemblyChrono()
    {
        return M_massAssemblyChrono;
    }

    //! Getter for the chrono of the assembly of the advection term
    chrono_type& advectionAssemblyChrono()
    {
        return M_advectionAssemblyChrono;
    }

    //! Getter for the chrono of the assembly of the diffusion term
    chrono_type& diffusionAssemblyChrono()
    {
        return M_diffusionAssemblyChrono;
    }

    //! Getter for the chrono of the setup of the assembler
    chrono_type& setupChrono()
    {
        return M_setupChrono;
    }

    //! Getter for the chrono of the assembly of the rhs (mass)
    chrono_type& massRhsAssemblyChrono()
    {
        return M_massRhsAssemblyChrono;
    }

    //@}

private:


    typedef CurrentFE                                    currentFE_type;
    typedef boost::scoped_ptr<currentFE_type>            currentFE_ptrType;

    typedef MatrixElemental                                      localMatrix_type;
    typedef boost::scoped_ptr<localMatrix_type>          localMatrix_ptrType;

    typedef VectorElemental                                      localVector_type;
    typedef boost::scoped_ptr<localVector_type>          localVector_ptrType;

    //! @name Private Methods
    //@{

    // Copy constructor is a no-sense.
    ADRAssembler (const ADRAssembler&);

    //@}

    // Finite element space for the unknown
    fespace_ptrType M_fespace;

    // Finite element space for the advection
    fespace_ptrType M_betaFESpace;

    // CurrentFE for the mass
    currentFE_ptrType M_massCFE;

    // CurrentFE for the diffusion
    currentFE_ptrType M_diffCFE;

    // CurrentFE for the advection (unknown)
    currentFE_ptrType M_advCFE;

    // CurrentFE for the advection (beta)
    currentFE_ptrType M_advBetaCFE;

    // CurrentFE for the mass rhs
    currentFE_ptrType M_massRhsCFE;


    // Local matrix for the mass
    localMatrix_ptrType M_localMass;

    // Local advection matrix
    localMatrix_ptrType M_localAdv;

    // Local matrix for the diffusion
    localMatrix_ptrType M_localDiff;

    // Local vector for the right hand side
    localVector_ptrType M_localMassRhs;

    // Chronos
    chrono_type M_diffusionAssemblyChrono;
    chrono_type M_advectionAssemblyChrono;
    chrono_type M_massAssemblyChrono;
    chrono_type M_setupChrono;
    chrono_type M_massRhsAssemblyChrono;

};



template<typename mesh_type, typename matrix_type, typename vector_type>
ADRAssembler< mesh_type, matrix_type, vector_type>::
ADRAssembler() :

    M_fespace(),
    M_betaFESpace(),

    M_massCFE(),
    M_diffCFE(),
    M_advCFE(),
    M_advBetaCFE(),
    M_massRhsCFE(),

    M_localMass(),
    M_localAdv(),
    M_localDiff(),
    M_localMassRhs(),

    M_diffusionAssemblyChrono(),
    M_advectionAssemblyChrono(),
    M_massAssemblyChrono(),
    M_setupChrono(),
    M_massRhsAssemblyChrono()
{}

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
setup ( const fespace_ptrType& fespace, const fespace_ptrType& betaFESpace )
{
    ASSERT (fespace != 0, " Empty FE Space cannot setup the ADR assembler ");
    ASSERT (betaFESpace != 0, "Empty beta FE Space cannot setup the ADR assembler ");

    M_setupChrono.start();

    setFespace (fespace);
    setBetaFespace (betaFESpace);

    M_setupChrono.stop();
}

// ===================================================
// Methods
// ===================================================

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
addMass (matrix_ptrType matrix, const Real& coefficient, const UInt& offsetLeft, const UInt& offsetUp)
{
    // Check that the fespace is set
    ASSERT (M_fespace != 0, "No FE space for assembling the mass!");

    M_massAssemblyChrono.start();

    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );
    const UInt fieldDim (M_fespace->fieldDim() );
    const UInt nbTotalDof (M_fespace->dof().numTotalDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the mass current FE
        M_massCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_PHI | UPDATE_WDET );

        // Clean the local matrix
        M_localMass->zero();

        // Local Mass
        AssemblyElemental::mass (*M_localMass, *M_massCFE, coefficient, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < fieldDim; ++iFieldDim)
        {
            assembleMatrix ( *matrix,
                             *M_localMass,
                             *M_massCFE,
                             *M_massCFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * nbTotalDof + offsetLeft, iFieldDim * nbTotalDof + offsetUp);
        }
    }

    M_massAssemblyChrono.stop();
}


template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
addAdvection (matrix_ptrType matrix, const vector_type& beta, const UInt& offsetLeft, const UInt& offsetUp)
{
    // Beta has to be repeated!
    if (beta.mapType() == Unique)
    {
        addAdvection (matrix, vector_type (beta, Repeated), offsetLeft, offsetUp);
        return;
    }

    // Check that the fespace is set
    ASSERT (M_fespace != 0, "No FE space for assembling the advection!");
    ASSERT (M_betaFESpace != 0, "No FE space (beta) for assembling the advection!");

    M_advectionAssemblyChrono.start();


    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );
    const UInt fieldDim (M_fespace->fieldDim() );
    const UInt betaFieldDim (M_betaFESpace->fieldDim() );
    const UInt nbTotalDof (M_fespace->dof().numTotalDof() );
    const UInt nbQuadPt (M_advCFE->nbQuadPt() );

    // Temporaries
    //Real localValue(0);
    std::vector< std::vector< Real > > localBetaValue (nbQuadPt, std::vector<Real> ( betaFieldDim, 0.0 ) );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the advection current FEs
        M_advCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_PHI | UPDATE_DPHI | UPDATE_WDET );
        M_advBetaCFE->update (M_fespace->mesh()->element (iterElement), UPDATE_PHI );

        // Clean the local matrix
        M_localAdv->zero();

        // Interpolate beta in the quadrature points
        AssemblyElemental::interpolate (localBetaValue, *M_advBetaCFE, betaFieldDim, M_betaFESpace->dof(), iterElement, beta);

        // Assemble the advection
        AssemblyElemental::advection (*M_localAdv, *M_advCFE, 1.0, localBetaValue, fieldDim);


        // Assembly
        for (UInt iFieldDim (0); iFieldDim < fieldDim; ++iFieldDim)
        {
            assembleMatrix ( *matrix,
                             *M_localAdv,
                             *M_advCFE,
                             *M_advCFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * nbTotalDof + offsetLeft, iFieldDim * nbTotalDof + offsetUp );
        }
    }

    M_advectionAssemblyChrono.stop();
}

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
addDiffusion (matrix_ptrType matrix, const Real& coefficient, const UInt& offsetLeft, const UInt& offsetUp)
{
    // Check that the fespace is set
    ASSERT (M_fespace != 0, "No FE space for assembling the diffusion!");

    M_diffusionAssemblyChrono.start();

    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );
    const UInt fieldDim (M_fespace->fieldDim() );
    const UInt nbTotalDof (M_fespace->dof().numTotalDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_diffCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local matrix
        M_localDiff->zero();

        // local stiffness
        AssemblyElemental::stiffness (*M_localDiff, *M_diffCFE, coefficient, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < fieldDim; ++iFieldDim)
        {
            assembleMatrix ( *matrix,
                             *M_localDiff,
                             *M_diffCFE,
                             *M_diffCFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * nbTotalDof + offsetLeft, iFieldDim * nbTotalDof + offsetUp );
        }
    }

    M_diffusionAssemblyChrono.stop();
}

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
addMassRhs (vector_type& rhs, const vector_type& f)
{
    // f has to be repeated!
    if (f.mapType() == Unique)
    {
        addMassRhs (rhs, vector_type (f, Repeated) );
        return;
    }

    // Check that the fespace is set
    ASSERT (M_fespace != 0, "No FE space for assembling the right hand side (mass)!");

    M_massRhsAssemblyChrono.start();

    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );
    const UInt fieldDim (M_fespace->fieldDim() );
    const UInt nbFEDof (M_massRhsCFE->nbFEDof() );
    const UInt nbQuadPt (M_massRhsCFE->nbQuadPt() );
    const UInt nbTotalDof (M_fespace->dof().numTotalDof() );

    // Temporaries
    Real localValue (0.0);
    std::vector<Real> fValues (nbQuadPt, 0.0);

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_massRhsCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_PHI | UPDATE_WDET );

        // Clean the local matrix
        M_localMassRhs->zero();

        // Assemble the local diffusion
        for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
        {
            localVector_type::vector_view localView = M_localMassRhs->block (iterFDim);

            // Compute the value of f in the quadrature nodes
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                fValues[iQuadPt] = 0.0;
                for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
                {
                    fValues[iQuadPt] +=
                        f[ M_fespace->dof().localToGlobalMap (iterElement, iDof) + iterFDim * nbTotalDof]
                        * M_massRhsCFE->phi (iDof, iQuadPt);
                }
            }

            // Loop over the basis functions
            for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
            {
                localValue = 0.0;

                //Loop on the quadrature nodes
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                {
                    localValue += fValues[iQuadPt]
                                  * M_massRhsCFE->phi (iDof, iQuadPt)
                                  * M_massRhsCFE->wDetJacobian (iQuadPt);
                }

                // Add on the local matrix
                localView (iDof) = localValue;
            }
        }

        // Here add in the global rhs
        for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
        {
            assembleVector ( rhs,
                             iterElement,
                             *M_localMassRhs,
                             nbFEDof,
                             M_fespace->dof(),
                             iterFDim,
                             iterFDim * M_fespace->dof().numTotalDof() );
        }
    }

    M_massRhsAssemblyChrono.stop();
}

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
addMassRhs (vector_type& rhs, const function_type& f, const Real& t)
{
    // Check that the fespace is set
    ASSERT (M_fespace != 0, "No FE space for assembling the right hand side (mass)!");

    M_massRhsAssemblyChrono.start();

    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );
    const UInt fieldDim (M_fespace->fieldDim() );
    const UInt nbFEDof (M_massRhsCFE->nbFEDof() );
    const UInt nbQuadPt (M_massRhsCFE->nbQuadPt() );

    // Temporaries
    Real localValue (0.0);
    std::vector<Real> fValues (nbQuadPt, 0.0);

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_massRhsCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_QUAD_NODES | UPDATE_PHI | UPDATE_WDET );

        // Clean the local matrix
        M_localMassRhs->zero();

        // Assemble the local diffusion
        for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
        {
            localVector_type::vector_view localView = M_localMassRhs->block (iterFDim);

            // Compute the value of f in the quadrature nodes
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                fValues[iQuadPt] = f (t,
                                      M_massRhsCFE->quadNode (iQuadPt, 0),
                                      M_massRhsCFE->quadNode (iQuadPt, 1),
                                      M_massRhsCFE->quadNode (iQuadPt, 2),
                                      iterFDim);
            }

            // Loop over the basis functions
            for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
            {
                localValue = 0.0;

                //Loop on the quadrature nodes
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
                {
                    localValue += fValues[iQuadPt]
                                  * M_massRhsCFE->phi (iDof, iQuadPt)
                                  * M_massRhsCFE->wDetJacobian (iQuadPt);
                }

                // Add on the local matrix
                localView (iDof) = localValue;
            }
        }

        // Here add in the global rhs
        for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
        {
            assembleVector ( rhs,
                             iterElement,
                             *M_localMassRhs,
                             nbFEDof,
                             M_fespace->dof(),
                             iterFDim,
                             iterFDim * M_fespace->dof().numTotalDof() );
        }
    }

    M_massRhsAssemblyChrono.stop();
}


// ===================================================
// Set Methods
// ===================================================

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
setFespace (const fespace_ptrType& fespace)
{
    ASSERT (fespace != 0, " Setting the FE space for the unknown to 0 is not permitted ");

    M_fespace = fespace;

    M_massCFE.reset (new currentFE_type (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_localMass.reset (new localMatrix_type (M_fespace->fe().nbFEDof(),
                                             M_fespace->fieldDim(),
                                             M_fespace->fieldDim() ) );

    M_advCFE.reset (new currentFE_type (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_localAdv.reset (new localMatrix_type (M_fespace->fe().nbFEDof(),
                                            M_fespace->fieldDim(),
                                            M_fespace->fieldDim() ) );

    M_diffCFE.reset (new currentFE_type (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_localDiff.reset (new localMatrix_type (M_fespace->fe().nbFEDof(),
                                             M_fespace->fieldDim(),
                                             M_fespace->fieldDim() ) );

    M_massRhsCFE.reset (new currentFE_type (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_localMassRhs.reset (new localVector_type (M_fespace->fe().nbFEDof(), M_fespace->fieldDim() ) );
}

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
setBetaFespace (const fespace_ptrType& betaFESpace)
{
    ASSERT (M_fespace != 0, " No FE space for the unknown! Use setFespace before setBetaFespace!");
    ASSERT (M_advCFE != 0, " No current FE set for the advection of the unknown! Internal error.");
    M_betaFESpace = betaFESpace;

    M_advBetaCFE.reset (new currentFE_type (M_betaFESpace->refFE(), M_fespace->fe().geoMap(), M_advCFE->quadRule() ) );
}


} // Namespace LifeV

#endif /* ADRASSEMBLER_H */
