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
    @brief This file contains the definition of the LocalDofPattern class.

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef _LOCAL_DOF_PATTERN_HH
#define _LOCAL_DOF_PATTERN_HH

#include <life/lifecore/life.hpp>

#include <utility>

namespace LifeV
{

//! Local pattern type
/*!
  This enum allows to distinguish the normal standard local pattern, which is a full pattern involving all
  degrees of freedom to special patterns. It is stored in LocalDofPattern for later use by Dof
*/

enum DofPatternType {STANDARD_PATTERN = 1,
                     P1ISOP2_SEG_PATTERN = 2,
                     P1ISOP2_TRIA_PATTERN = 3
                    };

//! LocalDofPattern - A class to store the "couplings" between the basis functions
/*!
  The aim of this class is to store the way the basis functions couple one with each other. This might seem useless,
  however, some "advanced" finite elements require this structure.

  For example, consider the P1-iso-P2 element in 2D. This finite element is composed of 6 basis functions, based on the
  nodes with the same numerotation as the P2 element. The reference triangle is split into 4 subtriangles using the
  nodes on the faces. Each basis function is build such that it is 1 on its node, 0 on all the other nodes and such that
  it is linear in each subtriangle.

  @see in "Numerical Approximation of Partial Differential Equations" by A. Quarteroni and A. Valli, p.311, for an
  illustration and further informations.

  With this definition of the P1-iso-P2 finite element, we see that the basis functions 1 and 2 have no common support.
  So, they are not directly coupled.

  In order to represent the couplings between the basis functions, we can use a matrix \f$ C \f$ such that
  \f$ C_{ij} = 1 \f$ if the basis functions \f$ i \f$ and \f$ j \f$ have a common support, otherwise it is \f$ 0 \f$.
  This matrix is symmetric.

  For the P1-iso-P2 element, this matrix would be:
  \f[
  \begin{array}{|c||cccccc|}
  \hline
   & 1 & 2 & 3 & 4 & 5 & 6 \\
  \hline
  \hline
  1 & 1 & 0 & 0 & 1 & 0 & 1 \\
  2 & 0 & 1 & 0 & 1 & 1 & 0 \\
  3 & 0 & 0 & 1 & 0 & 1 & 1 \\
  4 & 1 & 1 & 0 & 1 & 1 & 1 \\
  5 & 0 & 1 & 1 & 1 & 1 & 1 \\
  6 & 1 & 0 & 1 & 1 & 1 & 1 \\
  \hline
  \end{array}
  \f]

  When references to diagonal or upper part are made, it is with respect to this matrix. For most of the finite
  elements, this matrix is full (lagrangian FE, lagrangian FE with bubbles,...).

  Instead of this representation with a matrix, we usually prefer to get the DoFs that are coupled in a list. This
  is the implemented in this class and the list of pairs (patternFirst(i),patternSecond(i)) representes all the DoFs
  that are coupled.


  Note: The documentation of this class (and some improvements) has been done by Samuel Quinodoz (15.01.2010),
  but its original implementation was prior the documentation and no name of author or date was available.
 */


class LocalDofPattern
{

public:

    //! @name Constructor & Destructor
    //@{

    //! Full constructor for 3D elements
    LocalDofPattern( const UInt& nbLocalDof, const UInt& nbDofPerVertex, const UInt& nbDofPerEdge,
                     const UInt& nbDofPerFace, const UInt& nbDofPerVolume, const DofPatternType& patternType );

    //! Full constructor for 2D elements
    LocalDofPattern( const UInt& nbLocalDof, const UInt& nbDofPerVertex, const UInt& nbDofPerEdge,
                     const UInt& nbDofPerFace, const DofPatternType& patternType );

    //! Full constructor for 1D elements
    LocalDofPattern( const UInt& nbLocalDof, const UInt& nbDofPerVertex, const UInt& nbDofPerEdge,
                     const DofPatternType& patternType );

    //! Simple copy constructor
    LocalDofPattern( const LocalDofPattern& localDofPattern);

    //! Empty destructor
    ~LocalDofPattern()
    {};

    //@}


    //! @name Methods
    //@{

    //!  patternFirst(i): row index in the element matrix of the i-th term of the pattern (the index starts from 0, not from 1 !).
    const UInt& patternFirst(const UInt& i ) const
    {
        ASSERT_BD( i < M_nbPattern );
        return M_pattern[i].first;
    }

    //! patternSecond(i): column index in the element matrix of the i-th term of the pattern (the index starts from 0, not from 1 !).
    const UInt& patternSecond(const UInt& i ) const
    {
        ASSERT_BD( i < M_nbPattern );
        return M_pattern[i].second;
    }

    //! The showMe method for the pattern
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Get Methods
    //@{

    //! Number of non-zero terms in the element matrix
    const UInt& nbPattern() const
    {
        return M_nbPattern;
    }

    //! Number of diagonal terms in the element matrix
    const UInt& nbDiag() const
    {
        return M_nbDiag;
    }

    //! Number of upper terms in the element matrix
    const UInt& nbUpper() const
    {
        return M_nbUpper;
    }

    //! Return the number of local degrees of freedom
    const UInt& nbLocalDof() const
    {
        return M_nbLocalDof;
    };

    //! Return the number of degrees of freedom located on the vertices (0D structures)
    const UInt& nbDofPerVertex() const
    {
        return M_nbDofPerDimEntity[0];
    };

    //! Return the number of degrees of freedom located on the edges (1D structures)
    const UInt& nbDofPerEdge() const
    {
        ASSERT(M_dim >=1, "No edge available for that dimension");
        return M_nbDofPerDimEntity[1];
    };

    //! Return the number of degrees of freedom located on the faces (2D structures).
    /*!Beware that in the 2D case, the face of a triangle is the triangle itself
      (use edges or vertices if you want to access substructures).
     */
    const UInt& nbDofPerFace() const
    {
        ASSERT(M_dim >=2, "No face available for that dimension");
        return M_nbDofPerDimEntity[2];
    };

    //! Return the number of degrees of freedom located in the volume (3D structures).
    const UInt& nbDofPerVolume() const
    {
        ASSERT(M_dim >=3, "No volume available for that dimension");
        return M_nbDofPerDimEntity[3];
    };

    //! Return the number of degrees of freedom located per structDim object.
    /*! For example, if we want to access the vertices, structDim should be 0,
      if we want the edges, then it should be 1,...
     */
    const UInt& nbDofPerDimStrut(const UInt & structDim) const
    {
        ASSERT(structDim <= M_dim, "No structure with this dimension");
        return M_nbDofPerDimEntity[structDim];
    };

    //! Return the number of degrees of freedom located per structCodim object.
    /*! The codimension of a structure is the full dimension of the element
      minus the dimension of the structure. For example, in 3D, faces have codimension
      1, edges 2 and vertices 3. This method could be usefull to code "dimension-free" code.
      (for example, IP is built on edges in 2D, faces in 3D, so on objects with codimension 1).
     */
    const UInt& nbDofPerCodimStrut(const UInt & structCodim) const
    {
        ASSERT(structCodim <= M_dim, "No structure with this codimension");
        return M_nbDofPerDimEntity[M_dim-structCodim];
    };


    //@}


private:

    //! @name Private Methods
    //@{

    //! Default constructor disabled (because there is no setup/set method)
    LocalDofPattern();

    //! Method to setup the standard pattern, i.e. with all degrees of freedom coupled.
    void setupStandardPattern();

    //! Method for the P1isoP2 pattern for the segments (1D)
    void setupP1isoP2SegPattern();

    //! Method for the P1isoP2 pattern for the triangles (2D)
    void setupP1isoP2TriaPattern();

    //@}


    //! dimension of the element (3 for a tetrahedra for example).
    UInt M_dim;

    //! Total number of degrees of freedom (equal to refEle::nbDof)
    UInt M_nbLocalDof;

    //! Number of degrees of freedom per geometric entity
    /*! In this vector, we store all the number of degrees of freedom
      sorted by the dimension of the object where they lie.
      For example, the number of DoF per vertex (dimension 0)
      is stored as the 0th element. Then, for edges (dimension 1), it is
      stored in the 1st element of the vector,...
      This enables a n-dimensional implementation (not only 3D)
    */
    std::vector< UInt> M_nbDofPerDimEntity;


    //! Type of the pattern stored
    DofPatternType M_patternType;

    //! Pairs of couplings to appear in the pattern
    std::vector< std::pair< UInt,UInt > > M_pattern;

    //! Number of non-zero terms in the element matrix
    UInt M_nbPattern;

    //! Number of diagonal terms in the element matrix
    UInt M_nbDiag;

    //! Number of upper terms in the element matrix
    UInt M_nbUpper;

};


}
#endif

