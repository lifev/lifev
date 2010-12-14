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
    @brief File containing the CurrentFE class

    @author Jean-Frederic Gerbeau
    @date 00-04-2002

    @contributor V. Martin
                 Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef CURRENTFE_H
#define CURRENTFE_H 1


#include <life/lifecore/life.hpp>

#include <life/lifefem/geoMap.hpp>

#include <life/lifefem/refFEScalar.hpp>
#include <life/lifefem/refFEHdiv.hpp>
#include <life/lifefem/refFEHybrid.hpp>

#include <life/lifefem/quadRule.hpp>

#include <boost/multi_array.hpp>

#include <iostream>
#include <fstream>

namespace LifeV
{

/*! \page update_procedure Update of the current finite element

  \section general_section General issues

  When the update method of the LifeV::currentFE class is called, several arrays of values are updated with respect to the current cell (and the current quadrature). However, in order to avoid redundant computations, intermediate values are also computed. The next figure shows all the possible values that can be updated.

  \image html update_complete_scheme.png "Complete tree of the values. The arrow between two values means that the first value require the second value to be already computed."

  However, one could want only a part of the values to be updated (we maybe do not want the second derivative for a simple elliptic problem). This is why a set of flag has been set up. With this flags, one selects the values that he needs (represented on the top of the last figure) and all the needed values are updated, without redundancy. For example, if one calles the update method with the flag UPDATE_DPHI, only the values visible on the next figure will be updated.

  \image html update_dphi_scheme.png "Values updated when only the derivatives of the basis functions are required."


  Of course, it is possible to combine the different flags to update several values in the same time.

  \section flags_section Flags explaination

  In this section, we explain how flags work. This can be of interest if one wants to add a new flag.

  \subsection flag_sub1 What is a flag?

  At first sight, we might think that a flag is a complicated class implemented to do exaclty what we want, with overloaded operators... Actually, it is much simpler: a LifeV::flag_Type is just an unsigned integer.

  \subsection How to define a flag?

  The flags use the binary representation of the integers to work. This enables a very fast definition and use of the flags. To understand it, let us make a simple example. Suppose that we can update three quantities A,B and C.

  The first step is to define a "primitive" flag for each of these quantities. These flags are defined as powers of 2. Here, we will define

  \code
  flag_Type UPDATE_A(1);
  flag_Type UPDATE_B(2);
  flag_Type UPDATE_C(4);
  \endcode

  We need here powers of 2 because this makes the binary representation of the "primitive" flags simple: UPDATE_A is 001, UPDATE_B is 010 and UPDATE_C is 100 (if the integers are coded on 3 bits, otherwise there are many zeros before). The fact that we want A to be update is then represented by a 1 in the third position, for B it is the second position and for C the first position.

  So now, if we want to build a flag that updates A and C, we want it to be in binary 101, that is 5. So, the flag UPDATE_AC will be defined as:

  \code
  flag_Type UPDATE_AC(5);
  \endcode

  \subsection flag_sub2 How to combine flags?

  With the last example, one could think that is it just a matter of addition, but it is not. Suppose that we want to combine the flags UPDATE_A and UPDATE_AC. If we add the corresponding values, we will get 1+5=6 that is represented in binary by 110. This would mean update B and C, what is not what we want!

  To combine flags, we use the binary operator | that do exactly the job that we want: 1|5=5.

  \subsection flag_sub3 How to detect flags?

  Now that we can build every flag that we want, we want to detect quickly the different flags. This is achievied by using the binary operator & and the "primitive" flags.

  Suppose that we want to know if A has to be updated. Then, we perform that operation "& UPDATE_A" on the incoming flag. If the result is zero, then we do not need to update it: for example, UPDATE_B & UPDATE_A = 0 . Otherwise, we have to update it: for example, UPDATE_AC & UPDATE_A = 1.


  \section list_section Possible flags

  There are several flags defined for you. Here is a list of the possible flags that are usually used:

  <table>
   <tr>
    <th> Flag name </th> <th> Effect </th>
   </tr>
   <tr>
    <th> UPDATE_QUAD_NODES </th> <th> Update everything needed to know the position of the quadrature nodes in the current cell </th>
   </tr>
   <tr>
    <th> UPDATE_PHI </th> <th> Update everything needed to know the values of the basis functions in the quadrature nodes (in the current cell) </th>
   </tr>
   <tr>
    <th> UPDATE_DPHI </th> <th>  Update everything needed to know the values of the derivatives of the basis functions in the quadrature nodes (in the current cell) </th>
   </tr>
   <tr>
    <th> UPDATE_D2PHI </th> <th>  Update everything needed to know the values of the second derivatives of the basis functions in the quadrature nodes (in the current cell) </th>
   </tr>
   <tr>
    <th> UPDATE_WDET </th> <th>  Update everything needed to know the determinant of the transformation multiplied by the weights of the quadrature. </th>
   </tr>

  </table>

  Besides this usual flags, there are a couple of "primitive" flags, that update only a particular element in the currentFE structure. Be sure to know what you are doing before using them.

*/



// This is a 32 bits type, so we can store up to 32 different flags
typedef unsigned int flag_Type;

const flag_Type UPDATE_ONLY_CELL_NODES(1);
const flag_Type UPDATE_ONLY_QUAD_NODES(2);
const flag_Type UPDATE_ONLY_PHI(4);
const flag_Type UPDATE_ONLY_DPHI_GEO_MAP(8);
const flag_Type UPDATE_ONLY_JACOBIAN(16);
const flag_Type UPDATE_ONLY_T_INVERSE_JACOBIAN(32);
const flag_Type UPDATE_ONLY_W_DET_JACOBIAN(64);
const flag_Type UPDATE_ONLY_DPHI_REF(128);
const flag_Type UPDATE_ONLY_DPHI(256);
const flag_Type UPDATE_ONLY_D2PHI_REF(512);
const flag_Type UPDATE_ONLY_D2PHI(1024);
const flag_Type UPDATE_ONLY_PHI_VECT(2048);
const flag_Type UPDATE_ONLY_DIV_PHI_REF(4096);
const flag_Type UPDATE_ONLY_DET_JACOBIAN(8192);


const flag_Type UPDATE_QUAD_NODES(UPDATE_ONLY_CELL_NODES
                                  |UPDATE_ONLY_QUAD_NODES);

const flag_Type UPDATE_PHI(UPDATE_ONLY_PHI);

const flag_Type UPDATE_DPHI(UPDATE_ONLY_CELL_NODES
                            |UPDATE_ONLY_DPHI_GEO_MAP
                            |UPDATE_ONLY_JACOBIAN
                            |UPDATE_ONLY_T_INVERSE_JACOBIAN
                            |UPDATE_ONLY_DPHI_REF
                            |UPDATE_ONLY_DPHI);

const flag_Type UPDATE_D2PHI(UPDATE_ONLY_CELL_NODES
                             |UPDATE_ONLY_DPHI_GEO_MAP
                             |UPDATE_ONLY_JACOBIAN
                             |UPDATE_ONLY_T_INVERSE_JACOBIAN
                             |UPDATE_ONLY_D2PHI_REF
                             |UPDATE_ONLY_D2PHI);

const flag_Type UPDATE_WDET(UPDATE_ONLY_CELL_NODES
                            |UPDATE_ONLY_DPHI_GEO_MAP
                            |UPDATE_ONLY_JACOBIAN
                            |UPDATE_ONLY_DET_JACOBIAN
                            |UPDATE_ONLY_W_DET_JACOBIAN);

const flag_Type UPDATE_PHI_VECT(UPDATE_ONLY_CELL_NODES
                                |UPDATE_ONLY_DPHI_GEO_MAP
                                |UPDATE_ONLY_JACOBIAN
                                |UPDATE_ONLY_DET_JACOBIAN
                                |UPDATE_ONLY_PHI_VECT);

const flag_Type UPDATE_DIV_PHI(UPDATE_ONLY_DIV_PHI_REF);



//! currentFE - A primordial class for the assembly of the local matrices/vectors retaining the values on the real cells.
/*!
  This class is an essential piece of any finite element solver implemented in LifeV, as it is used each time one wants to assemble the matrix or the right hand side associated to a given discrete problem (in variational form). The role of the class is, given a cell (or similar data), to compute the values of the basis function, their derivatives, the location of the quadrature nodes,...

  During the assembly procedure, the code of the solver will typically loops on the cells of a given mesh. Each time a cell is selected, a local matrix is built (using only the degrees of freedom of that cell). This local matrix represents a local contribution that is added in the global matrix (concerning all the degrees of freedom of the problem).

  <b>Example</b>: suppose that we are assembling the local matrix of a Laplacian problem. The local contributions are given by \f$ \int_{K} \nabla \phi_i \cdot \nabla \phi_j \f$ where \f$K\f$ is the considered cell. As we use numerical quadrature, to build the local matrix, we have to provide:

  <ol>
  <li> The values of the <b>derivatives</b> of the basis functions in the quadrature nodes
  <li> The <b>weights</b> associated to the quadrature and adapted to the considered cell
  </ol>

  The CurrentFE has then to compute these values. To precise that we want to be able to access these values, we use flags. Here we need the two flags (see \ref update_procedure "this page" for more informations on the flags)

  <ol>
  <li> For the derivative of the basis functions : UPDATE_DPHI
  <li> For the modified weights : UPDATE_WDET
  </ol>

  The code that we would typically use is:

  \code
  my_currentFE.update(my_mesh.volumeList(i), UPDATE_DPHI | UPDATE_WDET);
  \endcode

  This line of code asks the CurrentFE to compute the values specified with the flags and corresponding to the cell given by the i-th volume of the desired mesh. Now, these values are accessible using simple calls:

  \code
  Real my_value = my_currentFE.dphi(0,0,1);
  Real my_weight = my_currentFE.wDetJacobian(1);
  \endcode

  These two lines just get some values already stored in the CurrentFE (so it is quite cheap!). Usually, one here uses the elementary operation already defined in the file lifefem/elemOper.hpp.


  <b>Note</b>: One can change the quadrature of the CurrentFE at any moment, but this is a quite expensive operation and all updates previously done are lost. Use this feature only when it is really needed!

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 13 Apr 2010
    @version 2.0

    \todo Put all the remaining public member in private
    \todo Put ifdef for the checks (definition of the booleans)
    \todo Change the old update for wrappers to the new update
    \todo CXXFLAG to disable the boost asserts: -DBOOST_DISABLE_ASSERTS?

*/
class CurrentFE
{

public:

    //! @name Constructor & Destructor
    //@{

    //! Full constructor with default quadrature
    /*!
      @param refFE Reference finite element used
      @param geoMap Geometric mapping used
      @param qr Quadrature rule for the computations
     */
    CurrentFE( const RefFE& refFE, const GeoMap& geoMap, const QuadRule& qr );

    //! Constructor without quadrature specification
    /*!
      If you use this constructor, you have to use the method CurrentFE::setQuadRule before
      using any update method!

      @param refFE Reference finite element used
      @param geoMap Geometric mapping used
     */
    CurrentFE( const RefFE& refFE, const GeoMap& geoMap);

    //! Destructor
    ~CurrentFE() { delete M_quadRule;}

    //@}


    //! @name Methods
    //@{

    //! Update method using a given element.
    /*!
      This method is the most important one in this class. It allows the user to update different values on
      the selected cell, like the values of the derivatives of the basis functions,...
      @param geoele The cell that we are looking at
      @param upFlag The flag to explain the quantities that we want to update
     */
    template<typename GeoElementType>
    void update(const GeoElementType& geoele, const flag_Type& upFlag);

    //! Update method using only point coordinates. It used the flags, as defined in \ref update_procedure "this page".
    void update(const std::vector<std::vector<Real> >& pts, const flag_Type& upFlag);

    //! Update method using only point coordinates. It used the flags, as defined in \ref update_procedure "this page".
    void update(const std::vector<GeoVector>& pts, const flag_Type& upFlag);

    //! Return the measure of the current element
    Real measure() const;

    //! return the barycenter of the element
    void barycenter( Real& x, Real& y, Real& z ) const;

    //! return the diameter of the element in the 1-norm
    Real diameter() const;

    //! return the diameter of the element in the 2-norm
    Real diameter2() const;

    /*!
      compute the coordinate (x,y,z)= F(xi,eta,zeta), (F: geo mappping)
      where (xi,eta,zeta) are the coor in the ref element
      and   (x,y,z) the coor in the current element
      (if the code is compiled in 2D mode then z=0 and zeta is disgarded)
    */
    void coorMap( Real& x, Real& y, Real& z, const Real & xi, const Real & eta, const Real & zeta ) const;

    /*!
      return the coordinates in the current element of the point
      P given in the reference frame.
    */
    GeoVector coorMap(const GeoVector& P) const;

    //! Export the quadrature rule on the current FE
    /*!
      This method can be used to position of the quadrature nodes in the
      considered element using VTK format. The quadrature point must be updated.
      It differs from the quadRule::VTKexport in the sens that here, the quadrature
      is mapped on the current element, while it is still in the reference element
      for quadRule::VTKexport.
     */
    void quadRuleVTKexport( const std::string& filename) const;

    void __attribute__ ((__deprecated__)) QuadRuleVTKexport( const std::string& filename) const
    {
        return quadRuleVTKexport(filename);
    }


    //@}


    //! @name Set Methods
    //@{

    //! Setter for the quadrature rule
    /*! This method can be used to change the quadrature rule used. Notice that it is an
      expensive method. Use it only when it is really needed.
      @param newQuadRule The new quadrature rule
     */
    void setQuadRule(const QuadRule& newQuadRule);

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the ID of the current cell
    inline const UInt& currentId() const
    {
        return M_currentId;
    }

    //! Getter for the local ID of the current cell
    inline const UInt& currentLocalId() const
    {
        return M_currentLocalId;
    }

    //! Getter for the number of quadrature nodes
    inline const UInt& nbQuadPt() const
    {
        return M_nbQuadPt;
    }

    //! Getter for the number of nodes
    inline const UInt& nbFEDof() const
    {
        return nbNode;
    }

    //! Getter for the reference FE
    inline const RefFE& refFE() const
    {
        return *M_refFE;
    };

    //! Getter for the GeoMap reference
    inline const GeoMap& geoMap() const
    {
        return *M_geoMap;
    }

    //! Getter for the quadrature rule
    inline const QuadRule& quadRule() const
    {
        return *M_quadRule;
    };

    //! Getter for the number of entries in the pattern
    inline const UInt& nbPattern() const
    {
        return M_nbPattern;
    };

    //! Getter for the diagonal entries in the pattern
    inline const UInt& nbDiag() const
    {
        return M_nbDiag;
    };

    //! Getter for the number of entries in the pattern
    inline const UInt& nbUpper() const
    {
        return M_nbUpper;
    };

    //! Getter for the number of coordinates
    inline const UInt& nbCoor() const
    {
        return M_nbCoor;
    };

    //@}


    //! @name Get Methods for the FE values
    //@{

    //! Getter for the nodes of the cell
    inline const Real& cellNode(const UInt& node, const UInt& coordinate) const
    {
        ASSERT(M_cellNodesUpdated,"Cell nodes are not updated!");
        return M_cellNodes[node][coordinate];
    };

    //! Getter for the quadrature nodes
    inline const Real& quadNode(const UInt& node, const UInt& coordinate) const
    {
        ASSERT(M_quadNodesUpdated,"Quad nodes are not updated!");
        return M_quadNodes[node][coordinate];
    };

    //! Getter for the determinant of the jacobian in a given quadrature node
    inline const Real& detJacobian(const UInt& quadNode) const
    {
        ASSERT(M_detJacobianUpdated,"Jacobian determinant is not updated!");
        return M_detJacobian[quadNode];
    };

    //! Getter for the weighted jacobian determinant
    inline const Real& wDetJacobian(const UInt& quadNode) const
    {
        ASSERT(M_wDetJacobianUpdated,"Weighted jacobian determinant is not updated!");
        return M_wDetJacobian[quadNode];
    };

    //! Getter for the inverse of the jacobian
    inline const Real& tInverseJacobian(const UInt& element1, const UInt& element2, const UInt& quadNode) const
    {
        ASSERT(M_tInverseJacobianUpdated,"Inverse jacobian is not updated!");
        return M_tInverseJacobian[element1][element2][quadNode];
    };

    //! Getter for basis function values (scalar FE case)
    inline const Real& phi(const UInt& node, const UInt& quadNode) const
    {
        ASSERT(M_phiUpdated,"Function values are not updated!");
        return M_phi[node][0][quadNode];
    };

    //! Getter for basis function values (vectorial FE case)
    inline const Real& phi(const UInt& node, const UInt& component, const UInt& quadNode) const
    {
        ASSERT(M_phiVectUpdated,"Function values are not updated!");
        return M_phiVect[node][component][quadNode];
    };

    //! Getter for the derivatives of the basis functions
    inline const Real& dphi(const UInt& node, const UInt& derivative, const UInt& quadNode) const
    {
        ASSERT(M_dphiUpdated,"Basis derivatives are not updated!");
        return M_dphi[node][derivative][quadNode];
    };

    //! Getter for the second derivatives of the basis functions
    inline const Real& d2phi(const UInt& node, const UInt& derivative1, const UInt& derivative2, const UInt& quadNode) const
    {
        ASSERT(M_d2phiUpdated,"Basis second derivatives are not updated!");
        return M_d2phi[node][derivative1][derivative2][quadNode];
    };

    //! Getter for the divergence of a vectorial FE in the reference frame.
    inline const Real& divPhiRef(const UInt& node, const UInt& quadNode) const
    {
        ASSERT(M_divPhiRefUpdated,"Basis divergence are not updated!");
        return M_divPhiRef[node][quadNode];
    };

    //@}



    //! @name Old methods (for backward compatibility, avoid using them)
    //@{

    //! Old accessor, use cellNode instead.
    inline const Real& point(const UInt& node, const UInt& coordinate) const
    {
        ASSERT(M_cellNodesUpdated,"Cell nodes are not updated!");
        return M_cellNodes[node][coordinate];
    };

    //! Old accessor, use quadNode instead
    inline const Real& quadPt(const UInt& node, const UInt& coordinate) const
    {
        ASSERT(M_quadNodesUpdated,"Quad nodes are not updated!");
        return M_quadNodes[node][coordinate];
    };

    //! Old accessor, use wDetJacobian instead
    inline const Real& weightDet(const UInt& quadNode) const
    {
        ASSERT(M_wDetJacobianUpdated,"Weighted jacobian determinant is not updated!");
        return M_wDetJacobian[quadNode];
    };

    //! Getter for the determinant of the jacobian in a given quadrature node
    inline const Real& detJac(const UInt& quadNode) const
    {
        ASSERT(M_detJacobianUpdated,"Jacobian determinant is not updated!");
        return M_detJacobian[quadNode];
    };

    //! Old accessor, use iInverseJacobian instead
    inline const Real& tInvJac(const UInt& element1, const UInt& element2, const UInt& quadNode) const
    {
        ASSERT(M_tInverseJacobianUpdated,"Inverse jacobian is not updated!");
        return M_tInverseJacobian[element1][element2][quadNode];
    };

    //! Old accessor, use dphi instead
    inline const Real& phiDer(const UInt& node, const UInt& derivative, const UInt& quadNode) const
    {
        ASSERT(M_dphiUpdated,"Basis derivatives are not updated!");
        return M_dphi[node][derivative][quadNode];
    };

    //! Old accessor, use d2phi instead
    inline const Real& phiDer2(const UInt& node, const UInt& derivative1, const UInt& derivative2, const UInt& quadNode) const
    {
        ASSERT(M_d2phiUpdated,"Basis second derivatives are not updated!");
        return M_d2phi[node][derivative1][derivative2][quadNode];
    };

    //@}




private:
    CurrentFE( );
    CurrentFE( const CurrentFE& );

    //! Update the nodes of the cell to the current one.
    template<typename GeoElementType>
    void computeCellNodes(const GeoElementType& geoele);

    //! Update only the nodes of the cells to the current one.
    void computeCellNodes(const std::vector<std::vector< Real> >& pts);

    //! Update the location of the quadrature in the current cell.
    void computeQuadNodes();

    //! Compute the values of the basis functions in the quadrature nodes
    void computePhi();

    //! Compute the values of the derivatives of the mapping in the quadrature nodes
    void computeDphiGeoMap();

    //! Compute the jacobian in the quadrature nodes
    void computeJacobian();

    //! Compute the transposed inverse of the jacobian in the quadrature nodes
    void computeTInverseJacobian();

    //! Compute the determinant of the jacobian
    void computeDetJacobian();

    //! Compute the determinant of the jacobian in the quadrature nodes
    void computeWDetJacobian();

    //! Compute the value of the derivatives in the reference frame
    void computeDphiRef();

    //! Compute the value of the derivatives in the current element
    void computeDphi();

    //! Compute the value of the second derivatives in the reference frame
    void computeD2phiRef();

    //! Compute the value of the derivatives in the current element
    void computeD2phi();

    //! Compute the value of the vectorial FE using Piola transform
    void computePhiVect();

    // Constants
public:
    UInt nbNode;
private:
    const UInt M_nbCoor;
    const UInt M_nbDiag;
    const UInt M_nbUpper;
    const UInt M_nbPattern;

    UInt M_currentId;
    UInt M_currentLocalId;

    UInt M_nbGeoNode;
    UInt M_nbQuadPt;

    // Important structures

    const RefFE* M_refFE;
    const GeoMap* M_geoMap;
    QuadRule* M_quadRule;


    // Internal storage for the data

    // Nodes of the current cell
    boost::multi_array<Real,2> M_cellNodes;
    boost::multi_array<Real,2> M_quadNodes;

    boost::multi_array<Real,3> M_dphiGeoMap;
    boost::multi_array<Real,3> M_jacobian;
    boost::multi_array<Real,1> M_detJacobian;
    boost::multi_array<Real,1> M_wDetJacobian;
    boost::multi_array<Real,3> M_tInverseJacobian;

    boost::multi_array<Real,3> M_phi;
    boost::multi_array<Real,3> M_dphi;
    boost::multi_array<Real,4> M_d2phi;
    boost::multi_array<Real,3> M_phiVect;

    // M_phiRef is useless because M_phi is already the same.
    boost::multi_array<Real,3> M_dphiRef;
    boost::multi_array<Real,4> M_d2phiRef;
    boost::multi_array<Real,2> M_divPhiRef;

    // Check
    bool M_cellNodesUpdated;
    bool M_quadNodesUpdated;

    bool M_dphiGeoMapUpdated;
    bool M_jacobianUpdated;
    bool M_detJacobianUpdated;
    bool M_wDetJacobianUpdated;
    bool M_tInverseJacobianUpdated;

    bool M_phiUpdated;
    bool M_dphiUpdated;
    bool M_d2phiUpdated;
    bool M_phiVectUpdated;

    bool M_dphiRefUpdated;
    bool M_d2phiRefUpdated;
    bool M_divPhiRefUpdated;

// OLD FUNCTIONS

public:

    //! @name Old methods (for backward compatibility, avoid using them)
    //@{

    /*!
      compute the coordinate (xi,eta,zeta)=inv(F)(x,y,z)
    */
    void coorBackMap(const Real& x, const Real& y, const Real& z,
                     Real & xi, Real & eta, Real& zeta) const;

    /*!
      compute the jacobian at a given point : d x_compx / d zeta_compzeta
    */
    Real pointJacobian(const Real& hat_x, const Real& hat_y, const Real& hat_z,
                       int compx, int compzeta) const;

    /*!
      compute the inverse jacobian
     */

    Real pointInverseJacobian(const Real& hat_x, const Real& hat_y, const Real& hat_z,
                              int compx, int compzeta) const;

    /*!
      compute the determinant of the Jacobian at a given point
     */
    Real pointDetJacobian(const Real& hat_x, const Real& hat_y, const Real& hat_z) const;

    /*!  return (x,y,z) = the global coordinates of the quadrature point ig
      in the current element. \warning this function is almost obsolete since if
      you call the function updateFirstDerivQuadPt rather than updateFirstDeriv
      (for example), the coordinates of the quadrature points have already
      been computed and may be obtained via quadPt(ig,icoor). This is usually
      much less expensive since it avoids many calls to coorQuadPt
    */
    inline void coorQuadPt( Real& x, Real& y, Real& z, const int ig ) const
    {
        coorMap( x, y, z, M_quadRule->quadPointCoor( ig, 0 ), M_quadRule->quadPointCoor( ig, 1 ),
                 M_quadRule->quadPointCoor( ig, 2 ) );
    }
    //!  patternFirst(i): row index in the element matrix of the i-th term of the pattern
    inline int patternFirst( int i ) const
    {
        return M_refFE->patternFirst( i );
    }
    //! patternSecond(i): column index in the element matrix of the i-th term of the pattern
    inline int patternSecond( int i ) const
    {
        return M_refFE->patternSecond( i );
    }


    //---------------------------------------
    //! DECLARATION of the update methods
    //---------------------------------------
    /*!
      minimal update: we just identify the id of the current element
    */
    template <class GeoElementType>
    void update( const GeoElementType& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian on
      the current element
    */
    template <class GeoElementType>
    void updateJac( const GeoElementType& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian and quadPt
      on the current element
    */
    template <class GeoElementType>
    void updateJacQuadPt( const GeoElementType& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer on the current element
    */
    template <class GeoElementType>
    void updateFirstDeriv( const GeoElementType& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer and quadPt on the current element
    */
    template <class GeoElementType>
    void updateFirstDerivQuadPt( const GeoElementType& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer2 on the current element
    */
    template <class GeoElementType>
    void updateSecondDeriv( const GeoElementType& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer2 on the current element
    */
    template <class GeoElementType>
    void updateSecondDerivQuadPt( const GeoElementType& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer, phiDer2 on the current element
    */
    template <class GeoElementType>
    void updateFirstSecondDeriv( const GeoElementType& geoele );
    /*!
      compute the arrays detJac, weightDet, jacobian,
      tInvJac, phiDer, phiDer2 on the current element
    */
    template <class GeoElementType>
    void updateFirstSecondDerivQuadPt( const GeoElementType& geoele );

    //@}

};




template<typename GeoElement>
void CurrentFE::update(const GeoElement& geoele, const flag_Type& upFlag)
{
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();

    std::vector< std::vector <Real> > pts(M_nbGeoNode, std::vector<Real>(M_nbCoor));

    for ( UInt i(0); i < M_nbGeoNode; ++i )
    {
        for ( UInt icoor(0); icoor < M_nbCoor; ++icoor)
        {
            pts[i][icoor] = geoele.point(i+1).coordinate(icoor+1);
        }
    }
    update(pts,upFlag);
};


template<typename GeoElement>
void CurrentFE::computeCellNodes(const GeoElement& geoele)
{
    std::vector< std::vector <Real> > pts(M_nbGeoNode, std::vector<Real>(M_nbCoor));

    for ( UInt i(0); i < M_nbGeoNode; ++i )
    {
        for ( UInt icoor(0); icoor < M_nbCoor; ++icoor)
        {
            pts[i][icoor] = geoele.point(i+1).coordinate(icoor+1);
        }

    }
    computeCellNodes(pts);
};

//---------------------------------------
//! IMPLEMENTATION of the CurrentFE::update methods
//---------------------------------------

/*!
    minimal update: we just identify the id of the current element
  */
template <class GeoElementType>
void CurrentFE::update( const GeoElementType& geoele )
{
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
}

/*!
    compute the arrays detJac, weightDet, jacobian on
    the current element
*/
template <class GeoElementType>
void CurrentFE::updateJac( const GeoElementType& geoele )
{
    ASSERT(M_nbQuadPt!=0," No quadrature rule defined, cannot update!");
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
    //! compute the jacobian and its determinant...
    computeCellNodes(geoele);
    computeDphiGeoMap();
    computeJacobian();
    computeDetJacobian();
    computeWDetJacobian();
}

/*!
    compute the arrays detJac, weightDet, jacobian and quadPt
    on the current element
*/
template <class GeoElementType>
void CurrentFE::updateJacQuadPt( const GeoElementType& geoele )
{
    ASSERT(M_nbQuadPt!=0," No quadrature rule defined, cannot update!");
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
    // compute the jacobian and its determinant...
    computeCellNodes(geoele);
    computeDphiGeoMap();
    computeJacobian();
    computeDetJacobian();
    computeWDetJacobian();
    // and the coordinates of the quadrature points
    computeQuadNodes();
}

/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer on the current element
*/
template <class GeoElementType>
void CurrentFE::updateFirstDeriv( const GeoElementType& geoele )
{
    ASSERT(M_nbQuadPt!=0," No quadrature rule defined, cannot update!");
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
    //! compute the inverse jacobian...
    computeCellNodes(geoele);
    computeDphiGeoMap();
    computeJacobian();
    computeDetJacobian();
    computeWDetJacobian();
    computeTInverseJacobian();
    //! product InvJac by dPhiRef to compute phiDer
    computeDphiRef();
    computeDphi();
}

/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer and quadPt on the current element
*/
template <class GeoElementType>
void CurrentFE::updateFirstDerivQuadPt( const GeoElementType& geoele )
{
    ASSERT(M_nbQuadPt!=0," No quadrature rule defined, cannot update!");
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
    //! compute the inverse jacobian...
    computeCellNodes(geoele);
    computeDphiGeoMap();
    computeJacobian();
    computeDetJacobian();
    computeWDetJacobian();
    computeTInverseJacobian();
    //! product InvJac by dPhiRef to compute phiDer
    computeDphiRef();
    computeDphi();
    //! and the coordinates of the quadrature points
    computeQuadNodes();
}

// A. Veneziani, October 30, 2002
/*!
  compute the arrays detJac, weightDet, jacobian,
  tInvJac, phiDer2 on the current element
*/
template <class GeoElementType>
void CurrentFE::updateSecondDeriv( const GeoElementType& geoele )
{
    ASSERT(M_nbQuadPt!=0," No quadrature rule defined, cannot update!");
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
    //! compute the inverse jacobian...
    computeCellNodes(geoele);
    computeDphiGeoMap();
    computeJacobian();
    computeDetJacobian();
    computeWDetJacobian();
    computeTInverseJacobian();
    //! compute the second derivative
    computeD2phiRef();
    computeD2phi();
}

/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer2 on the current element
  */
template <class GeoElementType>
void CurrentFE::updateSecondDerivQuadPt( const GeoElementType& geoele )
{
    ASSERT(M_nbQuadPt!=0," No quadrature rule defined, cannot update!");
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
    //! compute the inverse jacobian...
    computeCellNodes(geoele);
    computeDphiGeoMap();
    computeJacobian();
    computeDetJacobian();
    computeWDetJacobian();
    computeTInverseJacobian();
    //! compute the second derivative
    computeD2phiRef();
    computeD2phi();
    //! and the coordinates of the quadrature points
    computeQuadNodes();
}
/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer, phiDer2 on the current element
  */
template <class GeoElementType>
void CurrentFE::updateFirstSecondDeriv( const GeoElementType& geoele )
{
    ASSERT(M_nbQuadPt!=0," No quadrature rule defined, cannot update!");
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
    //! compute the inverse jacobian...
    computeCellNodes(geoele);
    computeDphiGeoMap();
    computeJacobian();
    computeDetJacobian();
    computeWDetJacobian();
    computeTInverseJacobian();
    //! compute phiDer and phiDer2
    computeDphiRef();
    computeDphi();
    computeD2phiRef();
    computeD2phi();
}
/*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer, phiDer2 on the current element
  */
template <class GeoElementType>
void CurrentFE::updateFirstSecondDerivQuadPt( const GeoElementType& geoele )
{
    ASSERT(M_nbQuadPt!=0," No quadrature rule defined, cannot update!");
    M_currentId      = geoele.id();
    M_currentLocalId = geoele.localId();
    //! compute the inverse jacobian...
    computeCellNodes(geoele);
    computeDphiGeoMap();
    computeJacobian();
    computeDetJacobian();
    computeWDetJacobian();
    computeTInverseJacobian();
    //! compute phiDer and phiDer2
    computeDphiRef();
    computeD2phiRef();
    computeDphi();
    computeD2phi();
    //! and the coordinates of the quadrature points
    computeQuadNodes();
}

} // Namespace LifeV

#endif /* CURRENTFE_H */
