//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file BlockMatrix.hpp
    @brief file containing a class for handling a 2-blocks matrix with Robin-Neumann coupling

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 09 Jul 2010

 */

#ifndef BLOCKMATRIXRN_H
#define BLOCKMATRIXRN_H 1

#include <life/lifecore/life.hpp>
#include <lifemc/lifesolver/BlockMatrix.hpp>
#include <lifemc/lifesolver/RobinInterface.hpp>

namespace LifeV {

//! BlockMatrixRN - class for handling a 2-blocks matrix with Robin-Neumann coupling
/*!
    @author Paolo Crosetto

    This class derives both from BlockInterface, which is the base class for the block operators, and from RobinInterface,
    which is a class holding some general methods and attributes for the robin coupling.

    NOTE: this class has been tested for both the GE and GI time discretizations
    (method = monolithicGE and method = monolithicGI).
    NOTE: as the same block matrices are shared between the system matrix and the preconditioner, the preconditioner
    choices available are in principle automatically adapted to the RN case. The preconditioners tested for this case
    are the modular composedDN and the algebraic additive Schwarz AdditiveSchwarz.
 */
class BlockMatrixRN : public BlockMatrix, RobinInterface
{
public:

    //! @name Public Types
    //@{
    typedef BlockMatrix super;
    typedef RobinInterface  superRobin;
    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    BlockMatrixRN(UInt flag):
        super(flag),
        superRobin()
    {}

    //! Destructor
    ~BlockMatrixRN(){}
    //@}

    //! @name Public methods
    //@{

    //! Adds to the r.h.s. the part due to the Robin Coupling and runs BlockMatrix::GlobalAssemble()
    void GlobalAssemble();

    //! Sums all the blocks and the couplings into the system matrix, adds the robin coupling part
    void blockAssembling();


    //! Computes the coupling
    /*!
      computes all the coupling blocks specific for the chosen preconditioner. The coupling is handled
      through an augmented formulation, introducing new variables (multipliers). Needs as input: the global map of the problem,
      the FESpaces of the subproblems to be coupled with their offsets, a std::map holding the two numerations of the
      interface between the two subproblems (the numeration can be different but the nodes must be matching in each
      subdomain), an EpetraVector defined on the multipliers map containing the corresponding dof at the interface (NB: the multipliers map should be constructed from the second numeration in the std::map).
      Note that the FESpaces and the offsets have to be set before calling this method.
      @param map the map of the global problem
      @param locDofMap std::map with the correspondence between the interface dofs for the two different maps in
      the subproblems
      @param numerationInterface vector containing the correspondence of the Lagrange multipliers with the interface dofs
     */
    void coupler(map_shared_ptrtype map,
                 const std::map<ID, ID>& locDofMap,
                 const vector_ptrtype numerationInterface,
                 const Real& timeStep);

    //! sets the data relative to Robin (e.g. the coefficients \f$\alpha_f\f$ and \f$\alpha_s\f$).
    void setDataFromGetPot( const GetPot& data, const std::string& section );

    //! sets the matrix where the Robin contribution will be assembled (which have to passed from outside) and the
    /*! right hand side vector of the linear system, which will be updated with the Robin part.
     */
    void setRobin( matrix_ptrtype& matrix, vector_ptrtype& vec ){setRobinMatrix( matrix ); setRobinRhs( vec );}

    //! sets the matrix where the Robin contribution will be assembled
    void setRobin( matrix_ptrtype& matrix ){setRobinMatrix( matrix ); }
    //@}
private:

};

} // Namespace LifeV

#endif /* BLOCKMATRIXRN_H */
