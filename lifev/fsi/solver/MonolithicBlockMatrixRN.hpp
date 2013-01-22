/* -*- mode: c++ -*- */
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
    @file MonolithicBlockMatrixRN.hpp
    @brief file containing a class for handling a 2-blocks matrix with Robin-Neumann coupling

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 09 Jul 2010

 */

#ifndef BLOCKMATRIXRN_H
#define BLOCKMATRIXRN_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/fsi/solver/MonolithicBlockMatrix.hpp>
#include <lifev/fsi/solver/MonolithicRobinInterface.hpp>

namespace LifeV
{

//! MonolithicBlockMatrixRN - class for handling a 2-blocks matrix with Robin-Neumann coupling
/*!
    @author Paolo Crosetto

    This class derives both from MonolithicBlock, which is the base class for the block operators, and from MonolithicRobinInterface,
    which is a class holding some general methods and attributes for the robin coupling.

    NOTE: this class has been tested for both the GE and GI time discretizations
    (method = monolithicGE and method = monolithicGI).
    NOTE: as the same block matrices are shared between the system matrix and the preconditioner, the preconditioner
    choices available are in principle automatically adapted to the RN case. The preconditioners tested for this case
    are the modular composedDN and the algebraic additive Schwarz AdditiveSchwarz.
 */
class MonolithicBlockMatrixRN : public MonolithicBlockMatrix, private MonolithicRobinInterface
{
public:

    //! @name Public Types
    //@{
    typedef MonolithicBlockMatrix super_Type;
    typedef MonolithicRobinInterface  superRobin;
    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MonolithicBlockMatrixRN(const std::vector<Int>& flags /*UInt flag*/):
            super_Type(flags),
            superRobin()
    {}

    //! Destructor
    ~MonolithicBlockMatrixRN() {}
    //@}

    //! @name Public methods
    //@{


    //! sets the data relative to Robin (e.g. the coefficients \f$\alpha_f\f$ and \f$\alpha_s\f$).
    void setDataFromGetPot( const GetPot& data, const std::string& section );

    //! Adds to the r.h.s. the part due to the Robin Coupling and runs MonolithicBlockMatrix::GlobalAssemble()
    void GlobalAssemble();


    //! Computes the specific coupling for a block
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
      @param timeStep the time step
      @param coefficient coefficient, usually is the term multiplying the mass in the time discretization
      @param couplingFlag integer parameter identifying which block is coupled with which. See the method 'couplingMatrix' in @see MonolithicBlock class for a more detailed explanation.
     */
    void coupler(mapPtr_Type& map,
                 const std::map<ID, ID>& locDofMap,
                 const vectorPtr_Type& numerationInterface,
                 const Real& timeStep,
                 const Real& coefficient,
                 const Real& rescaleFactor,
                 UInt        couplingFlag);

    //!Computes the coupling
    void coupler(mapPtr_Type& map,
                 const std::map<ID, ID>& locDofMap,
                 const vectorPtr_Type& numerationInterface,
                 const Real& timeStep,
                 const Real& coefficient,
                 const Real& rescaleFactor);

    //! Sums all the blocks and the couplings into the system matrix, adds the robin coupling part
    void blockAssembling();

    //! sets the matrix where the Robin contribution will be assembled (which have to passed from outside) and the
    /*! right hand side vector of the linear system, which will be updated with the Robin part.
     */
    void setRobin( matrixPtr_Type& matrix, vectorPtr_Type& vec ){setRobinMatrix( matrix ); setRobinRhs( vec );}

    //! sets the matrix where the Robin contribution will be assembled
    void setRobin( matrixPtr_Type& matrix ){setRobinMatrix( matrix ); }
    //@}


    //!@name Factory Method
    //@{
    static MonolithicBlockMatrix*    createAdditiveSchwarzRN()
    {
        const Int couplings[] = { 15, 0, 16 };//to modify (15 to 7) to neglect the coupling (and solve Navier--Stokes)
        const std::vector<Int> couplingVector(couplings, couplings+3);
        return new MonolithicBlockMatrixRN(couplingVector);
    }
    //@}

private:

};

} // Namespace LifeV

#endif /* BLOCKMATRIXRN_H */
