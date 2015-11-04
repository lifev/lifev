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
    @file
    @brief This file contains a class implementing a preconditioner for the Fluid-Structure Interaction system
    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 08 Jun 2010

    We call the monolithic GCE matrix (with the fluid snd structure blocks C and N, couplings B and D):
    \f$
    A=\left(\begin{array}{cc}
    C&B\\
    D&N
    \end{array}\right)\f$.

    The preconditioner is obtained from a block Gauss-Seidel preconditioner
    \f$P_{DN}=
    \left(\begin{array}{cc}
    C&B\\
    0&N
    \end{array}\right)\f$,
    decomposing it into two factors
    \f$P_{DN}=P_1P_2
    \left(
    \begin{array}{cc}
    I&0\\
    0&N
    \end{array}\right)
    \left(\begin{array}{cc}
    C&B\\
    0&I
    \end{array}\right)\f$ and applying a preconditioning strategy (algebraic additive Schwarz \f$P_{AS}\f$)
    to each factor, so that \f$ P^{-1}=(P_{AS}(P_2))^{-1}(P_{AS}(P_1))^{-1}\f$.

    NOTE: this class is used also in the geometry implicit case, with an additional factor for the mesh motion, which is
    coupled to the solid block using the default coupling method MonolithicBlockComposed::coupler. In that case the
    preconditioner is decomposed in three factors, the fluid and structure ones being the same as for the GCE case.
    NOTE2: this class is also the base class for other types of preconditioners, like ComposedDN2, ComposedDND. In fact for
    instance it is used as F-S block in the preconditioners for the GI matrix in FSIFSIMonolithicGI
 */

#ifndef COMPOSEDDN_H
#define COMPOSEDDN_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/algorithm/PreconditionerComposed.hpp>

#include <lifev/fsi/solver/MonolithicBlockComposed.hpp>

namespace LifeV
{

//! MonolithicBlockComposedDN - Short description of the class
/*!
    @author Paolo Crosetto
    @see \cite CrosettoEtAl2009


 */
class MonolithicBlockComposedDN : public MonolithicBlockComposed
{
public:
    typedef MonolithicBlockComposed super_Type;

    MonolithicBlockComposedDN ( const std::vector<Int>& flag, const std::vector<Int>& order) :
        super_Type ( flag, order ),
        M_blockPrecs()
    {
    }

    //! @name public methods
    //@{

    //! sets the parameters related to M_blockPrecs from the data file
    /*!
      \param dataFile: GetPot data file
      \param section: string identifying the section in the data file
     */
    void setDataFromGetPot ( const GetPot& dataFile,
                             const std::string& section );

    //! Solves the preconditioned linear system
    /*!
      Provided the linear solver and the right hand side this method computes the solution and returns it into
      the result vector.
        @param rhs right hand side of the linear system
        @param result output result
        @param linearSolver the linear system
     */
    int     solveSystem ( const vector_Type& rhs, vector_Type& step, solverPtr_Type& linearSolver);

    //! Computes the coupling
    /*!
      computes all the coupling blocks specific for the chosen preconditioner. The coupling is handled
      through an augmented formulation, introducing new variables (multipliers). Needs as input: the global map of the problem,
      the FESpaces of the subproblems to be coupled with their offsets, a std::map holding the two numerations of the
      interface between the two subproblems (the numeration can be different but the nodes must be matching in each
      subdomain), an EpetraVector defined on the multipliers map containing the corresponding dof at the interface (NB: the multipliers map should be constructed from the second numeration in the std::map).

      In this case the coupling matrices are two:
      \f$
      C_1=
      \left(
      \begin{array}{cc}
      I&0\\
      0&0
      \end{array}
      \right)\f$
      and
      \f$
      C_2=
      \left(
      \begin{array}{cc}
      0&C\\
      0&I
      \end{array}
      \right)\f$

      Note that the FESpaces and the offsets have to be set before calling this method.
      @param map the map of the global problem
      @param locDofMap std::map with the correspondence between the interface dofs for the two different maps in
      the subproblems
      @param numerationInterface vector containing the correspondence of the Lagrange multipliers with the interface dofs
     */
    virtual void coupler (mapPtr_Type& map,
                          const std::map<ID, ID>& locDofMap,
                          const vectorPtr_Type& numerationInterface,
                          const Real& timeStep,
                          const Real& coefficient,
                          const Real& rescaleFactor);

    //!pushes back the preconditioner for a block
    /*!
      In this case M_blockPrecs is of type PreconditionerComposed, thus this method calls PreconditionerComposed::push_back(...)
      which computes the AAS preconditioner for the input matrix
      \param Mat: input matrix
     */
    virtual void    push_back_precs ( matrixPtr_Type& Mat);

    //! returns the true if the preconditioner has at leas one factor computed
    bool set()
    {
        return  M_blockPrecs->preconditionerCreated();
    }

    /*! copies the shared_ptr to the communicator in the member M_comm and builds the empty ComposedPreconditioneronditioner
    M_blockPrecs
    */
    void setComm ( std::shared_ptr<Epetra_Comm> comm )
    {
        M_comm = comm;
        M_blockPrecs.reset ( new PreconditionerComposed (M_comm) );
    }

    const std::vector<std::shared_ptr<Preconditioner> >& blockPrecs() const
    {
        return M_blockPrecs->composedPreconditionerPtr()->Operator();
    }

    //@}
    //!@name Factory Methods
    //@{

    static MonolithicBlock* createComposedDN()
    {
        const Int order[] = { MonolithicBlockComposed::solid, MonolithicBlockComposed::fluid};
        const Int couplingsDN[] = { 0, 7};
        const std::vector<Int> couplingVectorDN (couplingsDN, couplingsDN + 2);
        const std::vector<Int> orderVector (order, order + 2);
        return new MonolithicBlockComposedDN (couplingVectorDN, orderVector);
    }


    static MonolithicBlock* createComposedDN2()
    {
        const Int order[] = { MonolithicBlockComposed::fluid, MonolithicBlockComposed::solid};
        const Int couplingsDN2[] = { 8, 6};
        const std::vector<Int> couplingVectorDN2 (couplingsDN2, couplingsDN2 + 2);
        const std::vector<Int> orderVector (order, order + 2);
        return new MonolithicBlockComposedDN (couplingVectorDN2, orderVector);
    }

    static MonolithicBlock* createComposedDNGI()
    {
        //! Factorization in three:
        /*!
        - Solid: Neumann
        | I | 0 | 0 | 0 |
        |---+---+---+---|
        | 0 | S | 0 | 0 |
        |---+---+---+---|
        | 0 | 0 | I | 0 |
        |---+---+---+---|
        | 0 | 0 | 0 | I |
        - ALE: Dirichlet
        | I |  0 | 0 | 0 |
        |---+----+---+---|
        | 0 |  I | 0 | 0 |
        |---+----+---+---|
        | 0 |  0 | I | 0 |
        |---+----+---+---|
        | 0 | -C | 0 | H |
        - Fluid: Dirichlet
        | F | 0     | C | SD |
        |---+-------+---+----|
        | 0 | S     | 0 | 0  |
        |---+-------+---+----|
        | C | -C/dt | 0 | 0  |
        |---+-------+---+----|
        | 0 | 0     | 0 | I  |
             */
        const Int order[] = { MonolithicBlockComposed::solid, MonolithicBlockComposed::mesh, MonolithicBlockComposed::fluid };
        const Int couplingsDNGI[] = { 0/*solid*/, 7/*fluid*/, 16/*ALE*/ };
        const std::vector<Int> couplingVectorDNGI (couplingsDNGI, couplingsDNGI + 3);
        const std::vector<Int> orderVector (order, order + 3);
        return new MonolithicBlockComposedDN ( couplingVectorDNGI, orderVector );
    }


    static MonolithicBlock* createComposedDN2GI()
    {
        //! Factorization in three:
        /*!
        - Fluid: Dirichlet
               | F | 0 | 0 | SD |
               |---+---+---+----|
               | 0 | I | 0 |  0 |
               |---+---+---+----|
               | C | 0 | I |  0 |
               |---+---+---+----|
               | 0 | 0 | 0 |  I |
        - Solid: Neumann
               | I |     0 | 0 | 0 |
               |---+-------+---+---|
               | 0 |     S | C | 0 |
               |---+-------+---+---|
               | 0 | -C/dt | 0 | 0 |
               |---+-------+---+---|
               | 0 |     0 | 0 | I |
        - ALE: Dirichlet
               | I |  0 | 0 | 0 |
               |---+----+---+---|
               | 0 |  I | 0 | 0 |
               |---+----+---+---|
               | 0 |  0 | I | 0 |
               |---+----+---+---|
               | 0 | -C | 0 | H |
             */
        const Int order[] = {   MonolithicBlockComposed::fluid, MonolithicBlockComposed::solid, MonolithicBlockComposed::mesh };
        const Int couplingsDN2GI[] = { 8, 6, 16 };
        const std::vector<Int> couplingVectorDN2GI (couplingsDN2GI, couplingsDN2GI + 3);
        const std::vector<Int> orderVector (order, order + 3);
        return new MonolithicBlockComposedDN ( couplingVectorDN2GI, orderVector );
    }

    //@}

protected:

    //! @name Protected Methods
    //@{

    /*!
      Replaces the preconditioner in M_blockPrecs with another one that is already constructed
    */
    virtual void    replace_precs ( matrixPtr_Type& Mat, UInt position);
    //@}

    //! @name Protected Members
    //@{

    /*!
      Pointer to an PreconditionerComposed object containing the preconditioners for each block
    */
    std::shared_ptr<PreconditionerComposed>        M_blockPrecs;
    //@}

private:
    //    static bool                                      reg;

};

} // Namespace LifeV

#endif /* COMPOSEDDN_H */
