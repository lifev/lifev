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
  \include ../../testsuite/test_monolithic/fluidstructure.dox
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

    NOTE: this class is also the base class for other types of preconditioners, like ComposedDN2, ComposedDND. In fact for
    instance it is used as F-S block in the preconditioners for the GI matrix in MonolithicGI
 */

#ifndef COMPOSEDDN_H
#define COMPOSEDDN_H 1

#include <life/lifecore/life.hpp>

#include <lifemc/lifesolver/ComposedBlockOper.hpp>

#include <lifemc/lifealg/IfpackComposedPrec.hpp>

namespace LifeV {

//! ComposedDN - Short description of the class
/*!
    @author Paolo Crosetto
    @see \ref CDFQ


 */
class ComposedDN : public ComposedBlockOper
{
public:
    typedef ComposedBlockOper super;

    ComposedDN():
        super(),
        M_couplingFlag(7),
        M_blockPrecs(),
        M_uMap(),
        M_pMap(),
        M_dMap(),
        M_interfaceMap(),
        M_multipliers(0)
    {
        //M_bch.resize(2);
    }

    ComposedDN( Int flag, Int superFlag = 16 ):
        super( superFlag ),
        M_couplingFlag( flag ),
        M_blockPrecs(),
        M_uMap(),
        M_pMap(),
        M_dMap(),
        M_interfaceMap(),
        M_multipliers(0)
    {
    }

    //! @name public methods
    //@{
    //! Solves the preconditioned linear system
    /*!
      Provided the linear solver and the right hand side this method computes the solution and returns it into
      the result vector.
        @param rhs right hand side of the linear system
        @param result output result
        @param linearSolver the linear system
     */
    int     solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver);

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
    virtual void coupler(map_shared_ptrtype map,
                         const std::map<ID, ID>& locDofMap,
                         const vector_ptrtype numerationInterface,
                         const Real& timeStep);

    //!pushes back the preconditioner for a block
    /*!
      In this case M_blockPrecs is of type IfpackComposedPrec, thus this method calls IfpackComposedPrec::push_back(...)
      which computes the AAS preconditioner for the input matrix
      \param Mat: input matrix
     */
    virtual void    push_back_precs( matrix_ptrtype& Mat);

    //! sets the parameters related to M_blockPrecs from the data file
    /*!
      \param dataFile: GetPot data file
      \param section: string identifying the section in the data file
     */
    void setDataFromGetPot( const GetPot& dataFile,
                            const std::string& section );

    bool set(){return (bool) M_blockPrecs.get() && M_blockPrecs->getNumber();}

    void setComm( boost::shared_ptr<Epetra_Comm> comm )
    {
        M_comm = comm;
        M_blockPrecs.reset( new IfpackComposedPrec(M_comm));
    }

protected:

    virtual void    replace_precs( matrix_ptrtype& Mat, UInt position);
    boost::shared_ptr<IfpackComposedPrec>            M_blockPrecs;
    const UInt                                       M_couplingFlag;
    map_ptrtype                                      M_uMap;
    map_ptrtype                                      M_pMap;
    map_ptrtype                                      M_dMap;
    map_ptrtype                                      M_interfaceMap;
    UInt                                             M_multipliers;
private:
    //    static bool                                      reg;

};

} // Namespace LifeV

#endif /* COMPOSEDDN_H */
