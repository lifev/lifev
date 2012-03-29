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
    @brief File containing a calss for Neumann-Neumann preconditioners

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 02 Jul 2010

     File containing a calss for composed preconditioner of the following type:
    given the matrix \f$A=A_1A_2+A_3A_4f$, \f$A^{-1}\approxA_1^{-1}A_2^{-1}+A_3^{-1}A_4^{-1}\approx P_1^{-1}P_2^{-1}+P_3^{-1}P_4^{-1}\f$ then we compute the preconditioner
    \f$P^{-1}=_4^{-1}P_3^{-1}+P_2^{-1}P_1^{-1}\f$.
 */

#ifndef COMPOSEDDNND_H
#define COMPOSEDDNND_H 1

#include <lifev/core/LifeV.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedNN.hpp>

namespace LifeV
{

//! ComposedPrecDNND - Short description of the class
/*!
    @author Paolo Crosetto

    Class implementing a composed preconditioner of the following type:
    given the matrix \f$A= A_1 A_2+ A_3 A_4\approx  P_1 P_2+ P_3 P_4\f$ then we compute the preconditioner
    \f$P^{-1}=P_4^{-1}P_3^{-1}+P_2^{-1}P_1^{-1}\f$.
    In particular in this case we use for \f$A_1\f$ and \f$A_4\f$ Dirichlet problems, for \f$A_2\f$ and \f$A_3\f$ Neumann problems. The form of the decomposition is the following:
\f$
A=A_{(1)}+A_{(2)}
=
\left(
\begin{array}{ccc}
 2C& 0 & \mathcal I\\
 0&2N& \mathcal I\\
0&0&I
\end{array}
\right)
+
\left(
\begin{array}{ccc}
2C& 0 & 0\\
0 &2N&0\\
\tilde\Delta_1&\tilde\Delta_2&-I
\end{array}
\right).
\f$
So a block factorization for these two matrices can easily be computed,
in fact
\f$
A_{(1)}
=
\left(
\begin{array}{ccc}
 2C& 0 & \mathcal I\\
0 &I&0\\
2\tilde\Delta_1&0&0
\end{array}
\right)
\left(
\begin{array}{ccc}
I& 0 & 0 \\
0 &2 N&2\mathcal I\\
0&0&I
\end{array}
\right)
\f$
while
\f$
A_{(2)}
=
\left(
\begin{array}{ccc}
 2C& 0 & 0\\
0 &I&0\\
2\tilde\Delta_1&0&I
\end{array}
\right)
\left(
\begin{array}{ccc}
I& 0 & 0 \\
0 &2 N&2\mathcal I\\
0&\tilde\Delta_2&0
\end{array}
\right).
\f$
 */
class MonolithicBlockComposedDNND : public MonolithicBlockComposedNN
{
public:

    typedef MonolithicBlockComposedNN super_Type;

    //! @name Constructor and Destructor
    //@{

    MonolithicBlockComposedDNND( const std::vector<Int>& flag, const std::vector<Int>& order ):
            super_Type( flag, order )
    {
    }

    ~MonolithicBlockComposedDNND()
    {}

    //@}

    //! @name Public Methods
    //@{

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
      @param timestep the timestep chosen
     */
    void coupler(mapPtr_Type&      map,
                 const std::map<ID, ID>& locDofMap,
                 const vectorPtr_Type&    numerationInterface,
                 const Real& timeStep,
                 const Real& coefficient,
                 const Real& rescaleFactor);

    //@}
    //!@name Factory Method
    //@{

    static MonolithicBlock* createComposedDNND()
    {
        const Int couplingsDNND[] = { 8, 4, 2, 8, 1, 2 };
        const Int order[] = { MonolithicBlockComposedNN::fluid1, MonolithicBlockComposedNN::solid1, MonolithicBlockComposedNN::fluid2, MonolithicBlockComposedNN::solid2};
        const std::vector<Int> couplingVectorDNND(couplingsDNND, couplingsDNND+6);
        const std::vector<Int> orderVector(order, order+4);
        return new MonolithicBlockComposedDNND(couplingVectorDNND, orderVector);
    }

    //@}

private:

    //static bool reg;
};

} // Namespace LifeV

#endif /* COMPOSEDDNND_H */
