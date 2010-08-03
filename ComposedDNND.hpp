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
    @file
    @brief File containing a calss for Neumann-Neumann preconditioners

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 02 Jul 2010

     File containing a calss for composed preconditioner of the following type:
    given the matrix \f$A=2A_12A_2+2A_32A_4\approx 2P_12P_2+2P_32P_4\f$ then we compute the preconditioner
    \f$P^{-1}=\frac12P_4^{-1}\frac12P_3^{-1}+\frac12P_2^{-1}\frac12P_1^{-1}\f$.
 */

#ifndef COMPOSEDDNND_H
#define COMPOSEDDNND_H 1

#include <life/lifecore/life.hpp>
#include <lifemc/lifesolver/ComposedNN.hpp>

namespace LifeV {

//! ComposedPrecDNND - Short description of the class
/*!
    @author Paolo Crosetto

    Class implementing a composed preconditioner of the following type:
    given the matrix \f$A= \frac12A_1 \frac12A_2+ \frac12A_3 \frac12A_4\approx  \frac12P_1 \frac12P_2+ \frac12P_3 \frac12P_4\f$ then we compute the preconditioner
    \f$P^1=2P_4^{-1}2P_3^{-1}+2P_2^{-1}2P_1^{-1}\f$.
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
0 &\frac12N&0\\
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
 \frac12C& 0 & \mathcal I\\
0 &I&0\\
\frac12\tilde\Delta_1&0&0
\end{array}
\right)
\left(
\begin{array}{ccc}
I& 0 & 0 \\
0 &\frac12 N&\frac12\mathcal I\\
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
 \frac12C& 0 & 0\\
0 &I&0\\
\frac12\tilde\Delta_1&0&I
\end{array}
\right)
\left(
\begin{array}{ccc}
I& 0 & 0 \\
0 &\frac12 N&\frac12\mathcal I\\
0&\tilde\Delta_2&0
\end{array}
\right).
\f$
 */
class ComposedDNND : public ComposedNN
{
public:

    typedef ComposedNN super;

    ComposedDNND():
        super()
    {
    }

    ~ComposedDNND()
    {}


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
    void coupler(map_shared_ptrtype      map,
                         const std::map<ID, ID>& locDofMap,
                         const vector_ptrtype    numerationInterface,
                         const Real& timeStep);


    //! Sets the parameters needed by the preconditioner from data file
    /*!
        @param data GetPot object reading the text data file
        @param section string specifying the path in the data file where to find the options for the operator
     */
    void setDataFromGetPot( const GetPot&      dataFile,
                            const std::string& section ){}

private:

    //static bool reg;
};

} // Namespace LifeV

#endif /* COMPOSEDDNND_H */
