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
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 29 Jun 2010

    A more detailed description of the file (if necessary)
 */

#ifndef COMPOSEDNN_H
#define COMPOSEDNN_H 1

#include <lifemc/lifesolver/ComposedBlockOper.hpp>

namespace LifeV {

//! ComposedNN - Short description of the class
/*!
    @author Paolo Crosetto

    Class implementing a Neumann-Neumann composed preconditioner of the following type:
    given the matrix \f$A=A_1A_2+A_3A_4\approx P_1P_2+P_3P_4\f$ then we compute the preconditioner
    \f$P^1=P_4^{-1}P_3^{-1}+P_2^{-1}P_1^{-1}\f$.
    In particular in this case we use for \f$A_1\f$ and \f$A_2\f$ Dirichlet problems, for \f$A_3\f$ and \f$A_4\f$ Neumann problems.

    Notice that if \f$P^{-1}=(2A)^{-1}+(2A)^{-1}=A^{-1}\f$. Thus the factors that we push_back in the preconditioners should be as close as possible to \f$2A\f$
 */
class ComposedNN: public ComposedBlockOper
{
public:

    typedef ComposedBlockOper          super;
    typedef  ComposedPreconditioner<Ifpack_Preconditioner> composed_prec;


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    ComposedNN(const std::vector<Int>& flag):
        super( flag ),
        M_blockPrecs(),
        M_prec()
    {}


    //! Destructor
    ~ComposedNN(){}

    //@}



    //! @name Virtual methods
    //@{

    //! Solves the preconditioned linear system (used only when dealing with a preconditioner)
    /*!
      Provided the linear solver and the right hand side this method computes the preconditioners, builds the composed
      operator (of type ComposedPreconditioner) and solves the preconditioned linear system.
        @param rhs right hand side of the linear system
        @param result output result
        @param linearSolver the linear system
     */
    virtual int   solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver);



    //!Applies the correspondent boundary conditions to every block
    /*!
      note that this method must be called after blockAssembling(), that sums the coupling conditions to the blocks. For
      this type of preconditioners this method is overloaded. In fact in the diagonalization for the essential boundary
      conditions the value replaced on the diagonal must be 2 instead of 1.
      \param time: time
     */
    void applyBoundaryConditions(const Real& time, const UInt i);

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
    virtual void coupler( map_shared_ptrtype& map,
                          const std::map<ID, ID>& locDofMap,
                          const vector_ptrtype& numerationInterface,
                          const Real& timeStep);
    //@}

    //! @name Methods
    //@{

    //! Sets the parameters needed by the preconditioner from data file (creates the Ifpack list)
    /*!
        @param data GetPot object reading the text data file
        @param section string specifying the path in the data file where to find the options for the operator
     */
    void setDataFromGetPot( const GetPot& data, const std::string& section );

    //! Multiplies the block times 2 and calls super::push_back_matrix(...)
    /*!
      \param Mat: block matrix
      \param recompute: flag stating if the matrix need to be recomputed
     */
    void push_back_matrix( const  matrix_ptrtype& Mat, const  bool recompute );

    /*! Multiplies the block times 2 and calls super::replace_matrix(...) in the position "position" specified in
      input and in the shifted position "position"+2
      \param oper: input matrix
      \param position: position
     */
    void replace_matrix( const matrix_ptrtype& oper, UInt position );

    bool set(){return (bool) M_blockPrecs.get() && M_blockPrecs->getNumber();}
    //@}

protected:

    boost::shared_ptr<ComposedPreconditioner<ComposedPreconditioner<Ifpack_Preconditioner> > >            M_blockPrecs;
    Teuchos::ParameterList                                                 M_list;
    std::vector<boost::shared_ptr<Ifpack_Preconditioner> >                 M_prec;

private:

    boost::shared_ptr< composed_prec > M_firstCompPrec ;
    boost::shared_ptr< composed_prec > M_secondCompPrec;

};

} // Namespace LifeV

#endif /* COMPOSEDNN_H */
