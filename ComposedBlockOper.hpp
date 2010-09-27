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
    @file composedBlockOper.hpp
    @brief this file contains a class which is suited for handling a block-structured matrix that can be written as a
    multiplication of a variable number of factors. It contains a vector of pointers for each factor, BCHandler, FESpace
    and for each coupling part.

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 08 Jun 2010

 */

#ifndef COMPOSEDBLOCKOPER_H
#define COMPOSEDBLOCKOPER_H 1

#include <life/lifecore/life.hpp>
#include <lifemc/lifesolver/BlockInterface.hpp>

#include <boost/scoped_ptr.hpp>

namespace LifeV {

//! ComposedBlockOper - Class handling block-structured preconditioners
/*!
    @author Paolo Crosetto
    see \ref CDFQ

    Pure virtual class which is suited for handling a block-structured matrix that can be written as a
    multiplication of a variable number of factors. It contains a vector of pointers for each factor, BCHandler, FESpace
    and for each coupling part.

    It can be used e.g. to build preconditioners for the FSI system which are composed of different factors.
    It derives from the pure virtual class BlockInterface, it is still pure virtual and it is intended to be the base
    class for every preconditioner that can be split into factors. Notice that the vector of the preconditioned factors
    is missing, because it has to be implemented in the children. It could be for instance of type
    ComposedPreconditioner<T> or IfpackComposedPrec.
 */
class ComposedBlockOper            : public BlockInterface
{
public:

    enum Block { solid, fluid, mesh };

    //! @name Public Types
    //@{
    typedef BlockInterface                           super;
    typedef super::fespace_shared_ptrtype            fespace_ptrtype;
    typedef ComposedPreconditioner<Epetra_Operator>  operator_type;
    //@}

    //! @name Constructor and Destructor
    //@{
    ComposedBlockOper(const std::vector<Int>& flags, const std::vector<Block>& order):
        super(),
        M_recompute(),
        M_coupling(),
        M_couplingFlags(new std::vector<Int>(flags)),// here I copy, so that the input param can be destroyed
        M_blockReordering(new std::vector<Block>(order))
    {}


    ~ComposedBlockOper(){}
    //@}

    //! @name Pure virtual methods
    //@{
    //! Solves the preconditioned linear system
    /*!
      Provided the linear solver and the right hand side this method computes the solution and returns it into
      the result vector.
        @param rhs right hand side of the linear system
        @param result output result
        @param linearSolver the linear system
     */
    virtual int   solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver)=0;

    //! Sets the parameters needed by the preconditioner from data file
    /*!
        @param data GetPot object reading the text data file
        @param section string specifying the path in the data file where to find the options for the operator
     */
    virtual void  setDataFromGetPot(const GetPot& data, const std::string& section)=0;

    //@}

    //! pushes a block at the end of the vector
    /*!
      adds a new block
        @param Mat block matrix to push
        @param recompute flag stating wether the preconditioner for this block have to be recomputed at every time step
     */
    virtual void    push_back_matrix(const matrix_ptrtype& Mat, const bool recompute);

    //! replaces a block
    /*!
      replaces a block on a specified position in the vector
        @param Mat block matrix to push
        @param index position in the vector
     */
    virtual void replace_matrix( const matrix_ptrtype& oper, UInt position );//{M_blocks.replace(oper, position);}


    //! replaces a coupling block
    /*!
      replaces a coupling block on a specified position in the vector
        @param Mat block matrix to push
        @param index position in the vector
     */
    virtual void replace_coupling( const matrix_ptrtype& Mat, UInt index);


    //! runs GlobalAssemble on the blocks
    /*!
      closes and distributes all the blocks in the vector M_blocks, before computing the preconditioner.
     */
    void GlobalAssemble();

    //! returns true if the operator is set
    /*!
      returns the length of the vector M_blocks
    */
    virtual bool set()=0;

    //!sums the coupling matrices with the corresponding blocks
    /*!
      Everything (but the boundary conditions assembling) must have been set before calling this
    */
    virtual void blockAssembling();

    //! pushes a block at the end of the vector
    /*!
      adds a new block
        @param Mat block matrix to push
        @param recompute flag stating wether the preconditioner for this block have to be recomputed at every time step
     */
    virtual void addToCoupling( const matrix_ptrtype& Mat, UInt position);

    virtual void push_back_oper( ComposedBlockOper& Oper);
    //@}

    //!@name Getters
    //@{
    //! returns the vector of flags (by const reference).
    const std::vector<bool>& getRecompute(){return M_recompute;}

    //! returns the vector of pointers to the coupling blocks (by const reference).
    const std::vector<matrix_ptrtype> getCouplingVector(){return M_coupling;}

    //! adds a default coupling matrix for a specified block.
    /*!
      The default coupling matrix it is an identity matrix with zeros on the diagonal corresponding to the specified block.
      The couplings are specified through the vector M_couplingFlags, as usual. If the coupling block is not the last one
      in the M_coupling vector then it is inserted at the specified position.
      @param map: global EpetraMap of the problem
      @param locDofMap: std::map holding the connections between the coupling interface dofs
      @param numerationInterface: the numeration of the interface dofs
      @param timeStep: time step
      @param couplingBlock: UInt specifying the position of the coupling block to be added.
     */
    void coupler(map_shared_ptrtype& map,
                 const std::map<ID, ID>& locDofMap,
                 const vector_ptrtype& numerationInterface,
                 const Real& timeStep,
                 UInt couplingBlock
                 );

    virtual void push_back_coupling( matrix_ptrtype& coupling);

    //@}

protected:

    //! @name Protected Methods
    //@{

    //!sums the coupling matrix in the specified position with the corresponding block
    /*!
      Everything (but the boundary conditions assembling) must have been set before calling this
    */
    void blockAssembling(const UInt k);


    //! swaps the blocks
    /*!
      swaps two elements in the vectors of blocks, couplings, offsets, BC handlers, FE spaces, flags M_recompute
      \param i: element to swap
      \param j: element to swap
     */
    virtual void swap(const UInt i, const UInt j);
    //@}

    //! @name Protected Members
    //@{

    //! vector of flags saying if the matrix is to be recomputed every time
    std::vector<bool>                                           M_recompute;
    //! vector of coupling matrices
    std::vector<matrix_ptrtype>                                 M_coupling;

    //! vector of flags specifying the coupling strategy for each block.
    /*!
      In particular for each block the method couplingMatrix is called with the corresponding flag in this vector.
      The coupling vector is passed to the constructor when the preconditioner is registered in the factory. So
      each coupling vector defines a new preconditioner type. It is created statically before the registration, then it
      is copied into this scoped_ptr.
     */
    boost::scoped_ptr<std::vector<Int> >                        M_couplingFlags;

    //! vector of reordering for the different blocks.
    /*!the order in which the factors are allpied is specified by this
    vector. e.g. the fisrt block to be applied corresponds to the number M_blockReordering[0] in the vector
    M_blocks of blocks. This vector is assigned in the coupler method of each class.
    */
    boost::scoped_ptr<std::vector<Block> >                      M_blockReordering;
    //@}

private:


};

} // Namespace LifeV

#endif /* COMPOSEDBLOCKOPER_H */
