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
    @brief this file contains a class which is suited for handling a block-structured matrix that can be written as a
    multiplication of a variable number of factors. It contains a vector of pointers for each factor, BCHandler, FESpace
    and for each coupling part.

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 08 Jun 2010

 */

#ifndef MONOLITHICBLOCKCOMPOSED_H
#define MONOLITHICBLOCKCOMPOSED_H 1

#include <lifev/core/LifeV.hpp>
#include <Epetra_Operator.h>
#include <lifev/core/algorithm/ComposedOperator.hpp>
#include <lifev/fsi/solver/MonolithicBlock.hpp>

#include <boost/scoped_ptr.hpp>

namespace LifeV
{

//! MonolithicBlockComposed - Class handling block-structured preconditioners
/*!
    @author Paolo Crosetto
    see \cite CrosettoEtAl2009

    Pure virtual class which is suited for handling a block-structured matrix that can be written as a
    multiplication of a variable number of factors. It contains a vector of pointers for each factor, BCHandler, FESpace
    and for each coupling part.

    It can be used e.g. to build preconditioners for the FSI system which are composed of different factors.
    It derives from the pure virtual class MonolithicBlock, it is still pure virtual and it is intended to be the base
    class for every preconditioner that can be split into factors. Notice that the vector of the preconditioned factors
    is missing, because it has to be implemented in the children. It could be for instance of type
    ComposedOperator<T> or PreconditionerComposed.
 */
class MonolithicBlockComposed            : public MonolithicBlock
{
public:

    enum Block { solid, fluid, mesh };

    //! @name Public Types
    //@{
    typedef MonolithicBlock                           super_Type;
    typedef super_Type::fespacePtr_Type            fespacePtr_Type;
    typedef ComposedOperator<Epetra_Operator>  operatorPtr_Type;
    //@}

    //! @name Constructor and Destructor
    //@{
    //! Constructor
    /**
       The coupling and the order have to be specified in input. In particular the order specifies which block go
       first, second etc.
       \param flags: vector of flags specifying the type of coupling between the different blocks that we chose for this operator
       \param order: vector specifying the order of the blocks.
     */
    MonolithicBlockComposed (const std::vector<Int>& flags, const std::vector<Int>& order) :
        super_Type(),
        M_recompute (order.size() ),
        M_coupling(),
        M_couplingFlags (new std::vector<Int> (flags) ), // here I copy, so that the input param can be destroyed
        M_blockReordering (new std::vector<Int> (order) )
    {}


    ~MonolithicBlockComposed() {}
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
    virtual int   solveSystem ( const vector_Type& rhs, vector_Type& step, solverPtr_Type& linearSolver) = 0;

    //! Sets the parameters needed by the preconditioner from data file
    /*!
        @param data GetPot object reading the text data file
        @param section string specifying the path in the data file where to find the options for the operator
     */
    virtual void  setDataFromGetPot (const GetPot& data, const std::string& section) = 0;

    //! returns true if the operator is set
    /*!
      returns the length of the vector M_blocks
    */
    virtual bool set() = 0;

    //@}
    //!@name Public Methods
    //@{


    //! runs GlobalAssemble on the blocks
    /*!
      closes and distributes all the blocks in the vector M_blocks, before computing the preconditioner.
     */
    void GlobalAssemble();

    //!sums the coupling matrices with the corresponding blocks
    /*!
      Everything (but the boundary conditions assembling) must have been set before calling this
    */
    virtual void blockAssembling();


    //! adds a default coupling matrix for a specified block.
    /*!
      The default coupling matrix it is an identity matrix with zeros on the diagonal corresponding to the specified block.
      The couplings are specified through the vector M_couplingFlags, as usual. If the coupling block is not the last one
      in the M_coupling vector then it is inserted at the specified position.
      @param map: global MapEpetra of the problem
      @param locDofMap: std::map holding the connections between the coupling interface dofs
      @param numerationInterface: the numeration of the interface dofs
      @param timeStep: time step
      @param couplingBlock: UInt specifying the position of the coupling block to be added. Not used in this case, since the coupling flag for each block
      is contained in the vector M_couplingFlags. See MonolithicBlock::couplingMatrix to understand what the values
      for this flag correspond to.
      @param couplingBlock: flag specifying which block is considered (must not exceed the size of the vector of blocks,
      otherwise a std::bad_alloc exception is thrown). The coupling for each block is specified by a static vector passed to the constructor.
     */
    void coupler (mapPtr_Type& map,
                  const std::map<ID, ID>& locDofMap,
                  const vectorPtr_Type& numerationInterface,
                  const Real& timeStep,
                  const Real& coefficient,
                  const Real& rescaleFactor,
                  UInt couplingBlock
                 );


    //! pushes a block at the end of the vector
    /*!
      adds a new block
        @param Mat: block matrix to push
        @param recompute: flag stating wether the preconditioner for this block have to be recomputed at every time step
     */
    virtual void    push_back_matrix (const matrixPtr_Type& Mat, const bool recompute);

    //! Merges an input MonolithicBlockComposed operator with this one
    /*!
       pushes the operator vector of the input operator at the end of the operatoe vector of this instance, and does the same for the
       vector of coupling matrices.
      @param Oper: input operator
     */
    virtual void push_back_oper ( MonolithicBlockComposed& Oper);


    //! Pushes an extra coupling matrix at the end of the vector of coupling matrices
    /*!
      @param coupling: extra coupling matrix
     */
    virtual void push_back_coupling ( matrixPtr_Type& coupling);

    //! replaces a block
    /*!
      replaces a block on a specified position in the vector
        @param Mat: block matrix to push
        @param index: position in the vector
     */
    virtual void replace_matrix ( const matrixPtr_Type& oper, UInt position ); //{M_blocks.replace(oper, position);}


    //! replaces a coupling block
    /*!
      replaces a coupling block on a specified position in the vector
        @param Mat block matrix to push
        @param index position in the vector
     */
    virtual void replace_coupling ( const matrixPtr_Type& Mat, UInt index);

    //! pushes a block at the end of the vector
    /*!
      adds a new block
        @param Mat block matrix to push
        @param recompute flag stating wether the preconditioner for this block have to be recomputed at every time step
     */
    virtual void addToCoupling ( const matrixPtr_Type& Mat, UInt position);

    //!
    /*!
      adds an entry to the coupling matrix
        @param entry entry
        @param row row for the insertion
        @param col colon for the insertion
     */
    void addToCoupling ( const Real& entry , UInt row, UInt col, UInt position );

    //@}

    //!@name Get Methods
    //@{
    //! returns the vector of flags (by const reference).
    const std::vector<bool>& recompute()
    {
        return M_recompute;
    }

    //! returns the vector of pointers to the coupling blocks (by const reference).
    const std::vector<matrixPtr_Type>& couplingVector() const
    {
        return M_coupling;
    }

    //@}

    //!@name Set Methods
    //@{

    //! turns on/off the recomputation of the preconditioner for a specified factor
    void setRecompute ( UInt position, bool flag )
    {
        M_recompute[position] = flag;
    }

    const UInt whereIsBlock ( UInt position ) const;

    //@}
protected:

    //! @name Protected Methods
    //@{

    //!sums the coupling matrix in the specified position with the corresponding block
    /*!
      Everything (but the boundary conditions assembling) must have been set before calling this
    */
    void blockAssembling (const UInt k);


    //! swaps the blocks
    /*!
      swaps two elements in the vectors of blocks, couplings, offsets, BC handlers, FE spaces, flags M_recompute
      \param i: element to swap
      \param j: element to swap
     */
    virtual void swap (const UInt i, const UInt j);
    //@}

    //! @name Protected Members
    //@{

    //! vector of flags saying if the matrix is to be recomputed every time
    std::vector<bool>                                           M_recompute;
    //! vector of coupling matrices
    std::vector<matrixPtr_Type>                                 M_coupling;

    //! vector of flags specifying the coupling strategy for each block.
    /*!
      In particular for each block the method couplingMatrix is called with the corresponding flag in this vector.
      The coupling vector is passed to the constructor when the preconditioner is registered in the factory. So
      each coupling vector defines a new preconditioner type. It is created statically before the registration, then it
      is copied into this scoped_ptr.
     */
    std::unique_ptr<std::vector<Int> >                        M_couplingFlags;

    //! vector of reordering for the different blocks.
    /*!the order in which the factors are allpied is specified by this
    vector. e.g. the fisrt block to be applied corresponds to the number M_blockReordering[0] in the vector
    M_blocks of blocks. This vector is assigned in the coupler method of each class.
    */
    std::unique_ptr<std::vector<Int> >                      M_blockReordering;

    //@}

private:

};

} // Namespace LifeV

#endif /* MONOLITHICBLOCKCOMPOSED_H */
