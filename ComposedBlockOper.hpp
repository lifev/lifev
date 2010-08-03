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

namespace LifeV {

//! ComposedBlockOper - Short description of the class
/*!
    @author Paolo Crosetto
    see [Crosetto, Deparis, Fourestey, Quarteroni 2009]

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
    typedef BlockInterface                           super;
    typedef super::fespace_shared_ptrtype            fespace_ptrtype;
    typedef ComposedPreconditioner<Epetra_Operator>  operator_type;

    ComposedBlockOper():
        super(),
        M_coupling()
    {}

    ~ComposedBlockOper(){}

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

//     //! sets up a vector of raw pointers to the EpetraMaps of each block
//     /*!
//       The number of maps in this vector is not set a-priori. It will be specified in the children classes.
//      */
//     virtual void setBlockMaps(const UInt multipliers, ...)=0;
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

    //! runs GlobalAssemble on the blocks
    /*!
      closes and distributes all the blocks in the vector M_blocks, before computing the preconditioner.
     */
    void GlobalAssemble();

    //! returns true if the operator is set
    /*!
      returns the length of the vector M_blocks
    */
    virtual bool set(){return (bool) M_blocks.size();}

    //!sums the coupling matrices with the corresponding blocks
    /*!
      Everything (but the boundary conditions assembling) must have been set before calling this
    */
    virtual void blockAssembling();

protected:

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

    //! vector of flags saying if the matrix is to be recomputed every time
    std::vector<bool>                                           M_recompute;
    //! vector of coupling matrices
    std::vector<matrix_ptrtype>                                 M_coupling;

private:
    //    static bool                                                 reg;
};

} // Namespace LifeV

#endif /* COMPOSEDBLOCKOPER_H */
