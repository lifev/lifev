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
    @brief This file contains a pure virtual class for the linear operators with a block structure
    (i.e. block matrices and preconditioners). The specializations of this class can be used to handle many
    types of composed preconditioners and composed operators.

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 07 Jun 2010
 */

#ifndef BLOCKINTERFACE_H
#define BLOCKINTERFACE_H 1

#include <cstdarg>
#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>

#include <life/lifefem/FESpace.hpp>

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>
#include <lifemc/lifealg/ComposedOperator.hpp>
#include <life/lifefem/bcManage.hpp>

namespace LifeV {

//! BlockInterface - This is a pure virtual class for the linear operators with a block structure
/*!    (i.e. block matrices and preconditioners).
    @author Paolo Crosetto

    The specializations of this class can be used to handle many
    types of composed preconditioners and composed operators. It conteins a vector of pointers to EpetraMatrix
    holding the diferent blocks. These matrices should be zero everywhere but on the diagonal block corresponding
    to the corresponding problem. This class also handles the coupling through a vector of pointer to coupling matrices.
 */

class BlockInterface
{
public:

    //! @name Public Types
    //@{

    typedef EpetraVector                                               vector_type;
    typedef boost::shared_ptr< vector_type >                           vector_ptrtype;
    typedef EpetraMatrix< Real >                                       matrix_type;
    typedef boost::shared_ptr< matrix_type >                           matrix_ptrtype;
    typedef boost::shared_ptr< Epetra_Operator >                       epetra_operator_ptrtype;
    typedef boost::shared_ptr< EpetraPreconditioner >                  epetra_preconditioner_ptrtype;
    typedef matrix_type::matrix_type                                   epetra_matrix_type;
    typedef SolverTrilinos                                             solver_type;
    typedef boost::shared_ptr< SolverTrilinos >                        solver_ptrtype;
    typedef boost::shared_ptr< FESpace<RegionMesh3D<LinearTetra>, EpetraMap> >  fespace_shared_ptrtype;
    typedef fespace_shared_ptrtype                                     fespace_ptrtype;
    //    typedef FESpace<RegionMesh3D<LinearTetra>, EpetraMap>*                 fespace_ptrtype;
    typedef EpetraMap*                                                 map_ptrtype;
    typedef boost::shared_ptr< EpetraMap >                             map_shared_ptrtype;
    //typedef BCHandler*                                                 bchandler_ptrtype;
    typedef boost::shared_ptr< BCHandler >                             bchandler_ptrtype;
    //@}


    //! @name Constructor & Destructor
    //@{
    //! Empty Constructor
    BlockInterface():
        M_bch(),
        M_blocks(),
        M_FESpace(),
        M_comm()
    {}

//     BlockInterface( Int flag ):
//         M_bch(),
//         M_blocks(),
//         M_FESpace(),
//         M_comm(),
//         M_superCouplingFlag(flag)
//     {}

    //! Destructor
    ~BlockInterface()
{
//     free(M_offset);
//     free(M_FESpace);
}
    //@}



    //! @name Pure virtual methods
    //@{
    //! Solves the preconditioned linear system (used only when dealing with a preconditioner)
    /*!
      Provided the linear solver and the right hand side this method computes the solution and returns it into
      the result vector.
        @param rhs right hand side of the linear system
        @param result output result
        @param linearSolver the linear system
     */
    virtual int  solveSystem( const vector_type& rhs, vector_type& result, solver_ptrtype& linearSolver)=0;

    //! Sets the parameters needed by the preconditioner from data file
    /*!
        @param data GetPot object reading the text data file
        @param section string specifying the path in the data file where to find the options for the operator
     */
    virtual void setDataFromGetPot(const GetPot& data, const std::string& section)=0;

    //! pushes a block at the end of the vector
    /*!
      adds a new block
        @param Mat block matrix to push
        @param recompute flag stating wether the preconditioner for this block have to be recomputed at every time step
     */
    virtual void push_back_matrix( const matrix_ptrtype& Mat, const bool recompute ) =0;

    //!
    /*!
      adds a new block
        @param Mat block matrix to push
        @param recompute flag stating wether the preconditioner for this block have to be recomputed at every time step
     */
    virtual void addToCoupling( const matrix_ptrtype& Mat, UInt position ) =0;

    virtual void setRecompute(UInt position, bool flag){assert(false);}

    //!
    /*!
      adds a new block
        @param Mat block matrix to push
        @param position position of the matrix to which we want to add Mat
     */
    virtual void addToBlock( const matrix_ptrtype& Mat, UInt position );

    //! replaces a block
    /*!
      replaces a block on a specified position in the vector
        @param Mat block matrix to push
        @param index position in the vector
     */
    virtual void replace_matrix( const matrix_ptrtype& Mat, UInt index)=0;


    //! replaces a coupling block
    /*!
      replaces a block on a specified position in the vector
        @param Mat block matrix to push
        @param index position in the vector
     */
    virtual void replace_coupling( const matrix_ptrtype& Mat, UInt index)=0;

    //! runs GlobalAssemble on the blocks
    /*!
      closes and distributes all the matrices before computing the preconditioner
     */
    virtual void GlobalAssemble()=0;

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
    virtual void coupler(map_shared_ptrtype& map,
                         const std::map<ID, ID>& locDofMap,
                         const vector_ptrtype& numerationInterface,
                         const Real& timeStep)=0;



    //! Adds a default coupling block at a specified position
    /*!
      The coupling is handled
      through an augmented formulation, introducing new variables (multipliers). Needs as input: the global map of the problem,
      the two FESpaces of the subproblems to be coupled with their offsets, a std::map holding the two numerations of the
      interface between the two subproblems (the numeration can be different but the nodes must be matching in each
      subdomain), an EpetraVector defined on the multipliers map containing the corresponding dof at the interface (NB: the multipliers map should be constructed from the second numeration in the std::map)
      @param map the map of the global problem
      @param FESpace1 FESpace of the first problem
      @param offset1  offset for the first block in the global matrix
      @param FESpace2 FESpace of the second problem
      @param offset2  offset for the second block in the global matrix
      @param locDofMap std::map with the correspondence between the interface dofs for the two different maps in
      the subproblems
      @param numerationInterface vector containing the correspondence of the Lagrange multipliers with the interface dofs
      @param couplingBlock: flag specifying the block associated with the coupling
     */
    virtual void coupler(map_shared_ptrtype& map,
                          const std::map<ID, ID>& locDofMap,
                          const vector_ptrtype& numerationInterface,
                          const Real& timeStep,
                          UInt couplingBlock
                          )=0;

    //! returns true if the operator is set
    /*!
    */
    virtual bool set()=0;
    //@}

    //@name virtual methods
    //@{
    //!replace a block preconditioner
    /*!
      (only used if the operator is a preconditioner)
     */
    virtual void replace_precs ( const epetra_operator_ptrtype& Mat, UInt index)
    {ERROR_MSG("this method should not be implemented");}

    //!pushes back a block preconditioner
    /*!
      (only used if the operator is a preconditioner)
     */
    virtual void push_back_precs (const epetra_operator_ptrtype& Mat)
    {ERROR_MSG("this method should not be implemented");}
    //@}

    //!replaces a BCHandler
    /*!
      replaces a BCHandler in the vector of bcHandlers at a specified position
      \param bch: input BCHandler to replace
      \param position: position
     */
    virtual void replace_bch(bchandler_ptrtype& bch, UInt position){};

    //!Applies the correspondent boundary conditions to every block
    /*!
      note that this method must be called after blockAssembling(), that sums the coupling conditions to the blocks
      \param time: time
     */
    virtual void applyBoundaryConditions(const Real& time);

    //!Applies the correspondent boundary conditions to a specified block
    /*!
      note that this method must be called after blockAssembling(), that sums the coupling conditions to the blocks
      \param time: time
      \param block: number of the block
     */
    virtual void applyBoundaryConditions(const Real& time, const UInt block);

    //!sums the coupling matrices with the corresponding blocks
    /*!
      Everything (but the boundary conditions) must have been set before calling this
     */
    virtual void blockAssembling(){}


    //!resets the blocks (frees the shared pointers)
    /*!
     */
    virtual void resetBlocks()
    {
        M_blocks.clear();
    }

    //!sets the communicator
    /*!
     */
    virtual void setComm(boost::shared_ptr<Epetra_Comm> comm )
    {
        M_comm = comm;
    }

    //!resets the blocks, boundary conditions, FE spaces.
    /*!
     */
    virtual void reset()
    {
        M_blocks.clear();
        M_bch.clear();
        M_FESpace.clear();
        M_offset.clear();
    }

    //!Applies the robin preconditioners
    /*!
      Applies the robin preconditioners when needed, otherwise does nothing
     */
    virtual void setRobin(matrix_ptrtype& mat, vector_ptrtype& rhs){}


    //! builds the coupling matrix.
    /*!
      Computes the matrix that couples 2 (or more) blocks using an augmented formulation (Lagrange multipliers).
      It adds from 1 to 4 coupling blocks, depending on the flag 'coupling'. This flag is an Int, it works as
      the bash command chmod. It ranges from 1 to 15 for two blocks. The values from 15 to 31 account for the presence
      of a third block: in FSI these three blocks are the fluid block, the structure and the fluid mesh displacement
      blocks. For 'coupling > 31' the method is called recursively.
       \param bigMatrix: the coupling matrix to be built
       \param coupling: flag specifying what wether to consider all the coupling or to neglect one part
       \param problem1: FESpace of the first block
       \param offset1: offset of the first block
       \param problem2: FESpace of the second block
       \param offset2: offset of the second block
       \param locDofMap: std::map with the correspondence between the interface dofs for the two different maps in
       the subproblems
       \param numerationInterface vector containing the correspondence of the Lagrange multipliers with the interface dofs
       \param value value to insert in the coupling blocks
     */
    void couplingMatrix(matrix_ptrtype& bigMatrix,
                        Int coupling,
                        const std::vector<fespace_ptrtype>& problem,
                        const std::vector<UInt>& offset,
                        const std::map<ID, ID>& locDofMap,
                        const vector_ptrtype& numerationInterface,
                        const Real& timeStep=1.e-3,
                        const Real& value=1.); // not working with non-matching grids


    //!applies a Robin preconditioner
    /*!
      (only used if the operator is a particular preconditioner, that is AdditiveSchwarz or AdditiveSchwarzRN)
      applies a Robin preconditioner (mixing the two coupling conditions) to the input matrix, to the rhs, and to the internal
      preconditioner matrix M_blockPrecs.
      \param matrix: input matrix (the matrix of the linear system)
      \param rhsFull: right hand side of the linear system
     */
    void robinPreconditioner(matrix_ptrtype& matrix, vector_ptrtype& rhsFull){}

    //!sets the vector of raw pointer to the BCHandler
    /*!
      each entry of the vector correspond to a block in the same
      position in the vector of matrix pointers M_blocks. An arbitrary number of raw pointers can be passed.
      \param blocks: total number of blocks
     */
    void setConditions(std::vector<bchandler_ptrtype>& vec );


    //!sets the vector of raw pointer to the FESpaces
    /*!
      An arbitrary number of raw pointers can be passed.
      \param blocks: total number of input FESpaces
     */
    void setSpaces    ( std::vector<fespace_ptrtype>& vec );

    //!sets the vector of raw pointer to the offsets of the different blocks
    /*!
      An arbitrary number of raw pointers can be passed.
      \param blocks: total number of input offsets
     */
    void setOffsets    ( UInt blocks, .../*fespace_ptrtype ptr1, fespace_ptrtype ptr2*/);


    //!computes the Robin coupling matrix
    /*!
      Computes a matrix that mixes the  coupling conditions between the blocks:
       [0,0,0,0;0,0,0,alphaf;0,0,0,0;0,alphas,0,0].

      If the coupling
      conditions are of Dirichlet and Neumann type (e.g. continuity of velocity and stress) then the preconditioned
      system will have two Robin conditions instead.
      \param matrix: the input matrix that will hold the robin coupling;
      \param alphaf: parameter multiplying the b.c. of the first block;
      \param alphas: parameter multiplying the b.c. of the second block;
      \param coupling: flag specifying if we want to neglect part of the coupling (used for preconditioning purposes);
      \param map: global map of the whole problem;
      \param FESpace1: FESpace of the first block;
      \param offset1: offset of the first block in the global matrix;
      \param FESpace2: FESpace of the second block;
      \param offset2: offset of the second block in the global matrix;
      \param locDofMap: std::map with the correspondance between the numeration of the interface in the 2 FE spaces.
      \param numerationInterface:  vector containing the correspondence of the Lagrange multipliers with the interface dofs
     */
    void robinCoupling( BlockInterface::matrix_ptrtype& matrix,
                        Real& alphaf,
                        Real& alphas,
                        UInt coupling,
                        const BlockInterface::fespace_ptrtype& FESpace1,
                        const UInt& offset1,
                        const BlockInterface::fespace_ptrtype& FESpace2,
                        const UInt& offset2,
                        const std::map<ID, ID>& locDofMap,
                        const BlockInterface::vector_ptrtype& numerationInterface );

    virtual void push_back_oper( BlockInterface& Oper);

    virtual void push_back_coupling( matrix_ptrtype& coupling)=0;
    //@}

    //!@name Getters
    //@{
    //! returns the vector of pointers to the blocks (by const reference).
    const std::vector<matrix_ptrtype>&    getBlockVector(){return M_blocks;}

    //! returns the vector of pointers to the BCHandlers (by const reference).
    const std::vector<bchandler_ptrtype>& getBChVector(){return M_bch;}

    //! returns the vector of pointers to the FE spaces (by const reference).
    const std::vector<fespace_ptrtype>&   getFESpaceVector(){return M_FESpace;}

    //! returns the vector of the offsets (by const reference).
    const std::vector<UInt>&              getOffsetVector(){return M_offset;}
    //@}

protected:

    //! @name Protected Methods
    //@{

    //!sums the coupling matrix in the specified position with the corresponding block
    /*!
      Everything (but the boundary conditions) must have been set before calling this
     */
    virtual void blockAssembling(const UInt k){}

    //!swaps two boost::shared_ptr. The tamplate argument of the shared_ptr is templated
    /*!
      \param operFrom shared_ptr to be swapped
      \param operTo shared_ptr to be swapped
     */
    template <typename Operator>
    void
    swap(boost::shared_ptr<Operator>& operFrom, boost::shared_ptr<Operator>& OperTo);

    //!swaps two boost::shared_ptr. The tamplate argument of the shared_ptr is templated
    /*!
      \param operFrom shared_ptr to be swapped
      \param operTo shared_ptr to be swapped
     */
    template <typename Operator>
    void
    insert(std::vector<Operator>& operFrom, std::vector<Operator>& OperTo);
    //@}

    //! @name Protected Members
    //@{
    //ComposedOperator<Epetra_Operator>                      M_blocks;
    std::vector<matrix_ptrtype>                                  M_blocks;
    std::vector<bchandler_ptrtype>                               M_bch;
    std::vector<fespace_ptrtype>                                 M_FESpace;
    std::vector<UInt>                                            M_offset;
    vector_ptrtype                                               M_numerationInterface;
    boost::shared_ptr<Epetra_Comm>                               M_comm;
    //Int                                                          M_superCouplingFlag;
    //@}
    //boost::shared_ptr<OperatorPtr>                M_blockPrec;
private:

    //! @name Private Methods
    //@{

    //!private copy constructor: this class should not be copied.
    /*!
      If you need a copy you should implement it, so that it copies the shared pointer one by one, without copying
      the content.
     */
    BlockInterface( matrix_ptrtype& T);
    //@}
};


typedef singleton<factory<BlockInterface,  std::string> >     BlockPrecFactory;


    template <typename Operator>
    void
    BlockInterface::swap(boost::shared_ptr<Operator>& operFrom, boost::shared_ptr<Operator>& operTo)
    {
        boost::shared_ptr<Operator> tmp = operFrom;
        operFrom = operTo;
        operTo = tmp;
    }


    template <typename Operator>
    void
    BlockInterface::insert(std::vector<Operator>& vectorFrom, std::vector<Operator>& vectorTo)
    {
        vectorTo.insert(vectorTo.end(), vectorFrom.begin(), vectorFrom.end());
    }


} // Namespace LifeV

#endif /* BLOCKINTERFACE_H */
