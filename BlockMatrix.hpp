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
    @brief This file contains a class which handles matrices with a block structure.

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 08 Jun 2010
 */

#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H 1

#include <lifemc/lifesolver/BlockInterface.hpp>

namespace LifeV
{

//! BlockMatrix - class which handles matrices with a block structure.
/*!
    @author Paolo Crosetto

    This class handles a matrix with block structure. In particular a vector of shared_ptr
    points to the different blocks, while one matrix contains the coupling part. all the blocks and the coupling
    are summed into the matrix M_globalMatrix for the solution of the linear system. If a BlockMatrix is used as
    a preconditioner by default an algebraic additive Schwarz preconditioner is built on the matrix M_globalMatrix.
 */

class BlockMatrix            : public BlockInterface
{
public:

    //! @name Public Types
    //@{
    typedef  BlockInterface                 super_Type;
    typedef  super_Type::fespacePtr_Type  fespacePtr_Type;
    typedef  super_Type::vector_Type             vector_Type;
    typedef  super_Type::vectorPtr_Type          vectorPtr_Type;
    typedef  super_Type::solverPtr_Type          solverPtr_Type;
    typedef  super_Type::matrix_Type             matrix_Type;
    typedef  super_Type::matrixPtr_Type          matrixPtr_Type;
    typedef  super_Type::epetraOperatorPtr_Type  epetraOperatorPtr_Type;
    typedef  super_Type::mapPtr_Type             mapPtr_Type;
    typedef singleton<factory<BlockMatrix,  std::string> >     Factory;
    //@}


    //! @name Constructor & Destructor
    //@{

    BlockMatrix(UInt coupling):
            super_Type(),
            M_globalMatrix(),
            M_coupling(),
            M_interfaceMap(),
            M_interface(0),
            M_couplingFlag(coupling),
            M_numerationInterface()
    {}

    ~BlockMatrix() {}
    //@}

    //! @name Virtual methods
    //@{
    //! Sets the parameters needed by the class from data file
    /*!
        @param data GetPot object reading the text data file
        @param section string specifying the path in the data file where to find the options for the operator
     */
    virtual void  setDataFromGetPot( const GetPot& data, const std::string& section);

    //! runs GlobalAssemble on the blocks
    /*!
      closes and distributes all the blocks, sums all the blocks together into
      M_globalMatrix.
     */
    virtual void GlobalAssemble();

    //! Computes the coupling matrix
    /*!
      computes the coupling matrix for the system (or for the preconditioner). The coupling is handled
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
     */
    virtual void coupler(mapPtr_Type& map,
                         const std::map<ID, ID>& locDofMap,
                         const vectorPtr_Type& numerationInterface,
                         const Real& timeStep);


    //! Adds a coupling part to the already existing coupling matrix
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
      @param  unused flag kept for compliance with the base class
     */
    void coupler(mapPtr_Type& map,
                 const std::map<ID, ID>& locDofMap,
                 const vectorPtr_Type& numerationInterface,
                 const Real& timeStep,
                 UInt /*flag*/ );

    //! returns true if the operator has at least one block
    /*!
    */
    virtual bool  set() {return (bool) super_Type::M_blocks.size();}

    //@}
    //! @name Public methods
    //@{

    //! Solves the preconditioned linear system (used only when blockMatrix is used as a preconditioner)
    /*!
      Provided the linear solver and the right hand side this method computes the solution and returns it into
      the result vector.
      @param rhs right hand side of the linear system
      @param result output result
      @param linearSolver the linear system
    */
    int   solveSystem( const vector_Type& rhs, vector_Type& step, solverPtr_Type& linearSolver);


    //! pushes a block at the end of the vector
    /*!
      adds a new block
        @param Mat block matrix to push
        @param recompute flag stating wether the preconditioner for this block have to be recomputed at every time step.
        In this case it is not used since it is equal to the boolean specifying wether the whole preconditioner must be
        recomputed or not.
     */
    void  push_back_matrix( const matrixPtr_Type& Mat, bool /*recompute*/);

    //! replaces a block
    /*!
      replaces a block on a specified position in the vector
        @param Mat block matrix to push
        @param index position in the vector
     */
    void  replace_matrix( const matrixPtr_Type& Mat, UInt index);


    //! replaces the coupling matrix without copies
    /*!
      this method is deprecated, it is implemented for compatibility with the base class
      \param Mat replacing matrix
     */
    void replace_coupling( const matrixPtr_Type& /*Mat*/, UInt /*index*/)
    {
        // not used for matrices (only for preconditioners)
        /*M_coupling = Mat;*/
    }

    //! never used within this class
    /*!
      runs assert(false) when called. This method is used only when the preconditioner is a composed operator
     */
    void  replace_precs( const epetraOperatorPtr_Type& Mat, UInt index);

    //!sums the coupling matrix in the specified position with the global matrix
    /*!
      Everything (but the boundary conditions) must have been set before calling this
     */
    void blockAssembling( );

    //! returns the global matrix, with all the blocks and the coupling parts
    /*!
    */
    matrixPtr_Type& matrix( ){return M_globalMatrix;}

    //! multiplies the whole system times a matrix
    /*!
      applies a matrix robinCoupling to the system matrix M_globalMatrix, to the rpeconditioner prec passed as input,
      to the rhs passed as input. This is useful e.g. when we want to substitute the Dirichlet-Neumann coupling with a
      Robin-Robin one.
      \param robinCoupling matrix that multiplies the system
      \param prec preconditioner matrix
      \param rhs right hand side of the system
    */
    void applyPreconditioner(   matrixPtr_Type robinCoupling, matrixPtr_Type prec, vectorPtr_Type& rhs);

    //! multiplies two matrices
    /*!
      multiplies on the left the matrix prec times the matrix robinCoupling
      \param robinCoupling the matrix that multiplies
      \param prec the matrix multiplied (modified)
     */
    void applyPreconditioner( const matrixPtr_Type robinCoupling, matrixPtr_Type& prec );


    //! applies a matrix to the system
    /*!
      multiplies on the left the system matrix M_globalMatrix and the rhs times the matrix matrix.
      \param matrix the preconditioning matrix
      \param rhsFull the right hand side of the system
     */
    void applyPreconditioner( const matrixPtr_Type matrix, vectorPtr_Type& rhsFull);

    //!creates the map for the coupling
    /*!
      This method create the coupling map needed for two blocks. The coupling is handled adding a number of Lagrange
      multipliers equal to the interface degrees of freedom. This method generates the map for the vector of the Lagrange
      multipliers.
      \param interfaceMap: the map of the interface (not contiguous, the numeration of this map has "holes")
      \param locDofMap: the map between the two (e. g. in FSI the fluid (first) and the solid (second)) numerations of
      the interface
      \param subdomainMaxId: The maximum ID of the sub-map built from interfaceMap (is the number of dof in the subdomain
      used to buid the interfaceMap divided by the number of dimensions). The interface map used to build the
      vector M_numerationInterface contains once every node, not every interface degree of freedom (that is nDimensions *
      interfaceNodes)
      \param epetraWorldComm: The communicator
     */
    void createInterfaceMap( const EpetraMap& interfaceMap,
                             const std::map<ID, ID>& locDofMap,
                             const UInt subdomainMaxId,
                             const boost::shared_ptr<Epetra_Comm> epetraWorldComm );


    //! applies the b.c. to every block
    void applyBoundaryConditions(const Real& time);

    //! applies the b.c. to a specific block
    void applyBoundaryConditions(const Real& time, const UInt block );

    //! applies all the b.c.s to the global matrix, updates the rhs
    void applyBoundaryConditions(const Real& time, vectorPtr_Type& rhs);

    //! applies the b.c. to the a specific block, updates the rhs
    void applyBoundaryConditions(const Real& time, vectorPtr_Type& rhs, const UInt block);

    //! adds a block to the coupling matrix
    /*!
      @param Mat: block added
      @param position of the block, unused here because there is only one coupling matrix
     */
    void addToCoupling( const matrixPtr_Type& Mat, UInt /*position*/);

    //! adds a block to the coupling matrix
    /*!
      This method is specific for the BlockMatrix class, it is used e.g. to add the shape derivatives block in FSI
      to the global matrix
      @param Mat: block matrix to be added.
     */
    void addToGlobalMatrix( const matrixPtr_Type& Mat)
    {
        matrixPtr_Type tmp(new matrix_Type(M_globalMatrix->map()));
        *tmp += *M_globalMatrix;
        *tmp += *Mat;
        tmp->globalAssemble();
        M_globalMatrix = tmp;
    }

    //! adds a coupling block to the coupling matrix
    /*!
      This method is kept for compatibility with the base class. It calls the method addToCoupling.
     */
    void push_back_coupling( matrixPtr_Type& coupling)
    {
        addToCoupling(coupling, 0);
    }
    //@}

    //!@name Get Methods
    //@{
    //! returns the dimension of the interface
    /*!
      NOTE: it has to be multiplied times nDimensions to get the number of interface dofs
     */
    UInt interface(){return M_interface;}

    //! returns the map built for theLagrange multipliers
    mapPtr_Type interfaceMap() const { return M_interfaceMap; }

    //! returns the numeration of the interface
    /*!
      \param numeration: output vector
     */
    void numerationInterface( vectorPtr_Type& numeration ) { numeration =  M_numerationInterface; }

    //@}

    //!@name Factory Method
    //@{
    static BlockMatrix*    createAdditiveSchwarz()
    {
        return new BlockMatrix(15);
    }
    //@}


protected:

    //! @name Protected Members
    //@{
    matrixPtr_Type                              M_globalMatrix;
    matrixPtr_Type                              M_coupling;
    mapPtr_Type                          M_interfaceMap;
    UInt                                        M_interface;
    //@}

private:

    //! @name Private Members
    //@{
    const UInt                                  M_couplingFlag;
    vectorPtr_Type                              M_numerationInterface;
    //@}
};


} // Namespace LifeV

#endif /* BLOCKMATRIX_H */
