/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/**
 * \file Monolithic.hpp
 * \author Paolo Crosetto
 * \date 13-09-2008
 */
#ifndef _MONOLITHIC_HPP
#define _MONOLITHIC_HPP


#include <lifemc/lifesolver/FSIOperator.hpp>
//#include <Epetra_IntVector.h>

namespace LifeV
{

class Monolithic : public FSIOperator
{

public:


    typedef FSIOperator                            super;
    typedef FSIOperator::fluid_type::value_type::matrix_type   matrix_type;
    typedef boost::shared_ptr<matrix_type>                    matrix_ptrtype;


    Monolithic();
    ///constructor
    virtual ~Monolithic();
    ///destructor


    void   evalResidual(vector_type&        res,
                        const vector_type& _sol,
                        const int          _iter);
    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: monolithic solution
       \param iter: current nonLinRichardson (Newton) iteration
    */


    virtual void   solveJac(vector_type&       _muk,
                    const vector_type& _res,
                    const double       _linearRelTol);

    /**
       solves the Jacobian system
       \param _muk: output, solution at the current Newton step
       \param _res: nonlinear residual
       \param _linearRelTol: not used
    */

    virtual void updateSystem(const vector_type& _sol);

    /**
       updates the meshmotion, advances of a time step
       \param _sol: solution
    */


    void setDataFromGetPot( GetPot const& data );

    void monolithicToInterface(    vector_type& _lambdaSolid, const vector_type& _sol) ;

    /**
       restricts a vector with a monolithic map on the solid interface map
       \param _lambdasolid: vector on the solid interface
       \param _disp: monolithic vector
    */

    void monolithicToSolid(const vector_type& _disp, vector_type& _dispSolid);
    /**
       restricts a vector with a monolithic map on the solid map
       \param _dispSolid: vector on the solid domain
       \param _disp: monolithic vector
    */

    void setDispSolid(const vector_type& sol);
    /**
       sets the vector M_solid->dispSolid() to the monolithic solution M_solid->disp() in the solid nodes, to 0 in the fluid nodes
    */

    void monolithicToFluid(const vector_type& _disp, vector_type& _dispFluid);
    /**
        restricts a vector with a monolithic map on the fluid map
       \param _dispFluid: vector on the fluid domain
       \param _disp: monolithic vector
    */

    void iterateMesh(const vector_type& disp);
    /**
       iterates the mesh
       \param disp: monolithic solution
    */

    virtual void buildSystem();
    /**
      builds the constant part of the monolithic matrix
    */

    virtual void setup();
    /**
       assigns each mesh partition to the corresponding processor, builds the monolithic map
    */
    void shapeDerivatives(vector_ptrtype rhsShapeDerivatives, vector_ptrtype meshDeltaDisp, const vector_type& sol);
    /**
       calculates the terms due to the shape derivatives on the rhs of the monolithic system rhsShapeDerivatives
       given the mesh increment deltaDisp.
       \param rhsShapeDerivatives: output. Shape derivative terms.
       \param meshDeltaDisp: input. Mesh displacement increment.
    */

    //    bool isTangentProblem(){return M_isTangentProblem;}

    //    void setIsTangentProblem(bool isTangentProblem){M_isTangentProblem=isTangentProblem;}

    vector_type getRhsShapeDerivatives(){return *M_rhsShapeDerivatives;}

    //    EpetraMap& monolithicMap(){return M_monolithicMap;}
    const boost::shared_ptr<EpetraMap>& monolithicMap() const {return M_monolithicMap;}

    //    void addShapeDerivativesTerm(vector_ptrtype rhs, vector_ptrtype SDTerm);
    void setFluidPrecBC         (fluid_bchandler_type bc_fluid_prec);
    void setSolidPrecBC         (solid_bchandler_type bc_solid_prec);
    const UInt dimInterface() const {return nDimensions*M_interface ;}


protected:

        void addDiagonalEntries(Real entry, matrix_ptrtype matrix, EpetraMap& Map, UInt offset=0, bool replacs=false);
    /**
       \small adds a constant scalar entry to a diagonal block of a matrix
       \param entry: entry
       \param matrix: matrix
       \param Map: map specifying the location where to add diagonal entries
   */

        void addDiagonalEntries(Real entry, matrix_ptrtype matrix, std::map<ID, ID> const& Map, UInt offset=0, bool replacs=false);
    /**
       \small adds a constant scalar entry to a diagonal block of a matrix
       \param entry: entry
       \param matrix: matrix
       \param Map: map specifying the location where to add diagonal entries
   */

    virtual void couplingMatrix(matrix_ptrtype & bigMatrix, bool couplingSolid = true);
    /**
      \small  to the monolithic matrix.
      \param bigMatrix: matrix to be modified
    */
    void couplingRhs( vector_ptrtype rhs , vector_ptrtype un );
    /**
       \small adds the part due to coupling to the rhs
       \param rhs: right hand side
       \param un: current solution
    */

    void robinCoupling(matrix_ptrtype & IdentityMatrix, Real& alphaf, Real& alphas/*, bool inverse = false*/);
    /**
       \small adds the matrix [0,0,0,0;0,0,0,alphaf;0,0,0,0;0,alphas,0,0] to IdentityMatrix
       \param IdentityMatrix: a matrix, if it is an identity the method builds the Robin-Robin linear combination of boundary conditions
       \param alphaf: coefficient
       \param alphas: coefficient
       \param inverse if true the inverse of the linear combination is computed. the output matrix is [];
    */

    void robinCouplingInv(matrix_ptrtype & IdentityMatrix, Real& alphaf, Real& alphas/*, bool inverse*/); // not working with non-matching grids

    void zeroBlock( matrix_ptrtype matrixPtr,vector_type& colNumeration , const std::map<ID, ID>& map, UInt rowOffset=0, UInt colOffset=0);

    boost::shared_ptr<EpetraMap>                      M_monolithicMap;
    matrix_ptrtype                                    M_couplingMatrix;
    UInt                                              M_interface;
    EpetraMap                                         M_interfaceMap;///the solid interface map
    boost::shared_ptr<vector_type>                    M_beta;///ALE relative velocity (u-w)
    matrix_ptrtype                                    M_monolithicMatrix;
    //    bool                                              firstIter;

private:

    int                                               M_updateEvery;
    bool                                              M_semiImplicit;
    boost::shared_ptr<vector_type>                    M_numerationInterface;
    boost::shared_ptr<vector_type>                    M_rhsShapeDerivatives;
    /**
vector containing the numeration of the coupling part of the monolithic map
     */

    boost::shared_ptr<vector_type>                    M_rhsNew;
    bool                                              M_fullMonolithic;
    matrix_ptrtype                                    M_bigPrecPtr;
    bool                                              M_recomputePrec;
    fluid_bchandler_type                              M_BCh_uPrec;
    solid_bchandler_type                              M_BCh_dPrec;
    //bool                                              M_splitPrec;
    Real                                              M_entry;
    int                                               M_DDBlockPrec;
    matrix_ptrtype                                    M_robinCoupling;
    matrix_ptrtype                                    M_robinCouplingInv;
    Real                                              M_alphaf;
    Real                                              M_alphas;
    bool                                              M_isDiagonalBlockPrec;
    bool                                              M_diagonalScale;
    UInt                                              M_offset;
    UInt                                              M_solidAndFluidDim;
    //    boost::shared_ptr<Epetra_IntVector>               M_numerationInterfaceInt;
    //    bool                                              M_isTangentProblem;


};

}
#endif
