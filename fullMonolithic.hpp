/* -*- mode: c++ -*-

   This file is part of the LifeV library

   Author(s): Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   Date: 2008-09-18

   Copyright (C) 2008

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file fullMonolithic.hpp
   \author crosetto <Paolo Crosetto>
   \date 18/09/2008
*/
#ifndef _FULLMONOLITHIC_HPP
#define _FULLMONOLITHIC_HPP

#include <lifemc/lifesolver/Monolithic.hpp>

namespace LifeV
{
class Epetra_FullMonolithic;

/**
 * Class handling the nonlinear monolithic solver for FSI problems. The (exact or inexact)
 * Newton algorithm is used to solve the nonlinearity.
 * The block structure of the jacobian matrix is
 *\f$\left(\begin{array}{ccc}
 C&B&S\\
 D&N&0\\
 0&E&H
 \end{array}\right)\f$
 * where \f$N\f$ represents the solid block, \f$C\f$ the fluid block, \f$H\f$ is the harmonic extension block,
 * while the extra
 * diagonal blocks represent the coupling. The implementation of the stress continuity coupling condition
 * is obtained by means of an augmented formulation.
 * Different possible preconditioners are implemented.
 */
class fullMonolithic : public Monolithic
{
public:

    typedef Monolithic super;

    fullMonolithic();
    ///constructor


    virtual ~fullMonolithic(){}
    ///destructor

    /**
       constructs the matrix handling the coupling and sums it to matrix
       \param matrix: output matrix
       \param coupling:   flag handling the coupling of the different blocks. Chosing properly it's value, ranging from 0 to 31, one can decide which couplings to keep and which to neglect (it works like the chmod command in bash)
    */

    void                        setUp( const GetPot& dataFile );

    /**
       creates FEspace
    */
    void                       setupFEspace();


    void                       setupFluidSolid();


    void                       setDataFromGetPot( GetPot const& data_file );

    /**
       updates the meshmotion, advances of a time step
       \param displacement: solution
    */
    void                        updateSystem();

    void                      couplingMatrix(matrix_ptrtype & matrix,
                                             int coupling=31);
    /**
       applies the velocity continuity coupling to the right hand side of the equation
    */
    //void                      couplingRhs(vector_ptrtype rhs);

    void                        buildSystem ();

    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: fluid domain displacement solution
       \param iter: current nonLinRichardson (block Gauss Seidel for the tangent system) iteration
    */
    void                      evalResidual(vector_type&        res,
                                           const vector_type& _sol,
                                           const UInt          _iter);

    int                        setupBlockPrec(vector_type& rhs);


    //    void updateSystem(const vector_type& _sol);


    const EpetraMap&            mapWithoutMesh() const {return *M_mapWithoutMesh;}
    //    vector_type& solution(){return *M_un;}


    /**
       solves the Jacobian system
       \param _muk: output, solution at the current block GS step
       \param _res: linear residual
       \param _linearRelTol: not used
    */
    void                        solveJac(vector_type       &_muk,
                                         const vector_type &_res,
                                         const Real       _linearRelTol);

    //! @getters
    //!{
    const matrix_ptrtype        getMatrixPtr() const {return this->M_monolithicMatrix;}

    const vector_ptrtype        uk()      const {return M_uk;}

    //const vector_type&          meshVel() const;

    const vector_type&          meshDisp()const
    {
        return this->M_meshMotion->dispOld();
    }
    //}



    void                        initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                                            FSIOperator::solid_type::value_type::Function const& p0,
                                            FSIOperator::solid_type::value_type::Function const& d0,
                                            FSIOperator::solid_type::value_type::Function const& w0,
                                            FSIOperator::solid_type::value_type::Function const& df0 );

    void                        initializeMesh(vector_ptrtype fluid_dispOld);


private:

    void                        initialize( vector_type const& u0, vector_type const& p0, vector_type const& d0, vector_type const& df0);

    /**
       calculates the terms due to the shape derivatives on the rhs of the monolithic system rhsShapeDerivatives
       given the mesh increment deltaDisp.
       \param rhsShapeDerivatives: output. Shape derivative terms.
       \param meshDeltaDisp: input. Mesh displacement increment.
    */
    void                        shapeDerivatives(matrix_ptrtype sdMatrix, const vector_type& sol,  bool fullImplicit, bool convectiveTerm);


    boost::shared_ptr<EpetraMap>         M_mapWithoutMesh;
    vector_ptrtype                       M_uk;
    vector_ptrtype                       M_meshVel;
    bool                                 M_domainVelImplicit;
    bool                                 M_convectiveTermDer;
};

}
#endif
