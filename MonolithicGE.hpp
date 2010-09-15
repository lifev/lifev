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
    @date 26 Jul 2010

    A more detailed description of the file (if necessary)
 */

#ifndef MONOLITHICGE_H
#define MONOLITHICGE_H 1

#include <lifemc/lifesolver/BlockMatrix.hpp>
#include <lifemc/lifesolver/BlockMatrixRN.hpp>
#include <lifemc/lifesolver/ComposedDN.hpp>
#include <lifemc/lifesolver/ComposedNN.hpp>
#include <lifemc/lifesolver/ComposedDNND.hpp>

#include <life/lifecore/life.hpp>

#include <lifemc/lifesolver/Monolithic.hpp>

namespace LifeV {

//! monolithicGE - Monolithic Geometry-Explicit solver
/*!
  \include ../../testsuite/test_monolithic/fluidstructure.dox
  @author Paolo Crosetto
  @see \ref CDFQ


 Important parameters to set properly in the data file:
 - useShapeDerivatives: MUST be false, because in the GE approach the geometry is explicit;
 - domainVelImplicit: MUST be false, because in the GE approach the geometry is explicit;
 - convectiveTermDer: false if the convective term is linearized (\f$u^{n+1}\nabla(u^n-w^n)\f$),
 otherwise it can be either true (if we use the Newton method to solve the convective term nonlinearity) or false
 (fixed-point method). For the GCE must be false;
 - semiImplicit:  if true only one iteration of the nonlinear solver is performed. Otherwise
 the nonlinear iterations continue up to the specified tolerance. Set it to true for the GCE;
 - method: can be either monolithicGE, monolithicGI if the geometry is treated respectively explicitly or implicitly,
 or exactJacobians, fixedPoint for partitioned strategies;
 - blockOper: specifies the matrix type to be used for the linear system: if AdditiveSchwarz, the matrix is the standard
 ine for GE; if AdditiveSchwarzRN the coupling blocks are of Robin type instead of Dirichlet and Neumann. The parameters
 for the Robin coupling are alphaf and alphas in the data file. NOTE: this method has currently been tested only for
 alphas=0.
 - DDBlockPrec: specifies the possible preconditioners to use. Can be: AdditiveSchwarz, ComposedDN, ComposedDN2,
 ComposedNN, ComposedDNND.

 */
class MonolithicGE : public Monolithic
{
public:

    typedef Monolithic super;
    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MonolithicGE():
        super()
    {}

    //! Destructor
    ~MonolithicGE()
    {}

    //@}



    //! @name Methods
    //@{

    void setupFluidSolid();

    void setupDOF( void );

    void setupSystem( );

    void assembleFluidBlock(UInt iter);

    void getSolution                  (vector_ptrtype& sol){sol = M_un;}


    //@}

private:

    void createOperator( std::string& operType )
    {
        M_monolithicMatrix.reset(BlockMatrix::Factory::instance().createObject( operType ));
    }

    static bool reg;
    //@}
};

} // Namespace LifeV

#endif /* MONOLITHICGCE_H */
