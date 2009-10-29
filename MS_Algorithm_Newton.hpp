/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-10-26

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_Algorithm_Newton.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-10-26
 */

#ifndef __MS_Algorithm_Newton_H
#define __MS_Algorithm_Newton_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>

#include <lifemc/lifesolver/MS_Algorithm.hpp>

namespace LifeV {

//! MS_Algorithm_Newton_Operator - The Epetra_Operator implementation for the Newton Algorithm
/*!
 *  The MS_Algorithm_Newton_Operator is an implementation of Epetra_Operator
 *  for the Newton algorithm.
 *
 *  @author Cristiano Malossi
 */
class MS_Algorithm_Newton_Operator : public Epetra_Operator
{
public:

    typedef MS_PhysicalModel::VectorType          VectorType;
    typedef MS_PhysicalModel::Vector_ptrType      Vector_ptrType;

    //! @name Constructors & Destructor
    //@{

    MS_Algorithm_Newton_Operator():
        M_algorithm (),
        M_jacobianProduct()
    {}

    virtual ~MS_Algorithm_Newton_Operator() {};

    //! @name Set Methods
    //@{

    //! SetUseTranspose
    void SetAlgorithm( const MS_Algorithm *algorithm )
    {
        M_algorithm = algorithm;

        M_jacobianProduct.reset( new EpetraVector( *M_algorithm->GetCouplingVariables() ) );
        *M_jacobianProduct = 0.0;
    }

    //! SetUseTranspose
    void InitializeJacobianProduct( const MS_Algorithm *algorithm )
    {
        M_algorithm = algorithm;
    }

    //! SetUseTranspose
    int SetUseTranspose ( bool /*UseTranspose*/ )
    {
        std::cout << "********* MS_Algorithm_Newton_Operator : transpose not available\n";
        return (EXIT_FAILURE);
    }

    //@}


    //! @name Get Methods
    //@{

    const char* Label() const
    {
        return "MS_Algorithm_Newton_Operator";
    }

    //! UseTranspose
    bool UseTranspose() const
    {
        return false;
    }

    //! HasNormInf
    bool HasNormInf() const
    {
        return false;
    }

    //! Comm
    const Epetra_Comm& Comm() const
    {
        return *( M_algorithm->GetCommunicator() );
    }

    //! OperatorDomainMap
    const Epetra_Map& OperatorDomainMap() const
    {
        return M_algorithm->GetCouplingResiduals()->getEpetra_Map();
    }

    //! OperatorRangeMap
    const Epetra_Map& OperatorRangeMap() const
    {
        return M_algorithm->GetCouplingResiduals()->getEpetra_Map();
    }

    //@}


    //! @name Methods
    //@{

    int Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
    {
        VectorType MyXCopy( X, M_jacobianProduct->getMap_ptr(), M_jacobianProduct->getMaptype() );
        M_algorithm->GetMultiScaleProblem()->ExportJacobianProduct( MyXCopy, *M_jacobianProduct );

        Y = M_jacobianProduct->getEpetraVector();

        //Y.scale( 0.0 );
        //Y.update
        //Y = M_jacobianProduct->getEpetraVector();

        return ( EXIT_SUCCESS );
    }

    int ApplyInverse( const Epetra_MultiVector& /*X*/, Epetra_MultiVector& /*Y*/ ) const
    {
        std::cout << "********* MS_Algorithm_Newton_Operator : inverse not available\n";
        return ( EXIT_FAILURE );
    }

    double NormInf() const
    {
        std::cout << "********* MS_Algorithm_Newton_Operator : NormInf not available\n";
        return ( EXIT_FAILURE );
    }

    //@}

private:

    const MS_Algorithm*                   M_algorithm;
    Vector_ptrType                        M_jacobianProduct;
};


//! MS_Algorithm_Newton - The MultiScale Algorithm implementation of Newton
/*!
 *  The MS_Algorithm_Newton is an implementation of MS_Algorithm
 *  which implements the Newton method.
 *
 *  @author Cristiano Malossi
 */
class MS_Algorithm_Newton : public virtual MS_Algorithm
{
public:

    typedef MS_Algorithm                  super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Algorithm_Newton();

    //! Copy constructor
    /*!
     * \param algorithm - MS_Algorithm_Newton
     */
    MS_Algorithm_Newton( const MS_Algorithm_Newton& algorithm );

    //! Destructor
    ~MS_Algorithm_Newton() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param algorithm - MS_Algorithm_Newton
     */
    MS_Algorithm_Newton& operator=( const MS_Algorithm_Newton& algorithm );

    //@}


    //! @name MultiScale PhysicalModel Virtual Functions
    //@{

    //! Setup the data of the algorithm
    void SetupData( const GetPot& DataFile );

    //! Perform sub-iteration on the couplings
    void SubIterate( void );

    //! Display some information about the algorithm
    void ShowMe( void );

    //@}

protected:

    //! @name Protected Methods
    //@{

    //@}

    SolverTrilinos                           M_solver;
    MS_Algorithm_Newton_Operator             M_operator;

};

//! Factory create function
inline MS_Algorithm* createNewton()
{
    return new MS_Algorithm_Newton();
}

} // Namespace LifeV

#endif /* __MS_Algorithm_Newton_H */
