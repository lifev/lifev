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
  @brief Base class for ionic models

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 04 - 2014
 */

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>


namespace LifeV
{
// ===================================================
//! Constructors
// ===================================================
ElectroIonicModel::ElectroIonicModel() :
    M_numberOfEquations (0),
    M_numberOfGatingVariables (0),
    M_restingConditions (),
    M_membraneCapacitance (1.),
    M_appliedCurrent    (0.),
    M_appliedCurrentPtr(),
    M_pacingProtocol ()
{
}

ElectroIonicModel::ElectroIonicModel ( int n ) :
    M_numberOfEquations (n),
    M_numberOfGatingVariables (0),
    M_restingConditions (std::vector<Real> (M_numberOfEquations, 0.0) ),
    M_membraneCapacitance (1.),
    M_appliedCurrent    (0.),
    M_appliedCurrentPtr(),
    M_pacingProtocol ()
{
}

ElectroIonicModel::ElectroIonicModel ( int n, int g ) :
    M_numberOfEquations (n),
    M_numberOfGatingVariables (g),
    M_restingConditions (std::vector<Real> (M_numberOfEquations, 0.0) ),
    M_membraneCapacitance (1.),
    M_appliedCurrent    (0.),
    M_appliedCurrentPtr(),
    M_pacingProtocol ()
{
}

ElectroIonicModel::ElectroIonicModel ( const ElectroIonicModel& Ionic ) :
    M_numberOfEquations ( Ionic.Size() ),
    M_numberOfGatingVariables ( Ionic.numberOfGatingVariables() ),
    M_restingConditions ( Ionic.restingConditions() ),
    M_membraneCapacitance ( Ionic.M_membraneCapacitance ),
    M_appliedCurrent    ( Ionic.M_appliedCurrent ),
    M_pacingProtocol (Ionic.M_pacingProtocol)
{
    if (Ionic.M_appliedCurrentPtr)
    {
        M_appliedCurrentPtr.reset (new vector_Type (*Ionic.M_appliedCurrentPtr) );
    }
}

// ===================================================
//! Methods
// ===================================================
ElectroIonicModel& ElectroIonicModel::operator = ( const ElectroIonicModel& Ionic )
{
    M_numberOfEquations = Ionic.M_numberOfEquations;
    M_numberOfGatingVariables = Ionic.M_numberOfGatingVariables;
    M_restingConditions = Ionic.M_restingConditions;
    M_membraneCapacitance = Ionic.M_membraneCapacitance;
    M_appliedCurrent = Ionic.M_appliedCurrent;
    if (Ionic.M_appliedCurrent)
    {
        M_appliedCurrentPtr = Ionic.M_appliedCurrentPtr;
    }
    M_pacingProtocol = Ionic.M_pacingProtocol;

    return      *this;
}

std::vector< std::vector<Real> > ElectroIonicModel::getJac (const std::vector<Real>& v, Real h)
{
    std::vector< std::vector<Real> > J ( M_numberOfEquations, std::vector<Real> (M_numberOfEquations, 0.0) );
    std::vector<Real> f1 (M_numberOfEquations, 0.0);
    std::vector<Real> f2 (M_numberOfEquations, 0.0);
    std::vector<Real> y1 (M_numberOfEquations, 0.0);
    std::vector<Real> y2 (M_numberOfEquations, 0.0);

    for (int i = 0; i < M_numberOfEquations; i++)
    {
        for (int j = 0; j < M_numberOfEquations; j++)
        {
            y1[j] = v[j] + ( (double) (i == j) ) * h;
            y2[j] = v[j] - ( (double) (i == j) ) * h;
        }
        this->computeRhs (y1, f1);
        f1[0] += M_appliedCurrent;
        this->computeRhs (y2, f2);
        f2[0] += M_appliedCurrent;

        for (int j = 0; j < M_numberOfEquations; j++)
        {
            J[j][i] = (f1[j] - f2[j]) / (2.0 * h);
        }
    }

    return J;
}

MatrixEpetra<Real> ElectroIonicModel::getJac (const vector_Type& v, Real /*h*/)
{
    matrix_Type J (v.map(), M_numberOfEquations, false);
    //  vector<Real> f1(M_numberOfEquations,0.0);
    //  vector<Real> f2(M_numberOfEquations,0.0);
    //  vector< vector<Real> > df ( M_numberOfEquations, vector<Real> (M_numberOfEquations,0.0) );
    //  vector<Real> y1(M_numberOfEquations,0.0);
    //  vector<Real> y2(M_numberOfEquations,0.0);
    //  const Int* k = v.blockMap().MyGlobalElements();
    //
    //  int* Indices = new int[M_numberOfEquations];
    //  double* Values =  new double[M_numberOfEquations];
    //  for(int i=0; i<M_numberOfEquations; i++)
    //      Indices[i] = i;
    //
    //  for(int i=0; i<M_numberOfEquations; i++)
    //  {
    //      for(int j=0; j<M_numberOfEquations; j++)
    //      {
    //          y1[j] = v[ k[j] ] + ((double)(i==j))*h;
    //          y2[j] = v[ k[j] ] - ((double)(i==j))*h;
    //      }
    //      this->computeRhs(y1, f1);
    //      this->computeRhs(y2, f2);
    //
    //      for(int j=0; j<M_numberOfEquations; j++)
    //          df[j][i] = (f1[j]-f2[j])/(2.0*h);
    //  }
    //
    //
    //  for(int i=0; i<M_numberOfEquations; i++)
    //  {
    //      for(int j=0; j<M_numberOfEquations; j++)
    //          Values[j] = df[i][j];
    //      J.matrixPtr()->InsertGlobalValues (i, M_numberOfEquations, Values, Indices);
    //  }
    //
    //  J.globalAssemble();
    //
    //  delete[] Indices;
    //  delete[] Values;

    return J;
}


void ElectroIonicModel::computeGatingRhs (   const std::vector<vectorPtr_Type>& v,
                                             std::vector<vectorPtr_Type>& rhs )
{

    int nodes = ( * (v.at (1) ) ).epetraVector().MyLength();


    std::vector<Real>   localVec ( M_numberOfEquations, 0.0 );
    std::vector<Real>   localRhs ( M_numberOfEquations - 1, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = ( * (v.at (1) ) ).blockMap().GID (k);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            localVec.at (i) = ( * ( v.at (i) ) ) [j];
        }

        if (M_appliedCurrentPtr)
        {
            M_appliedCurrent = (*M_appliedCurrentPtr) [j];
        }
        else
        {
            M_appliedCurrent = 0.0;
        }

        computeGatingRhs ( localVec, localRhs );

        for ( int i = 1; i < M_numberOfEquations; i++ )
        {
            ( * ( rhs.at (i) ) ) [j] =  localRhs.at (i - 1);
        }

    }

}

void ElectroIonicModel::computeNonGatingRhs (   const std::vector<vectorPtr_Type>& v,
                                                std::vector<vectorPtr_Type>& rhs )
{

    int nodes = ( * (v.at (1) ) ).epetraVector().MyLength();


    std::vector<Real>   localVec ( M_numberOfEquations, 0.0 );
    int offset = 1 + M_numberOfGatingVariables;
    std::vector<Real>   localRhs ( M_numberOfEquations - offset, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = ( * (v.at (1) ) ).blockMap().GID (k);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            localVec.at (i) = ( * ( v.at (i) ) ) [j];
        }

        if (M_appliedCurrentPtr)
        {
            M_appliedCurrent = (*M_appliedCurrentPtr) [j];
        }
        else
        {
            M_appliedCurrent = 0.0;
        }

        computeNonGatingRhs ( localVec, localRhs );

        for ( int i = offset; i < M_numberOfEquations; i++ )
        {
            ( * ( rhs.at (i) ) ) [j] =  localRhs.at (i - offset);
        }

    }

}


void ElectroIonicModel::computeRhs (   const std::vector<vectorPtr_Type>& v,
                                       std::vector<vectorPtr_Type>& rhs )
{

    int nodes = ( * (v.at (1) ) ).epetraVector().MyLength();


    std::vector<Real>   localVec ( M_numberOfEquations, 0.0 );
    std::vector<Real>   localRhs ( M_numberOfEquations, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = ( * (v.at (1) ) ).blockMap().GID (k);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            localVec.at (i) = ( * ( v.at (i) ) ) [j];
        }


        if (M_appliedCurrentPtr)
        {
            M_appliedCurrent = (*M_appliedCurrentPtr) [j];
        }
        else
        {
            M_appliedCurrent = 0.0;
        }
        computeRhs ( localVec, localRhs );
        addAppliedCurrent (localRhs);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            ( * ( rhs.at (i) ) ) [j] =  localRhs.at (i);
        }

    }

}

void ElectroIonicModel::computePotentialRhsICI (   const std::vector<vectorPtr_Type>& v,
                                                   std::vector<vectorPtr_Type>& rhs,
                                                   matrix_Type&                    massMatrix  )
{
    int nodes = ( * (v.at (0) ) ).epetraVector().MyLength();

    ( * ( rhs.at (0) ) ) *= 0.0;
    std::vector<Real>   localVec ( M_numberOfEquations, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = ( * (v.at (0) ) ).blockMap().GID (k);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            localVec.at (i) = ( * ( v.at (i) ) ) [j];
        }

        if (M_appliedCurrentPtr)
        {
            M_appliedCurrent = (*M_appliedCurrentPtr) [j];
        }
        else
        {
            M_appliedCurrent = 0.0;
        }

        ( * ( rhs.at (0) ) ) [j] =  computeLocalPotentialRhs ( localVec ) + M_appliedCurrent;

    }

    ( * ( rhs.at (0) ) ) = massMatrix * ( * ( rhs.at (0) ) );

}


void ElectroIonicModel::computePotentialRhsSVI (   const std::vector<vectorPtr_Type>& v,
                                                   std::vector<vectorPtr_Type>& rhs,
                                                   FESpace<mesh_Type, MapEpetra>& uFESpace )
{

    std::vector<Real> U (M_numberOfEquations, 0.0);
    Real I (0.0);
    ( * ( rhs.at (0) ) ) *= 0.0;

    std::vector<vectorPtr_Type>      URepPtr;
    for ( int k = 0; k < M_numberOfEquations; k++ )
    {
        URepPtr.push_back ( vectorPtr_Type ( new VectorEpetra (  * ( v.at (k) )     , Repeated ) ) );
    }

    VectorEpetra    IappRep ( *M_appliedCurrentPtr, Repeated );

    std::vector<elvecPtr_Type>      elvecPtr;
    for ( int k = 0; k < M_numberOfEquations; k++ )
    {
        elvecPtr.push_back ( elvecPtr_Type ( new VectorElemental (  uFESpace.fe().nbFEDof(), 1  ) ) );
    }

    VectorElemental elvec_Iapp ( uFESpace.fe().nbFEDof(), 1 );
    VectorElemental elvec_Iion ( uFESpace.fe().nbFEDof(), 1 );

    for (UInt iVol = 0; iVol < uFESpace.mesh()->numVolumes(); ++iVol)
    {

        uFESpace.fe().updateJacQuadPt ( uFESpace.mesh()->volumeList ( iVol ) );


        for ( int k = 0; k < M_numberOfEquations; k++ )
        {
            ( * ( elvecPtr.at (k) ) ).zero();
        }
        elvec_Iapp.zero();
        elvec_Iion.zero();

        UInt eleIDu = uFESpace.fe().currentLocalId();
        UInt nbNode = ( UInt ) uFESpace.fe().nbFEDof();

        //! Filling local elvec_u with potential values in the nodes
        for ( UInt iNode = 0 ; iNode < nbNode ; iNode++ )
        {

            Int  ig = uFESpace.dof().localToGlobalMap ( eleIDu, iNode );

            for ( int k = 0; k < M_numberOfEquations; k++ )
            {
                ( * ( elvecPtr.at (k) ) ).vec() [iNode] = ( * ( URepPtr.at (k) ) ) [ig];
            }

            elvec_Iapp.vec() [ iNode ] = IappRep[ig];

        }

        //compute the local vector
        for ( UInt ig = 0; ig < uFESpace.fe().nbQuadPt(); ig++ )
        {

            for ( int k = 0; k < M_numberOfEquations; k++ )
            {
                U.at (k) = 0;
            }
            I = 0;

            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {

                for ( int k = 0; k < M_numberOfEquations; k++ )
                {
                    U.at (k) +=  ( * ( elvecPtr.at (k) ) ) (i) *  uFESpace.fe().phi ( i, ig );
                }

                I += elvec_Iapp (i) * uFESpace.fe().phi ( i, ig );

            }

            for ( UInt i = 0; i < uFESpace.fe().nbFEDof(); i++ )
            {

                elvec_Iion ( i ) += ( computeLocalPotentialRhs (U) + I ) * uFESpace.fe().phi ( i, ig ) * uFESpace.fe().weightDet ( ig );

            }

        }

        //assembly
        for ( UInt i = 0 ; i < uFESpace.fe().nbFEDof(); i++ )
        {
            Int  ig = uFESpace.dof().localToGlobalMap ( eleIDu, i );
            ( * ( rhs.at (0) ) ).sumIntoGlobalValues (ig,  elvec_Iion.vec() [i] );
        }
    }
    rhs.at (0) -> globalAssemble();


    if (M_appliedCurrentPtr)
    {
        M_appliedCurrentPtr -> setMapType (Unique);
    }

}

void ElectroIonicModel::computePotentialRhsSVI (   const std::vector<vectorPtr_Type>& v,
                                                   std::vector<vectorPtr_Type>& rhs,
                                                   FESpace<mesh_Type, MapEpetra>& uFESpace,
                                                   const QuadratureRule& qr)
{
    uFESpace.setQuadRule ( qr );
    computePotentialRhsSVI (v, rhs, uFESpace);
}


void ElectroIonicModel::computeGatingVariablesWithRushLarsen ( std::vector<vectorPtr_Type>& v, const Real dt )
{
    int nodes = ( * (v.at (0) ) ).epetraVector().MyLength();

    std::vector<Real>   localVec ( M_numberOfEquations, 0.0 );

    int j (0);

    for ( int k = 0; k < nodes; k++ )
    {

        j = ( * (v.at (0) ) ).blockMap().GID (k);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            localVec.at (i) = ( * ( v.at (i) ) ) [j];
        }

        computeGatingVariablesWithRushLarsen (localVec, dt);

        for ( int i = 0; i < M_numberOfEquations; i++ )
        {
            ( * ( v.at (i) ) ) [j] =  localVec.at (i);
        }

    }

}


void ElectroIonicModel::initialize ( std::vector<Real>& v )
{
    for (int i (0); i <  M_numberOfEquations; i++ )
    {
        v.at (i) = M_restingConditions.at (i);
    }
}


void ElectroIonicModel::initialize ( std::vector<vectorPtr_Type>& v )
{
    for (int i (0); i < M_numberOfEquations; i++ )
    {
        * ( v.at (i) ) = M_restingConditions.at (i);
    }
}


}



