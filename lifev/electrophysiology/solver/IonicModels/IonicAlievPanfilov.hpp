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
  @file IonicAlievPanfilov
  @brief Ionic model of Aliev-Panfilov

  Model from:
  Nash, Martyn P., and Alexander V. Panfilov.
  "Electromechanical model of excitable tissue to study reentrant cardiac arrhythmias."
  Progress in biophysics and molecular biology 85.2 (2004): 501-522.

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _IONICALIEVPANFILOV_H_
#define _IONICALIEVPANFILOV_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>




namespace LifeV
{
//! IonicModel - This class implements an ionic model.


class IonicAlievPanfilov : public virtual ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef ElectroIonicModel super;
    typedef std::shared_ptr<VectorEpetra> vectorPtr_Type;
    typedef RegionMesh<LinearTetra> mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicAlievPanfilov();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicAlievPanfilov ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicAlievPanfilov object
     */
    IonicAlievPanfilov ( const IonicAlievPanfilov& model );
    //! Destructor
    virtual ~IonicAlievPanfilov() {}

    //@}

    //! @name Overloads
    //@{

    IonicAlievPanfilov& operator= ( const IonicAlievPanfilov& model );

    //@}

    //! @name Setters and getters
    //@{

    //parameters getters and setters
    inline const Real& Mu1()        const
    {
        return M_mu1;
    }
    inline const Real& Mu2()        const
    {
        return M_mu2;
    }
    inline const Real& K()      const
    {
        return M_k;
    }
    inline const Real& A()      const
    {
        return M_a;
    }
    inline const Real& Epsilon()    const
    {
        return M_epsilon;
    }

    inline void setMu1      ( const Real& mu1 )
    {
        this->M_mu1       = mu1;
    }
    inline void setMu2      ( const Real& mu2 )
    {
        this->M_mu2       = mu2;
    }
    inline void setK            ( const Real& k )
    {
        this->M_k         = k;
    }
    inline void setA            ( const Real& a )
    {
        this->M_a         = a;
    }
    inline void setEpsilon  ( const Real& epsilon )
    {
        this->M_epsilon   = epsilon;
    }


    void setup ( Teuchos::ParameterList& parameterList );
    //@}

    //! @name Methods
    //@{

    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    //Compute the rhs on a mesh/ 3D case
    //    void computeRhs( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );
    //
    //    void computeRhs( const std::vector<vectorPtr_Type>& v, const VectorEpetra& Iapp, std::vector<vectorPtr_Type>& rhs );
    //

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v);

    //    void computePotentialRhs(     const std::vector<vectorPtr_Type>& v,
    //                      const VectorEpetra& Iapp,
    //                      std::vector<vectorPtr_Type>& rhs,
    //                      FESpace<mesh_Type, MapEpetra>& uFESpace );

    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_mu1;
    Real M_mu2;
    Real M_epsilon;
    Real M_k;
    Real M_a;


    //@}

}; // class IonicAlievPanfilov



inline ElectroIonicModel* createIonicAlievPanfilov()
{
    return new IonicAlievPanfilov();
}

namespace
{
static bool register_IonicAlievPanfilov = ElectroIonicModel::IonicModelFactory::instance().registerProduct ("AlievPanfilov", &createIonicAlievPanfilov );
}

}

#endif
