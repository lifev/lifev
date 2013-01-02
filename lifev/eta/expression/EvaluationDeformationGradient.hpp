//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains the definition of the EvaluationTranspose class.

     @date 08/2012
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EVALUTATION_DEFORMATIONGRADIENT_HPP
#define EVALUTATION_DEFORMATIONGRADIENT_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorElemental.hpp>

#include <lifev/eta/expression/ExpressionDeformationGradient.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for the transpose of another Evaluation
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing the transpose operation during the assembly

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <typename MeshType, typename MapType, UInt SpaceDim>
class EvaluationDeformationGradient
{
public:

    //! @name Public Types
    //@{

    //! Type of the values returned by this class
    typedef MatrixSmall<SpaceDim, SpaceDim>                                      return_Type;
    typedef VectorEpetra                                                         vector_Type;
    typedef boost::shared_ptr<vector_Type>                                       vectorPtr_Type;
    typedef boost::shared_ptr< ETFESpace<MeshType,MapType,SpaceDim,SpaceDim> >   ETFespacePtr_Type;
    //@}


    //! @name Static constants
    //@{

    //! Flag for the global current FE
	const static flag_Type S_globalUpdateFlag;

    //! Flag for the test current FE
	const static flag_Type S_testUpdateFlag;

    //! Flag for the solution current FE
	const static flag_Type S_solutionUpdateFlag;

    //@}


    //! @name Constructors, destructor
    //@{

	//! Copy constructor
	EvaluationDeformationGradient(const EvaluationDeformationGradient<MeshType,MapType,SpaceDim>& eval)
		:
        M_fespace( eval.M_fespace ),
        M_displacement( eval.M_displacement, Repeated ),
        M_quadrature(0),
        M_currentFE(eval.M_currentFE),
		M_F(eval.M_F)
    {
        if (eval.M_quadrature !=0)
		{
			M_quadrature = new QuadratureRule(*(eval.M_quadrature));
		}
    }

	//! Constructor from the corresponding expression
	explicit EvaluationDeformationGradient(const ExpressionDeformationGradient<MeshType,MapType,SpaceDim>& expression)
		:
        M_fespace( expression.fespace() ),
        M_displacement( expression.vector() ),
        M_offset( expression.offset() ),
        M_currentFE(M_fespace->refFE(),M_fespace->geoMap()),
        M_F( )
    {}

    //! Destructor
    ~EvaluationDeformationGradient()
    {
        if (M_quadrature !=0) delete M_quadrature;
    }

    //@}


    //! @name Methods
    //@{

    //! Computation of the deformation gradient matrix F
	void update(const UInt& iElement)
	{
        //Extracting the local displacement
        VectorElemental dk_loc(this->M_FESpace->refFE().nbFEDof(), SpaceDim);
        UInt totDOFs = this->M_FESpace->dof().numTotlaDof();

	    for ( UInt iNode = 0; iNode < ( UInt ) this->M_FESpace->refFE().nbFEDof(); iNode++ )
        {
            //Local infos on the element
            UInt  iloc = this->M_FESpace->refFE().patternFirst( iNode );
            UInt  eleID = this->M_FESpace->mesh()->volume( iElement ).localid();

            for ( UInt iComp = 0; iComp < SpaceDim; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobalMap( eleID, iloc ) + iComp * totDOFs + this->M_offset;
                dk_loc[iloc + iComp * this->M_FESpace->refFE().nbFEDof()] = M_displacement[ ig ];
            }
        }

        //Computing F
        Real s;

        //! loop on quadrature points (ig)
        for ( UInt ig = 0; ig < this->M_quadrature->nbQuadPt(); ig++ )
        {
            //updating the infos of the derivatives on the quadrature point
            M_currentFE.updateDphi( ig );

            //! loop on space coordinates (icoor)
            for ( UInt icoor = 0; icoor < SpaceDim; icoor++ )
            {
                //! loop  on space coordinates (jcoor)
                for ( UInt jcoor = 0; jcoor < SpaceDim; jcoor++ )
                {
                    s = 0.0;
                    for ( UInt i = 0; i < this->M_FESpace->refFE().nbFEDof(); i++ )
                    {
                        //! \grad u^k at a quadrature point
                        s += this->M_currentFE.dphi( i, jcoor, ig ) * dk_loc[ i + icoor * this->M_FESpace->fe().nbFEDof() ];
                    }
                    //! gradient of displacement
                    M_F[ ig ]( icoor )( jcoor ) = s;
                }
            }
        }

        //! loop on quadrature points (ig)
        for ( UInt ig = 0; ig < this->M_quadrature->nbQuadPt(); ig++ )
        {
            //! loop on space coordinates (icoor)
            for ( UInt  icoor = 0; icoor < SpaceDim; icoor++ )
            {
                //! deformation gradient Fk
                M_F[ ig ]( icoor )( icoor ) +=  1.0;
            }
        }


	}

    //! Show the values
	void showValues() const
	{
		std::cout << " Deformation Gradient Matrix : " << std::endl;

		for (UInt i(0); i<M_quadrature->nbQuadPt(); ++i)
		{
            std::cout << "Deformation Gradient at the node: " << i << std::endl;

            for( UInt j(0); j < 3; j++ )
            {
                for( UInt k(0); k < SpaceDim; k++ )
                {
                    std::cout << "Value ["<< j << " ][ "  << k << " ]: " << M_F[i]( j,k )<< std::endl;
                }
            }
        }
	}

    //! Display method
	static void display(std::ostream& out = std::cout)
	{
        out << " deformation Gradient ";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Setter for the global current FE
	template< typename CFEType >
	void setGlobalCFE(const CFEType* /*globalCFE*/)
	{}

    //! Setter for the test current FE
	template< typename CFEType >
	void setTestCFE(const CFEType* /*testCFE*/)
	{}

    //! Setter for the solution FE
	template< typename CFEType >
	void setSolutionCFE(const CFEType* /*solutionCFE*/)
	{}

    //! Setter for the quadrature rule
	void setQuadrature(const QuadratureRule& qr)
	{
		if (M_quadrature !=0) delete M_quadrature;
		M_quadrature = new QuadratureRule(qr);
		M_currentFE.setQuadratureRule(qr);
		M_F.resize(qr.nbQuadPt());
	}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
	return_Type value_q(const UInt& q) const
    {
        return M_F[ q ];
    }

    //! Getter for the value for a vector
	return_Type value_qi(const UInt& q, const UInt& /*i*/) const
    {
	    return M_F[ q ];
    }

    //! Getter for the value for a matrix
	return_Type value_qij(const UInt& q, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return M_F[ q ];
    }

    //@}

private:

    //! @name Private Methods
    //@{

	//! No default
	EvaluationDeformationGradient();

    //To compute the derivatives
	ETFespacePtr_Type M_fespace;
    //The displacement vector which has to have a repeated map
	vector_Type M_displacement;
    UInt M_offset;

    //Quadrature rule to know where the quadrature points are
	QuadratureRule* M_quadrature;

    //Current ETFE
	ETCurrentFE<SpaceDim,SpaceDim> M_currentFE;

    //DeformationGradient matrix
    std::vector<return_Type> M_F;
    //@}

};


template< typename MeshType, typename MapType, UInt spaceDim >
const flag_Type EvaluationDeformationGradient<MeshType,MapType,spaceDim>::S_globalUpdateFlag
 = ET_UPDATE_NONE;

template< typename MeshType, typename MapType, UInt spaceDim >
const flag_Type EvaluationDeformationGradient<MeshType,MapType,spaceDim>::S_testUpdateFlag
 = ET_UPDATE_NONE;

template< typename MeshType, typename MapType, UInt spaceDim >
const flag_Type EvaluationDeformationGradient<MeshType,MapType,spaceDim>::S_solutionUpdateFlag
 = ET_UPDATE_NONE;

} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
