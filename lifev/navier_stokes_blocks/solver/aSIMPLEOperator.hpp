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
    @brief SIMPLE preconditioner for Navier-Stokes equations.

    @author Davide Forti <davide.forti@epfl.ch>
    @contributor Umberto Villa
    @date 08-12-2014

    @maintainer Davide Forti <davide.forti@epfl.ch>
 */

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/navier_stokes_blocks/solver/NavierStokesPreconditionerOperator.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>



#ifndef _aSIMPLEOPERATOR_H_
#define _aSIMPLEOPERATOR_H_ 1

namespace LifeV
{
namespace Operators
{

class aSIMPLEOperator: public NavierStokesPreconditionerOperator
{
public:
    //! @name Public Types
    //@{

    typedef  Epetra_MultiVector                        vector_Type;
    typedef  std::shared_ptr<vector_Type>            vectorPtr_Type;
    typedef  Epetra_Map                                map_Type;
    typedef  std::shared_ptr<map_Type> 			   mapPtr_Type;
    typedef  LinearOperatorAlgebra                     super;
    typedef  Epetra_CrsMatrix                          matrix_Type;
    typedef  std::shared_ptr<matrix_Type>            matrixPtr_Type;
    typedef  MatrixEpetra<Real>                        matrixEpetra_Type;
    typedef  std::shared_ptr<matrixEpetra_Type>      matrixEpetraPtr_Type;
    typedef  Epetra_Vector                             lumpedMatrix_Type;
    typedef  std::shared_ptr<lumpedMatrix_Type>      lumpedMatrixPtr_Type;
    typedef  super::comm_Type                          comm_Type;
    typedef  super::commPtr_Type                       commPtr_Type;
    typedef  std::shared_ptr<Teuchos::ParameterList> parameterListPtr_Type;
    typedef  MapEpetra                                 mapEpetra_Type;
    typedef  std::shared_ptr<mapEpetra_Type>         mapEpetraPtr_Type;
    typedef  VectorEpetra                              VectorEpetra_Type;
    typedef  std::shared_ptr<VectorEpetra_Type>      VectorEpetraPtr_Type;
    //@}

    //! @name Constructors and Destructors
    //@{

    //! Empty constructor
    aSIMPLEOperator();

    //! Destructor
    virtual ~aSIMPLEOperator();

    //@}

    //! @name SetUp
    //@{

    //! SetUp - case without stabilization
    /*!
	 *  @param F block(0,0) of NS matrix
	 *  @param B block(1,0) of NS matrix
	 *  @param Btranspose block(0,1) of NS matrix
     */
    void setUp(const matrixEpetraPtr_Type & F,
               const matrixEpetraPtr_Type & B,
               const matrixEpetraPtr_Type & Btranspose);

    //! SetUp - case with stabilization
    /*!
     *  @param F block(0,0) of NS matrix
     *  @param B block(1,0) of NS matrix
     *  @param Btranspose block(0,1) of NS matrix
     *  @param D block(1,1) of NS matrix
     */
    void setUp ( const matrixEpetraPtr_Type & F,
   				 const matrixEpetraPtr_Type & B,
				 const matrixEpetraPtr_Type & Btranspose,
				 const matrixEpetraPtr_Type & D );


    //! @name Setters
    //@{

    //! \warning Transpose of this operator is not supported
    int SetUseTranspose(bool UseTranspose){M_useTranspose = UseTranspose; return 0;}

    //! set the domain map
    void setDomainMap(const std::shared_ptr<BlockEpetra_Map> & domainMap){M_operatorDomainMap = domainMap;}

    //! set the range map
    void setRangeMap(const std::shared_ptr<BlockEpetra_Map> & rangeMap){M_operatorRangeMap = rangeMap;}

    //! Set the momentum preconditioner options
    void setMomentumOptions(const parameterListPtr_Type & _oList);

    //! Set the Schur Complement preconditioner options
    void setSchurOptions(const parameterListPtr_Type & _oList);

    //@}


    //! @name Methods
    //@{

    //! \warning No method \c Apply defined for this operator. It return an error code.
    int Apply(const vector_Type &/*X*/, vector_Type &/*Y*/) const {return -1;};

    //! Returns the High Order Yosida approximation of the inverse pressure Schur Complement applied to \c (Xu, Xp).
    int ApplyInverse( VectorEpetra_Type const& X_velocity,
                      VectorEpetra_Type const& X_pressure,
                      VectorEpetra_Type & Y_velocity,
                      VectorEpetra_Type & Y_pressure) const;

    //! Returns the High Order Yosida approximation of the inverse pressure Schur Complement applied to \c X.
    int ApplyInverse(const vector_Type &X, vector_Type &Y) const;

    //! \warning Infinity norm not defined for this operator
    double NormInf() const {return -1.0;}

    //! Updates the momentum preconditioner operator
    void updateApproximatedMomentumOperator();

    //! Updates the Schur Complement preconditioner operator
    void updateApproximatedSchurComplementOperator();

    //@}

    // @name Attribute access functions
    //@{
    //! Return a character string describing the operator
    const char * Label() const {return M_label.c_str();}
    //! Return the current UseTranspose setting \warning Not Supported Yet.
    bool UseTranspose() const {return M_useTranspose;}
    //! Return false.
    bool HasNormInf() const {return false;}
    //! return a reference to the Epetra_Comm communicator associated with this operator
    const comm_Type & Comm() const {return *M_comm;}
    //! Returns the Epetra_Map object associated with the domain of this operator
    const map_Type & OperatorDomainMap() const {return *(M_operatorDomainMap->monolithicMap());}
    //! Returns the Epetra_Map object associated with the range of this operator
    const map_Type & OperatorRangeMap() const {return *(M_operatorRangeMap->monolithicMap());}
    //@}

    //! @name Getters
    //@{

    //! Show information about the class
    void showMe();

    //! Return the list of options being used
    void setOptions(const Teuchos::ParameterList& solversOptions);

    //! Return the block(0,0)
    matrixEpetraPtr_Type const& F() const { return M_F; }

    //! Return the block(0,1)
    matrixEpetraPtr_Type const& B() const { return M_B; }

    //! Return the block(1,0)
    matrixEpetraPtr_Type const& Btranspose() const { return M_Btranspose; }

    //@}

private:

    //! Create the domain and the range maps
    void setMaps();

    //! create the matrix B*diag(F)^-1*Btranspose
    void buildShurComplement();

    std::shared_ptr<BlockEpetra_Map> M_operatorDomainMap;
    //! Range Map
    std::shared_ptr<BlockEpetra_Map> M_operatorRangeMap;

    matrixEpetraPtr_Type M_F;

    matrixEpetraPtr_Type M_B;

    matrixEpetraPtr_Type M_Btranspose;

    matrixEpetraPtr_Type M_D;

    matrixEpetraPtr_Type M_schurComplement;

    //! Communicator
    commPtr_Type M_comm;

    bool M_useTranspose;

    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedMomentumOperator;

    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedSchurComplementOperator;

    parameterListPtr_Type M_momentumOptions;

    parameterListPtr_Type M_schurOptions;

    mapEpetraPtr_Type M_monolithicMap;

    std::shared_ptr<Epetra_Vector> M_invD;

    matrixEpetraPtr_Type M_DBT;

    //! Label
    const std::string M_label;

    //! Vectors needed for the apply inverse
    std::shared_ptr<VectorEpetra_Type> M_Z;

    std::shared_ptr<VectorEpetra_Type> M_X_velocity;
    std::shared_ptr<VectorEpetra_Type> M_X_pressure;
    std::shared_ptr<VectorEpetra_Type> M_Y_velocity;
    std::shared_ptr<VectorEpetra_Type> M_Y_pressure;
    
    std::shared_ptr<mapEpetra_Type> M_domainDBT;
    std::shared_ptr<mapEpetra_Type> M_rangeDBT;

    bool M_useStabilization;

};

//! Factory create function
inline NavierStokesPreconditionerOperator * create_aSIMPLE()
{
    return new aSIMPLEOperator ();
}
namespace
{
static bool S_register_aSimple = NSPreconditionerFactory::instance().registerProduct ( "SIMPLE", &create_aSIMPLE );
}

} /* end namespace Operators */
} //end namespace
#endif
