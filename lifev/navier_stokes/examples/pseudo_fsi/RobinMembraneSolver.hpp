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
    @brief

    @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @author Claudia Colciago <claudia.colciago@epfl.ch>
    @date 08-03-2011
 */
#ifndef __RobinMembrane_H
#define __RobinMembrane_H 1

#include "OseenSolverBoundaryDerivative.hpp"
#include <lifev/fsi/solver/HarmonicExtensionSolver.hpp>
#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
#include <lifev/core/fem/BCManageNormal.hpp>
#include "ud_functions.hpp"

#include <lifev/core/fem/DOFInterface3Dto3D.hpp>


enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
 * \class ALE
 * \brief 2D/3D robinMembrane Simulation class
 *
 *  @author Christophe Prud'homme
 *  @see
 */
namespace LifeV {

class RobinMembraneSolver
//     :
//     public LifeV::Application
{
public:


    /** @name Typedefs
     */
    //@{
  
    typedef RegionMesh<LinearTetra>              mesh_Type;
    typedef OseenSolver< RegionMesh<LinearTetra> >::matrix_Type  matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                  matrixPtr_Type;
    typedef OseenSolver< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;
 
    typedef boost::shared_ptr< TimeAdvance< vector_Type > >                timeAdvance_Type;  
    typedef boost::shared_ptr<DOFInterface3Dto3D> dofInterface;
    typedef boost::shared_ptr<MapEpetra> mapEpetra;
    typedef std::map<ID,ID>::const_iterator Iterator;
  typedef OseenSolverBoundaryDerivative< mesh_Type , SolverAztecOO > solver_Type; 

//    typedef LifeV::Application super;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    RobinMembraneSolver( int argc,
	    char** argv );

    ~RobinMembraneSolver()
   
        {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{
  
    //! initialize test
    void initialize();

    //! run test
    void run();
  
    //! transfer vector from he to fluid
    void transferMeshMotionOnFluid( const vector_Type& _vec1, vector_Type& _vec2 );  

    //! interpolate velocity from he to fluid
    void interpolateVelocity( const vector_Type& _vec1, vector_Type& _vec2 ); 

    //! transfer vector from fluid to Interface
    void transferFluidOnInterface(const vector_Type &_vec1, vector_Type &_vec2);

    //! create interface map between from fluid to he
    void createInterfaceMap();

    //! inizialize points for computation of the normals
    void initializeBCNormal(BCManageNormal<matrix_Type>& bcNormal, 
			    const mesh_Type& mesh, 
			    const DOF& dof,
			    CurrentBoundaryFE& currentBdFE,
			    vector_Type& normals);

    //! compute the projection of vel in the normal direction
    void computeVelocityNormalComponent( vector_Type& vel, 
					 const vector_Type& normals, 
					 const MapEpetra& scalarMap,
					 const DOF& dof ); 
    //@}

    //! export global indices
    void exportGID( vector_Type& vectorGID );

    //! import information about gaussian and mean curvature on the surface
    void importCurvatureData( vector_Type& gaussCurve, 
			      vector_Type& meanCurve,
			      GetPot dataFile);
    void
    addRHSRobinGeneral( vector_Type& rightHandSide,
			    const mesh_Type& mesh,
			    const DOF& dof,
			    const BCBase& boundaryCond,
			    CurrentBoundaryFE& currentBdFE,
			    CurrentFE& currentFE,
			    const Real time,
			    UInt offset );


   void 
   initialize(GetPot const& data_file, 
	      vectorPtr_Type velAndPressure, 
	      vectorPtr_Type eta);

private:
    struct Private;
    boost::shared_ptr<Private> M_d;
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > M_uFESpace;
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > M_pFESpace;
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > M_mFESpace;
    boost::shared_ptr<solver_Type> M_fluid;

};

}
#endif /* __RobinMembraneSolver_H */
