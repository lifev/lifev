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
#ifndef ET_ROBIN_MEMBRANE_SOLVER
#define ET_ROBIN_MEMBRANE_SOLVER

#include <lifev/core/algorithm/SolverAztecOO.hpp>

//#include <lifev/core/array/MatrixBlockMonolithicEpetra.hpp>
//#include <lifev/core/array/VectorBlockMonolithicEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <lifev/core/util/LifeChrono.hpp>

//#include "OseenSolverBoundaryDerivative.hpp"
#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
//#include <lifev/core/mesh/RegionMesh3D.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/fem/DOFInterface3Dto3D.hpp>
#include <lifev/navier_stokes/solver/OseenAssembler.hpp>


enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
 * \class ALE
 * \brief 2D/3D robinMembrane Simulation class
 *
 *  @author Claudia Colciago
 *  @see
 */
namespace LifeV
{

class ETRobinMembraneSolver
    //     :
    //     public LifeV::Application
{
public:


    /** @name Typedefs
     */
    //@{

    typedef RegionMesh<LinearTetra>                                mesh_type;
    //typedef MatrixBlockMonolithicEpetra<Real>                      matrix_block_type;
    //typedef VectorBlockMonolithicEpetra                            vector_block_type;
    typedef MatrixEpetraStructured<Real>                      matrix_block_type;
    typedef VectorEpetraStructured                          vector_block_type;
    typedef MatrixEpetra<Real>                                     matrix_type;
    typedef VectorEpetra                                           vector_type;
    typedef std::shared_ptr<matrix_type>                         matrixPtr_type;
    typedef std::shared_ptr<mesh_type>                           meshPtr_type;
    typedef std::shared_ptr<vector_type>                         vectorPtr_type;

    typedef std::shared_ptr< TimeAdvance< vector_type > >        timeAdvancePtr_type;
    typedef OseenAssembler< mesh_type , matrix_type, vector_type > assembler_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{
    ETRobinMembraneSolver ( int argc,
                            char** argv );

    ~ETRobinMembraneSolver()

    {}

    //@}

    /** @name  Methods
     */
    //@{

    //! initialize test
    void initialize();

    //! create interface map
    //void createInterfaceMap ( vectorPtr_type checkVector, meshPtr_type& mesh , const DOF& dof);
    void createInterfaceMap ( std::map<ID, ID> const& locDofMap , const DOF& dof );

    //! run test
    void run();

private:

    struct Private;
    std::shared_ptr<Private>                                    M_d;
    std::shared_ptr<FESpace<mesh_type, MapEpetra> >             M_uFESpace;
    std::shared_ptr<FESpace<mesh_type, MapEpetra> >             M_uCompFESpace;
    std::shared_ptr<FESpace<mesh_type, MapEpetra> >             M_pFESpace;
    std::shared_ptr<FESpace<mesh_type, MapEpetra> >             M_mFESpace;
    std::shared_ptr<ETFESpace<mesh_type, MapEpetra , 3 , 3 > >  M_ETuFESpace;
    std::shared_ptr<ETFESpace<mesh_type, MapEpetra , 3 , 1 > >  M_ETpFESpace;
    std::shared_ptr<assembler_type>                             M_pseudoFSI;
    std::shared_ptr<MapEpetra>                                  M_interfaceMap;

};

}
#endif /* ET_ROBIN_MEMBRANE_SOLVER */
