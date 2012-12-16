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

/**
 * @file    test_q2_mesh.cpp
 * @brief   Fair comparison between a Q1 mesh and a Q2 mesh (2D).
 * @author  Simone Pezzuto <simone.pezzuto@mail.polimi.it>
**/

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <boost/shared_ptr.hpp>

#include <lifev/core/LifeV.cpp>
#include <lifev/core/filter/ParserGmsh.hpp>
#include <lifev/core/mesh/ConvertBareMesh.hpp>

#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/core/fem/ReferenceFEScalar.hpp>
#include <lifev/core/fem/CurrentFE.hpp>

#include <cstdlib>
#include <iostream>
#include <cmath>

template <typename S>
struct tester
{
  typedef typename LifeV::BareMesh<S>    baremesh_t;
  typedef typename LifeV::RegionMesh<S>  mesh_t;
  typedef boost::shared_ptr<Epetra_Comm> comm_t;
  typedef LifeV::GeometricMap            geomap_t;
  typedef LifeV::ReferenceFE             reffem_t;

  typedef typename mesh_t::elements_Type   elms_t;
  typedef typename elms_t::const_iterator  elms_citer;

  static LifeV::Real volume(std::string& filename,
                            const reffem_t& reffem,
                            const geomap_t& geomap,
                            comm_t comm)
  {  
    bool ilead = (comm->MyPID() == 0);

    baremesh_t baremesh;

    if (!LifeV::MeshIO::ReadGmshFile(filename, baremesh, 1, ilead))
      return false;

    // Convert the mesh
    boost::shared_ptr<mesh_t> mesh (new mesh_t(comm));
    LifeV::convertBareMesh (baremesh, *mesh, ilead);

    // The current finite element
    LifeV::CurrentFE curFE (reffem, geomap, LifeV::quadRuleQuad4pt);
    // List of all the elements
    elms_t& elms = mesh->elementList();

    LifeV::Real area(0.0);

    // For every element ...
    for (elms_citer it = elms.begin(); it != elms.end(); ++it)
    {
      // Need to update only the weighted determinant
      curFE.update(*it, LifeV::UPDATE_WDET);

      const LifeV::UInt nbQuadPt(curFE.nbQuadPt());
      // Loop on the quadrature nodes
      for (LifeV::UInt iQuadPt(0); iQuadPt < nbQuadPt; ++iQuadPt)
      {
        area += curFE.wDetJacobian(iQuadPt);
      }
    }
    return area;
  }
};

int main (int argc, char* argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  boost::shared_ptr<Epetra_Comm> comm (new Epetra_SerialComm);
#endif

  // Parse command-line
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " [mesh_2D_Q1] [mesh_2D_Q2]" << std::endl;
    return EXIT_SUCCESS;
  }

  // Files
  std::string file1 = argv[1];
  std::string file2 = argv[2];

  bool ilead = (comm->MyPID() == 0);

  // Exact area
  LifeV::Real area_exact = 3. * M_PI;

  if (ilead) std::cout << "\n ====== LINEAR MESH ====== \n\n";

  LifeV::Real v1 = tester<LifeV::LinearQuad>::volume
      (file1, LifeV::feQuadQ1, LifeV::geoBilinearQuad, comm);

  if (ilead) std::cout << "\n :: Total area = "
                       << std::setprecision(10) << v1
                       << "\n :: Exact area = "
                       << area_exact
                       << "\n ::      Error = "
                       << std::setprecision(3)
                       << std::abs(v1-area_exact)/area_exact * 100 << " %"
                       << std::endl;

  if (ilead) std::cout << "\n ====== QUADRATIC MESH ====== \n\n";

  LifeV::Real v2 = tester<LifeV::QuadraticQuad>::volume
      (file2, LifeV::feQuadQ2, LifeV::geoBiquadraticQuad, comm);

  if (ilead) std::cout << "\n :: Total area = "
                       << std::setprecision(10) << v2
                       << "\n :: Exact area = "
                       << area_exact
                       << "\n ::      Error = "
                       << std::setprecision(3)
                       << std::abs(v2-area_exact)/area_exact * 100 << " %\n"
                       << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
