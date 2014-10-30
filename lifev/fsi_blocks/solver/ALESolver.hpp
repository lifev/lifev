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
    @brief Classes to hold algorithms for the mesh motion, for instance, involved in a ALE formulation.
    @author G. Fourestey
    @date 01-10-2007

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>

    This file contains classes which may be used to compute the extension inside the reference domain of a given
    displacement at a specified interface

*/

#ifndef ALESOLVER_H_
#define ALESOLVER_H_

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/fem/DOF.hpp>
#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/util/Displayer.hpp>

namespace LifeV
{
/*!
  \class ALESolver

  Base class which provides the harmonic extension of a given displacement on a specified part
  of the mesh boundary

  In order to deal with harmonic extensions, we have to provide a mesh (to  be moved), the parameters
  involved in the laplacian discretization: viscosity, quadrature rules and boundary conditions.
  This class contains a   PhysVectUnknown objet wich will hold the extension of the interface displacement.
  The constructor of the class built the global matrix of the discretized laplacian. The extension of the
  displacement is computed by calling the public method update. Finally, this extension can be recovered by
  calling method getDisplacement.

*/

class ALESolver
{
public:

    //! @name Public Types
    //@{

	typedef RegionMesh<LinearTetra> mesh_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructors for an harmonics extensions
    /*!
      \param mmFESpace the FEspace that describes the problem
      \param comm  the Epetra_Comm to be used for communication
    */

    ALESolver ( FESpace<mesh_Type, MapEpetra>& mmFESpace, boost::shared_ptr<Epetra_Comm>  comm);

    //! Constructors for an harmonics extensions with offset
    /*!
      \param mmFESpace the FEspace that describes the problem
      \param comm  the Epetra_Comm to be used for communication
      \param localMap use localMap instead of M_FESpace.map()
      \param offset use this offset to fill the matrix (both: row and column offset)
    */
    ALESolver ( FESpace<mesh_Type, MapEpetra>&      mmFESpace,
                              boost::shared_ptr<Epetra_Comm> comm,
                              MapEpetra&                     localMap,
                              UInt                           offset = 0
                            );

    //! virtual destructor
    virtual ~ALESolver() {};

    //@}

    //! @name Methods
    //@{

    //! Set up data from GetPot
    /*!
        @param dataFile GetPot object
     */
    void setUp ( const GetPot& dataFile );

    //! Update convective term, boundary condition and solve the linearized ns system
    /*!
        @param bcHandler BC handler
     */
    void iterate (BCHandler& BCh);

    //! returns wheter this processor is the leader.
    bool isLeader() const
    {
        return M_displayer.comm()->MyPID() == 0;
    }

    //! manually rescale the system matrix by dt
    void rescaleMatrix (Real& dt)
    {
        *M_matrHE *= dt;
    }

    matrixPtr_Type const& matrix() const
    {
    	return M_matrHE;
    }

    //! Adds the system matrix to the argument
    void addSystemMatrixTo (matrixPtr_Type matr) const
    {
        *matr += *M_matrHE;
    }

    //! Apply boundary conditions.
    /*!
        @param rightHandSide
        @param bcHandler
     */
    void applyBoundaryConditions (vector_Type& rhs, BCHandler& BCh);

    void computeMatrix();
    void updateDispDiff();

    //@}


    //! @name Get Methods
    //@{
    vector_Type const& disp()     const
    {
        return *M_disp;
    }
    vector_Type& disp()
    {
        return *M_disp;
    }

    MapEpetra const& getMap() const
    {
        return M_localMap;
    }

    FESpace<mesh_Type, MapEpetra> const& mFESpace() const
    {
        return M_FESpace;
    }

    const boost::shared_ptr<Epetra_Comm>& comm() const
    {
        return M_displayer.comm();
    }
    //@}

private:

    //! Finite Element Space

    FESpace<mesh_Type, MapEpetra>&      M_FESpace;

    //! local map
    MapEpetra                      M_localMap;

    //! The matrix holding the values
    matrixPtr_Type                 M_matrHE;

    Displayer                      M_displayer;
    int                            M_me;
    bool                           M_verbose;

    //! Elementary matrix : 3 blocks
    MatrixElemental                        M_elmat;

    //! The actual extension of the displacement
    vectorPtr_Type                    M_disp;

    //! Auxiliary vector holding the second right hand of the system
    vectorPtr_Type                    M_secondRHS;

    //! Diffusion coefficient for the laplacian operator
    Real                           M_diffusion;

    UInt                           M_offset;
};

} // namespace LifeV
#endif //  ALESOLVER_H_
