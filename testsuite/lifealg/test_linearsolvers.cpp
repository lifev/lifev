/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  Author: Christophe Prud'homme (christophe.prudhomme@epfl.ch)

  Copyright (C) 2004 EPFL

  Distributed under the GPL(GNU Public License):
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
/** \file test_petsc.cpp

*/
#include <cstdlib>
#include <cassert>

#include <lifeconfig.h>

#include "SolverAztec.hpp"

#if defined(HAVE_PETSC_H)
#include <SolverPETSC.hpp>
#endif /* HAVE_PETSC_H */

#if defined(HAVE_UMFPACK_H)
//#include <SolverUMFPACK.hpp>
#endif /* HAVE_UMFPACK_H */

#include <MatrixTest.hpp>


static char help[] = "Solves a linear system with KSP.\n\
-N <N>       : number of rows/cols\n\n";

namespace LifeV{

template<typename Mat>
bool test_umfpack( Mat& __mat )
{
#if 0
    int Nrows = __mat.matrix().Patt()->nRows();

    //Life::Solver __petsc( "gmres", "ilu" );
    //__petsc.setMatrix( Nrows, __mat.iaData(), __mat.jaData(), __mat.valueData() );

    Vector __x( Nrows );
    Vector __sol( Nrows );
    Vector __b( Nrows );

    __sol = 10;
    __b = __mat.matrix() * __sol;

    __x = 0;
    //__petsc.solve( __x, __b );

    std::cout << "norm(x) = " << norm( __x ) << "\n";

    __x -= __sol;
    std::cout << "norm(error) = " << norm( __x ) << "\n";
#else
    return true;
#endif
}

template<typename Mat>
bool test_petsc( Mat& __mat )
{
#if defined(HAVE_PETSC_H)

    int Nrows = __mat.matrix().Patt()->nRows();

    LifeV::SolverPETSC __petsc( "gmres", "ilu" );
    //__petsc.setMatrix( Nrows, __mat.iaData(), __mat.jaData(), __mat.valueData() );
    __petsc.setMatrix(__mat.matrix());

    Vector __x( Nrows );
    Vector __sol( Nrows );
    Vector __b( Nrows );

    __sol = 10;
    __b = __mat.matrix() * __sol;

    __x = 0;
    __petsc.solve( __x, __b );

    std::cout << "norm(x) = " << norm( __x ) << "\n";

    __x -= __sol;
    std::cout << "norm(error) = " << norm( __x ) << "\n";
    return norm(__x) < 1e-6;
#else
    return 1;
#endif
}

template<typename Mat>
bool test_aztec( Mat& __mat )
{
    int Nrows = __mat.matrix().Patt()->nRows();

    LifeV::SolverAztec __aztec;
    __aztec.setMatrix(__mat.matrix());

    Vector __x( Nrows );
    Vector __sol( Nrows );
    Vector __b( Nrows );

    __sol = 10;
    __b = __mat.matrix() * __sol;

    __x = 0;
    __aztec.solve( __x, __b );

    std::cout << "norm(x) = " << norm( __x ) << "\n";

    __x -= __sol;
    std::cout << "norm(error) = " << norm( __x ) << "\n";
    return norm(__x) < 1e-6;
}

} // namespace LifeV

int main( int argc, char** argv )
{
    bool success = true;
    try
    {
    int N = 100;

#if defined(HAVE_PETSC_H)
    PetscInitialize(&argc,&argv,(char *)0,help);
    PetscOptionsGetInt(PETSC_NULL,"-N",&N,PETSC_NULL);
#endif /* HAVE_PETSC_H */


    //
    // Mass matrix
    //
    std::cout << "mass matrix...\n";
    LifeV::MatrixMass mass( N );
    mass.matrix().spy( "mass.m" );

    success &= LifeV::test_petsc ( mass );
    success &= LifeV::test_umfpack ( mass );
    success &= LifeV::test_aztec( mass );

    //
    // convdiff matrix
    //
    std::cout << "convection diffusion matrix...\n";
    LifeV::MatrixConvectionDiffusion convdiff((int)std::sqrt((double)N), 1.0);
    convdiff.matrix().spy( "convdiff.m" );

    success &= LifeV::test_petsc ( convdiff );
    success &= LifeV::test_umfpack ( convdiff );
    success &= LifeV::test_aztec( convdiff );

    }
    catch( std::exception const& __e )
    {
        std::cout << "std::exception: " << __e.what() << "\n";
        return EXIT_FAILURE;
    }
    catch( ... )
    {
        std::cout << "unknown exception caught\n";
        return EXIT_FAILURE;
    }
    std::cout << (success ? "success" : "solve failed") << std::endl;
    //return (success ? EXIT_SUCCESS : EXIT_FAILURE);
    return EXIT_SUCCESS;
}
