/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file boundaryConditions.hpp
*/

#ifndef BC_HPP
#define BC_HPP

#include "life/lifecore/life.hpp"
#include "ud_functions.hpp"
#include "life/lifefem/bcHandler.hpp"
#include "life/lifefem/bcFunction.hpp"

namespace LifeV
{

typedef FSISolver::bchandler_type bchandler_type;


bchandler_type BCh_mesh()
{
    // Boundary condition for the mesh
    Debug( 10000 ) << "Boundary condition for the mesh\n";

    BCFunctionBase bcf(fZero);
    FSISolver::bchandler_type BCh_mesh( new BCHandler );
    BCh_mesh->addBC("Top",       3, Essential, Full, bcf,   3);
    BCh_mesh->addBC("Base",      2, Essential, Full, bcf,   3);
    BCh_mesh->addBC("Edges",    20, Essential, Full, bcf,   3);

    return BCh_mesh;
}

bchandler_type BCh_u()
{
    // Boundary conditions for the fluid velocity
    Debug( 10000 ) << "Boundary condition for the fluid\n";
    BCFunctionBase in_flow(u2);
    FSISolver::bchandler_type BCh_u( new BCHandler );
    BCFunctionBase bcf(fZero);

    BCh_u->addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    BCh_u->addBC("Edges",  20, Essential, Full, bcf,     3);

    return BCh_u;
}

bchandler_type BCh_d()
{
    // Boundary conditions for the solid displacement
    Debug( 10000 ) << "Boundary condition for the solid\n";
    FSISolver::bchandler_type BCh_d( new BCHandler );
    BCFunctionBase bcf(fZero);

    BCh_d->addBC("Top",       3, Essential, Full, bcf,  3);
    BCh_d->addBC("Base",      2, Essential, Full, bcf,  3);

    return BCh_d;
}

}

#endif
