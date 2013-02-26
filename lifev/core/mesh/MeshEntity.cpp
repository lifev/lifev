/*
 * MeshEntity.cpp
 *
 *  Created on: 11 Sep 2011
 *      Author: formaggia
 */
#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/MeshEntity.hpp>

namespace LifeV
{

namespace EntityFlags
{

std::string name ( const flag_Type& flag )
{
    switch ( flag )
    {
        case PHYSICAL_BOUNDARY:
            return "PHYSICAL_BOUNDARY";
        case INTERNAL_INTERFACE:
            return "INTERNAL_INTERFACE";
        case SUBDOMAIN_INTERFACE:
            return "SUBDOMAIN_INTERFACE";
        case OVERLAP:
            return "OVERLAP";
        case CUTTED:
            return "CUTTED";
        case VERTEX:
            return "VERTEX";
        case GHOST_ENTITY:
            return "GHOST_ENTITY";
        default:
            return "FLAG_NOT_FOUND";
    }
}

}// namespace EntityFlags

void MeshEntity::showMe ( std::ostream& output) const
{
    output << " Global ID : " << M_id << " -- " << " Local ID " << M_localId << std::endl;
    Flag::showMe ( M_flag, output );
    output << std::endl;
}

}// namespace LifeV


