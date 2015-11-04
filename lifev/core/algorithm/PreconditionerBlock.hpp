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
    @file PreconditionerBlock.cpp
    @brief PreconditionerBlock

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 08-10-2010
 */

#ifndef PRECONDITIONERBLOCK_HPP
#define PRECONDITIONERBLOCK_HPP 1

#include <vector>
#include <boost/shared_ptr.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/array/MapEpetra.hpp>

namespace LifeV
{

//! PreconditionerBlock
/*!
 *  @author Gwenol Grandperrin
 *
 *  The PreconditionerBlock class is an abstract class which
 *  defines the interfaces for a typical block preconditioner
 */
class PreconditionerBlock:
    public Preconditioner
{
public:

    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    PreconditionerBlock ( const std::shared_ptr<Epetra_Comm>& comm = std::shared_ptr<Epetra_Comm>() );

    /** Copy constructor*/
    PreconditionerBlock ( PreconditionerBlock& P, const std::shared_ptr<Epetra_Comm>& comm = std::shared_ptr<Epetra_Comm>() );

    //! default virtual destructor
    virtual ~PreconditionerBlock();

    //@}


    /** @name  Methods
     */
    virtual int numBlocksRows() const = 0;
    virtual int numBlocksCols() const = 0;

protected:

};

void buildBlockGIDs ( std::vector<std::vector<int> >& gids, const MapEpetra& map, const std::vector<int>& blockSizes );

} // namespace LifeV

#endif /* PRECONDITIONERBLOCK_HPP */
