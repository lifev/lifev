/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-10-08

  Copyright (C) 2010 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file BlockPreconditioner.hpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-10-08
 */


#ifndef _BLOCKPRECONDITIONER_HPP_
#define _BLOCKPRECONDITIONER_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/EpetraMap.hpp>

namespace LifeV {

//! BlockPreconditioner
/*!
 *  @author Gwenol Grandperrin
 *
 *  The BlockPreconditioner class is an abstract class which
 *  defines the interfaces for a typical block preconditioner
 */
class BlockPreconditioner:
        public EpetraPreconditioner
{
public:

    /** @name Typedefs
     */
    //@{
    typedef EpetraPreconditioner                 super;

    typedef Epetra_Operator                      prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;

    typedef super::operator_raw_type             operator_raw_type;
    typedef super::operator_type                 operator_type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    BlockPreconditioner(const boost::shared_ptr<Epetra_Comm>& comm = boost::shared_ptr<Epetra_Comm>());

    /** Copy constructor*/
    BlockPreconditioner(BlockPreconditioner& P, const boost::shared_ptr<Epetra_Comm>& comm = boost::shared_ptr<Epetra_Comm>() );

    //! default virtual destructor
    virtual ~BlockPreconditioner();

    //@}


    /** @name  Methods
     */
    virtual int numBlockRow() const=0;
    virtual int numBlockCol() const=0;

protected:

};

void buildBlockGIDs(std::vector<std::vector<int> > & gids,const EpetraMap & map,const std::vector<int>& blockSizes);

} // namespace LifeV

#endif
