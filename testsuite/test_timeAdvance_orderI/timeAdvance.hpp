/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  F. Nobile  <fabio.nobile@polimi.it>
              M. Pozzoli    <matteo1.pozzoli@mail.polimi.it>
              C. Vergara    <christian.vergara@polimi.it>
       Date: 2010-02-20

  Copyright (C) 2009 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file timeAdvance.hpp
   \author F. Nobile, M. Pozzoli, C. Vergara
   \date 2010-02-15
 */





#ifndef __timeAdvance_H
#define __timeAdvance_H





// ===================================================
//! Includes
// ===================================================
#include <life/lifecore/life.hpp>

/*!
 * \class problem
 * \brief LifeV problem test case
 *
 *  @author F. Nobile, M. Pozzoli, C. Vergara
 *  @see
 */
class problem
//
{
public:

    /** @name Constructors, destructor
     */
    //@{
    problem( int          argc,
             char**                argv,
             boost::shared_ptr<Epetra_Comm>        structComm );

    ~problem()
    {}

    //@}

    /** @name  Methods
     */
    //@{

    void run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> members;
};

#endif /* __timeAdvance_H */
