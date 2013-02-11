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
   \file timeAdvance.hpp
   \date 2010-02-15
  Author(s):   F. Nobile  <fabio.nobile@polimi.it>
              M. Pozzoli    <matteo1.pozzoli@mail.polimi.it>
              C. Vergara    <christian.vergara@polimi.it>
 */





#ifndef __timeAdvance_H
#define __timeAdvance_H





// ===================================================
//! Includes
// ===================================================
#include <lifev/core/LifeV.hpp>


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

    problem ( int argc,
              char** argv,
              boost::shared_ptr<Epetra_Comm>        structComm);

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
