/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  A. Fumagalli <alessio.fumagalli@mail.polimi.it>
       Date: 2010-05-24

  Copyright (C) 2010 EPFL, Politecnico di Milano

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
   @file darcy.hpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2009-03-24
 */


#ifndef __darcy_H
#define __darcy_H 1


// ===================================================
//! Includes
// ===================================================

#include <life/lifecore/application.hpp>
#include <life/lifecore/life.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
	#include <Epetra_MpiComm.h>
	#include <mpi.h>
#else
	#include <Epetra_SerialComm.h>
#endif

#include <life/lifemesh/structuredMesh3D.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/darcySolver.hpp>



/*!
 @class darcy
 @brief LifeV Darcy test case

 @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
 */
class darcy
//     :
//     public LifeV::Application
{
public:

    /*! @name Constructors and destructor
     */
    //@{

    darcy( int argc,
           char** argv,
           LifeV::AboutData const& ad,
           LifeV::po::options_description const& od );

    ~darcy()
        {}

    //@}

    /*! @name  Methods
     */
    //@{

    void run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> Members;
};

#endif /* __darcy_H */
