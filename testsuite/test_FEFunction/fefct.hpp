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
   @file fefct.hpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2011-05-02
 */


#ifndef __fefct_H
#define __fefct_H 1

// ===================================================
//! Includes
// ===================================================

#include <life/lifecore/LifeV.hpp>

#include <life/lifearray/MapEpetra.hpp>

#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>

#include <life/lifefem/FEFunction.hpp>

#ifdef HAVE_HDF5
#include <life/lifefilters/ExporterHDF5.hpp>
#endif
#include <life/lifefilters/ExporterEmpty.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>

/*!
 @class fefct
 @brief LifeV Darcy test case

 @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
 */
class fefct
//     :
//     public LifeV::Application
{
public:

    /*! @name Constructors and destructor
     */
    //@{

    fefct( int argc,
           char** argv );

    ~fefct()
    {}

    //@}

    /*! @name  Methods
     */
    //@{

    LifeV::Real run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> Members;
};

#endif /* __fefct_H */
