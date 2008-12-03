/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s):  <simone.deparis@epfl.ch>
       Date: 2008-08-08

  Copyright (C) 2008 EPFL

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
   \file ensightToHdf5.hpp
   \author  <simone.deparis@epfl.ch>
   \date 2008-08-08
 */


#ifndef __EnsightToHdf5_H
#define __EnsightToHdf5_H 1

#include <life/lifecore/application.hpp>

class EnsightToHdf5

{
public:


    /** @name Typedefs
     */
    //@{

//    typedef LifeV::Application super;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    EnsightToHdf5( int argc,
                   char** argv,
                   LifeV::AboutData const& ad,
                   LifeV::po::options_description const& od );

    ~EnsightToHdf5()
        {}

    void run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> d;
};

#endif /* __EnsightToHdf5_H */
