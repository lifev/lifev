/*
  This file is part of the LifeV library.

  Author: Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  Copyright (C) 2004 EPFL

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
#include <debug.hpp>
int
main()
{

    LifeV::Debug() << "Hello\n";


    //
    // To see this message setup the DEBUG variable to 2000
    // export DEBUG=2000
    // then execute test_debug
    //
    LifeV::Debug(2000) << "AREA 2000 is now enabled\n";

}

