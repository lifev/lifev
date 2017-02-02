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
   \file ensightToHdf5.hpp
   \author  <simone.deparis@epfl.ch>
   \date 2008-08-08
 */


#ifndef __EnsightToHdf5_H
#define __EnsightToHdf5_H 1

#include <boost/shared_ptr.hpp>

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

    EnsightToHdf5 ( int argc,
                    char** argv );

    ~EnsightToHdf5()
    {}

    void run();

    //@}


private:
    struct Private;
    std::shared_ptr<Private> d;
};

#endif /* __EnsightToHdf5_H */
