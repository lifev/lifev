/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  L. Iapichino  <laura.iapichino@epfl.ch>
              C. Malossi    <cristiano.malossi@epfl.ch>
              A. Manzoni    <andrea.manzoni@epfl.ch>
       Date: 2009-03-24

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
   \file laplacian.hpp
   \author L. Iapichino <laura.iapichino@epfl.ch>, C. Malossi <cristiano.malossi@epfl.ch>, A. Manzoni <andrea.manzoni@epfl.ch>
   \date 2009-03-24
 */


#ifndef __laplacian_H
#define __laplacian_H 1

// ===================================================
//! Includes
// ===================================================
#include <life/lifecore/life.hpp>





/*!
 * \class laplacian
 * \brief LifeV Laplacian test case
 *
 *  @author L. Iapichino, C. Malossi, A. Manzoni
 *  @see
 */
class laplacian
//     :
//     public LifeV::Application
{
public:

    /** @name Constructors, destructor
     */
    //@{

    laplacian( int argc,
               char** argv );

    ~laplacian()
    {}

    //@}

    /** @name  Methods
     */
    //@{

    void run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> Members;
};

#endif /* __laplacian_H */
