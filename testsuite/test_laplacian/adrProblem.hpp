/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  L. Iapichino  <laura.iapichino@epfl.ch>
              C. Malossi    <cristiano.malossi@epfl.ch>
              A. Manzoni    <andrea.manzoni@epfl.ch>
              T. Passerini  <tiziano@mathcs.emory.edu>

  Copyright (C) 2009-2010 EPFL, Emory University

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
   \file ADRProblem.hpp
   \author L. Iapichino <laura.iapichino@epfl.ch>
           C. Malossi <cristiano.malossi@epfl.ch>
           A. Manzoni <andrea.manzoni@epfl.ch>
           T. Passerini <tiziano@mathcs.emory.edu>
 */





#ifndef __ADRProblem_H
#define __ADRProblem_H 1





// ===================================================
//! Includes
// ===================================================
#include <life/lifecore/life.hpp>





/*!
 * \class ADRProblem
 * \brief LifeV advection-diffusion-reaction test case
 *
 *  @author L. Iapichino, C. Malossi, A. Manzoni, T. Passerini
 *
 */
class ADRProblem
{
public:

    /** @name Constructors, destructor
     */
    //@{

    ADRProblem( int argc,
                char** argv );

    ~ADRProblem()
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

#endif /* __ADRProblem_H */
