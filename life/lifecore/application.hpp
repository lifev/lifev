/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-03-17

  Copyright (C) 2005 EPFL

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
   \file application.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-03-17
 */
#ifndef __application_H
#define __application_H 1

#include <boost/program_options.hpp>
#include <life/lifecore/about.hpp>

namespace LifeV
{
namespace po = boost::program_options;

/**
 * \class Application
 * \brief provides information about the Application
 *
 * @author Christophe Prud'homme
 */
class Application
{
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     *
     *
     * @param argc
     * @param argv
     * @param od
     */
    Application( int argc, char** argv, AboutData const& ad );

    /**
     *
     *
     * @param argc
     * @param argv
     * @param od
     */
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    /**
     *
     *
     */
    Application( Application const & );

    ~Application();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * get the options description
     *
     *
     * @return the options description
     */
    po::options_description const& optionsDescription() const { return _M_desc; }

    /**
     * get the variable map
     *
     *
     * @return the variable map
     */
    po::variables_map const& vm() const { return _M_vm; }

    /**
     * get the about data of the application
     *
     *
     * @return the about data ofthe application
     * @see AboutData
     */
    AboutData const& about() const { return _M_about; }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{



    //@}



protected:

private:

    /**
     * process the generic options passed to the command line
     *
     */
    void processGenericOptions() const;

private:

    AboutData _M_about;

    po::options_description _M_desc;
    po::variables_map _M_vm;

};

}
#endif /* __Application_H */
