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
   \file application.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-03-17
 */
#include <iostream>
#include <iomanip>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

namespace LifeV
{
Application::Application( int argc, char** argv, AboutData const& ad )
    :
    _M_about( ad ),
    _M_desc( "Allowed options" ),
    _M_vm()
{
    po::options_description generic( "Generic options" );
    generic.add_options()
        ("authors,a", "prints the authors list")
        ("copyright,c", "prints the copyright statement")
        ("help,h", "prints this help message")
        ("license,l", "prints the license text")
        ("version,v", "prints the version")
        ("verbose,V", "verbose mode")
        ;
    _M_desc.add( generic );

    po::store(po::parse_command_line(argc, argv, _M_desc), _M_vm);
    po::notify(_M_vm);

    processGenericOptions();

}
Application::Application( int argc, char** argv, AboutData const& ad, po::options_description const& od  )
    :
    _M_about( ad ),
    _M_desc( "Allowed options" ),
    _M_vm()
{
    po::options_description generic( "Generic options" );
    generic.add_options()
        ("authors,a", "prints the authors list")
        ("copyright,c", "prints the copyright statement")
        ("help,h", "prints this help message")
        ("license,l", "prints the license text")
        ("version,v", "prints the version")
        ("verbose,V", "verbose mode")
        ;
    _M_desc.add( generic ).add( od );

    po::store(po::parse_command_line(argc, argv, _M_desc), _M_vm);
    po::notify(_M_vm);

    processGenericOptions();

}
Application::Application( Application const& __app )
    :
    _M_about( __app._M_about ),
    _M_desc( __app._M_desc ),
    _M_vm( __app._M_vm )
{}
Application::~Application()
{}
void
Application::processGenericOptions() const
{

    if ( _M_vm.count( "verbose" ) ||
         _M_vm.count( "help" ) ||
         _M_vm.count( "version" ) ||
         _M_vm.count( "copyright" ) ||
         _M_vm.count( "license" ) ||
         _M_vm.count( "authors" ) )
        std::cout << _M_about.appName() << ": " << _M_about.shortDescription() <<  "\n";

    if ( _M_vm.count( "help" ) )
        std::cout << _M_desc << "\n";

    if ( _M_vm.count( "version" ) )
        std::cout << " version : " << _M_about.version() << "\n";
    if ( _M_vm.count( "copyright" ) )
        std::cout << " copyright : " << _M_about.copyrightStatement() << "\n";
    if ( _M_vm.count( "license" ) )
        std::cout << " license : " << _M_about.license() << "\n";
    if ( _M_vm.count( "authors" ) )
    {
        std::cout << std::setw( 30 )
                  << "Author Name"
                  << " " << std::setw( 15 )
                  << "Task"
                  << " " << std::setw( 40 )
                  << "Email Address"
                  << "\n";
        std::cout << std::setw( 85+3 ) << std::setfill( '-' ) << "\n" << std::setfill( ' ' );
        std::for_each( _M_about.authors().begin(),
                       _M_about.authors().end(),
                       std::cout
                       << std::setw( 30 )
                       << lambda::bind( &AboutPerson::name,
                                        lambda::_1 )
                       << " " << std::setw( 15 )
                       << lambda::bind( &AboutPerson::task,
                                        lambda::_1 )
                       << " " << std::setw( 40 )
                       << lambda::bind( &AboutPerson::emailAddress,
                                        lambda::_1 )
                       << "\n");
    }
}
}
