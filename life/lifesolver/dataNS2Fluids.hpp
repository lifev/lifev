/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Daniele A. Di Pietro <dipietro@unibg.it>
       Date: 31-1-2005

  Copyright (C) 2005 Università degli Studi di Bergamo

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
   \file dataNS2Fluids.hpp
   \author Daniele A. Di Pietro <dipietro@unibg.it>
   \date 31-1-2005
 */

#ifndef _DATANS2FLUIDS_H_
#define _DATANS2FLUIDS_H_

#include <life/lifesolver/dataNavierStokes.hpp>

namespace LifeV {

    /*!
      \class dataNS2Fluids
      \brief Class which holds data for two fluid Navier-Stokes solvers.

      @author Daniele Antonio Di Pietro <dipietro@unibg.it>
      @see
    */

    template<typename MeshType>
    class DataNS2Fluids
        :
        public DataMesh<MeshType>,
        public DataTime
    {
    public:

        /** @name Typedefs
         */
        //@{
        typedef MeshType mesh_type;
        typedef enum {fluid1 = 1, fluid2 = 2} fluid_type;
        //@}

        /** @name Constructors and destructors
         */
        //@{
        DataNS2Fluids(const GetPot& datafile)
            :
            DataMesh<mesh_type>( datafile, "navier-stokes/discretization" ),
            DataTime( datafile, "navier-stokes/time" )
        {
            // Physics

            _M_rho_1 = datafile( "navier-stokes/physics/fluid_1/density", 1. );
            _M_mu_1 = datafile( "navier-stokes/physics/fluid_1/viscosity", 1. );

            _M_rho_2 = datafile("navier-stokes/physics/fluid_2/density", 1.);
            _M_mu_2 = datafile("navier-stokes/physics/fluid_2/viscosity", 1.);

            // Miscellaneous

            _M_verbose = datafile( "navier-stokes/miscellaneous/verbose", 1 );
        }

        ~DataNS2Fluids() {}
        //@}

        /** @name Accessors
         */
        //@{

        /**
           \Return the density of fluid_id-th fluid
        */

        inline Real density(fluid_type fluid_id) const {
            return (fluid_id == 1) ? _M_rho_1 : _M_rho_2;
        }

        /**
           \Return the viscosity of fluid_id-th fluid
        */

        inline Real viscosity(fluid_type fluid_id) const {
            return (fluid_id == 1) ? _M_mu_1 : _M_mu_2;
        }

        /**
           \Return 1 if verbose mode is set, 0 otherwise
        */
        bool verbose() const {
            return _M_verbose;
        }
        //@}

        /** @name Methods
         */
        //@{
        template<typename _MeshType>
        friend std::ostream& operator<<(std::ostream&, DataNS2Fluids<_MeshType>&);

        //! Fake member functions introduced for compliance with the new
        //! interface requirements of NavierStokesHandler class

        UInt computeMeanValuesPerSection() const {
            return 0;
        }

        UInt NbZSections() const {
            return 0;
        }

        Real ToleranceSection() const {
            return 0.;
        }

        Real XSectionFrontier() const {
            return 0.;
        }

        Real ZSectionInit() const {
            return 0.;
        }

        Real ZSectionFinal() const {
            return 0;
        }

        UInt NbPolygonEdges() const {
            return 0;
        }

        //@}
    private:
        //! Physics

        Real _M_rho_1;
        Real _M_mu_1;
        Real _M_rho_2;
        Real _M_mu_2;

        //! Miscellaneous

        bool _M_verbose;

    }; // class DataNS2Fluids

    template<typename _MeshType>
    std::ostream& operator<<(std::ostream& ostr, DataNS2Fluids<_MeshType>& d) {
        typedef typename DataNS2Fluids<_MeshType>::fluid_type fluidType;
        fluidType f1 = DataNS2Fluids<_MeshType>::fluid1;
        fluidType f2 = DataNS2Fluids<_MeshType>::fluid2;

        ostr << "==================================================" << std::endl;
        ostr << "Data for two-fluid Navier-Stokes problem" << std::endl;
        ostr << "==================================================" << std::endl;
        ostr << "Fluid 1" << std::endl;
        ostr << "--------------------------------------------------" << std::endl;
        ostr << " rho_1 = " << d.density(f1) << std::endl;
        ostr << " mu_1  = " << d.viscosity(f1) << std::endl;
        ostr << "--------------------------------------------------" << std::endl;
        ostr << "Fluid 2" << std::endl;
        ostr << "--------------------------------------------------" << std::endl;
        ostr << " rho_2 = " << d.density(f2) << std::endl;
        ostr << " mu_2  = " << d.viscosity(f2) << std::endl;
        ostr << "--------------------------------------------------" << std::endl;
        ostr << "Time span" << std::endl;
        ostr << "--------------------------------------------------" << std::endl;
        ostr << " t0    = " << d.getInitialTime() << std::endl;
        ostr << " T     = " << d.getEndTime() << std::endl;
        ostr << "--------------------------------------------------" << std::endl;
        ostr << "Miscellaneous" << std::endl;
        ostr << "--------------------------------------------------" << std::endl;
        ostr << " verbose ";
        d.verbose() ? ostr << "yes" << std::endl : ostr << "no" << std::endl;
        ostr << "==================================================" << std::endl;

        return ostr;
    } // operator <<

} // namespace LifeV

#endif
