/* -*- mode: c++ -*-
   This program is part of the LifeV library

   Authors: M. Kern (Inria, Estime) -- initial version
            V. Martin -- adaptation to LifeV
            C. Prud'homme -- portability fixes, ie when grace is not installed

   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file gracePlot.hpp
  \author Michel Kern (Inria, Estime)
  \date 07/2004
  \version 1.0

  \brief This file contains an interface for graceplot.
         Useful for plotting 1D solutions.
     See http://graceplot.sourceforge.net/

     When grace is not installed/available, then a dummy class that won't do
     anything is instantiated, a warning is written on the standard error output.
*/

#ifndef _GRACEPLOT_H_
#define _GRACEPLOT_H_

#include <lifeconfig.h>

#include <string>
#include <vector>

#if defined(HAVE_GRACE_NP_H)
#include <grace_np.h>
#endif /* HAVE_GRACE_NP_H */

#include <life/lifecore/GetPot.hpp>


#include <life/lifearray/RNM.hpp>
#include <life/lifearray/vecUnknown.hpp>
#include <life/lifemesh/basicOneDMesh.hpp>


namespace LifeV
{
typedef KN<double> Rn;


#if defined(HAVE_GRACE_NP_H)
class GracePlot {

public:
    GracePlot();
    GracePlot( GetPot const& );

    ~GracePlot()
        {
            if (  doPlot() )
            {
                GraceClose();
            }
        };

    bool doPlot() const { return _M_do_plot; }

    void Title(std::string title)
        {
            if (  doPlot() )
            {
                GracePrintf("with g0");
                GracePrintf("title \"%s\"", title.data());
            }
        }

    void Legend(std::string label)
        {
            if (  doPlot() )
            {
                GracePrintf("with g0");
                GracePrintf("s0 legend \"%s\"", label.data());
            }
        }

    void Xlabel(std::string xlabel)
        {
            if (  doPlot() )
            {
                GracePrintf("with g0");
                GracePrintf("xaxis label \"%s\"", xlabel.data());
            }
        }

    void Ylabel(std::string ylabel)
        {
            if (  doPlot() )
            {
                GracePrintf("with g0");
                GracePrintf("yaxis label \"%s\"", ylabel.data());
            }
        }

    void Plot(const Rn& x, const Rn& y);

    void Plot(const std::vector< Point1D >& x,
              const ScalUnknown<Vector>& y);

    void Sleep(double delay)
        {
            if (  doPlot() )
            {
                GracePrintf("sleep %g", delay);
            }
        }

    void setPlot( bool __do_plot )
        {
            _M_do_plot = __do_plot;
            if (  doPlot() )
                GraceInit();
            else
                GraceClose();
        }

private:
    void GraceInit();
private:
    bool _M_do_plot;
};
#else
/*!\class GracePlot
  \brief dummy class for grace plots

  Dummy class for grace plots when grace is not available.

  @see http://graceplot.sourceforge.net/ for further information about grace
*/
class GracePlot {

public:
    GracePlot()
        {
            std::cerr << "================================================================================\n"
                      << "WARNING: you are using the dummy graceplot class\n"
                      << "because grace is not installed on your system\n"
                      << "see http://graceplot.sourceforge.net/ to install grace\n"
                      << "================================================================================\n";
        }
    GracePlot( GetPot const& )
        {
            std::cerr << "================================================================================\n"
                      << "WARNING: you are using the dummy graceplot class\n"
                      << "because grace is not installed on your system\n"
                      << "see http://graceplot.sourceforge.net/ to install grace\n"
                      << "================================================================================\n";
        }
    ~GracePlot(){ };

    bool doPlot() const { return false; }

    void Title(std::string /*title*/) {}

    void Legend(std::string /*label*/) {}

    void Xlabel(std::string /*xlabel*/) {}

    void Ylabel(std::string /*ylabel*/) {}

    void Plot(const Rn& x, const Rn& y);

    void Plot(const std::vector< Point1D >& /*x*/,
              const ScalUnknown<Vector>& /*y*/) {}

    void Sleep(double /*delay*/) {}
    //  void Hold(bool onoff);
    //  void Clear();

};
#endif /* HAVE_GRACE_NP_H */
}


#endif
