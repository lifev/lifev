/* -*- mode: c++ -*-
   This program is part of the LifeV library
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
  \file gracePlot.cpp
  \author Michel Kern (Inria, Estime) (copied and adapted by V. Martin)
  \date 07/2004
  \version 1.0
*/

#include <cmath>
#include <sstream>

#include <life/lifefilters/gracePlot.hpp>

namespace LifeV
{
#if defined(HAVE_GRACE_NP_H)
void GracePlot::GraceInit()
{
    if ( GraceOpenVA("xmgrace", 2048, "-nosafe", "-noask", NULL)== -1) {
        std::cerr << "Can't run grace" << std::endl;
        exit (EXIT_FAILURE);
    }
    GracePrintf("g0 on");                     /* Activate graph 0 */
    GracePrintf("with g0");                   /* reset the current graph to graph 0 */
    GracePrintf("view 0.1, 0.1, 1.2, 0.9");
}
GracePlot::GracePlot()
    :
    _M_do_plot( false )
{
    std::cout << "do_plot: " << std::boolalpha << _M_do_plot << "\n";
}

GracePlot::GracePlot( GetPot const& __data )
    :
    _M_do_plot( __data( "miscellaneous/show_graceplot",  false ) )
{
    std::cout << "do_plot: " << std::boolalpha << _M_do_plot << "\n";
    if ( doPlot() )
    {
        GraceInit();
    }
}

void GracePlot::Plot(const Rn& x, const Rn& y)
{
    if (  doPlot() )
    {
        int n = x.N();
        ASSERT( y.N() == n,
                "Plot: x and y should have same size." );

        GracePrintf("with g0");
        GracePrintf("kill s0");
        for (int i=0; i<n; i++)
            GracePrintf ("g0.s0 point %g, %g", x(i), y(i));
        GracePrintf(" ");


        GracePrintf("autoscale");
        GracePrintf("redraw");
    }
}

void GracePlot::Plot(const std::vector< Point1D >& x,
                     const ScalUnknown<Vector>& y)
{
    if (  _M_do_plot )
    {
        UInt n = x.size();
        ASSERT( y.size() == n,
                "Plot: x and y should have same size." );

        GracePrintf("with g0");
        GracePrintf("kill s0");
        for (UInt ii=0; ii<n; ii++)
            GracePrintf ("g0.s0 point %g, %g", x[ii].x(), y(ii));
        GracePrintf(" ");


        GracePrintf("autoscale");
        GracePrintf("redraw");
    }
}


#ifdef TEST
double D=0.2;
double a=2.;
double b=0.1;

double f(double x) {return sin(2*M_PI*x)*sin(2*M_PI*x);}

double solex2(double x, double t) {
    return exp((b-D)*t)*sin(x-a*t);
}

double solex1(double x, double t) {
    double tmp=sqrt(a*a+4*b*D);
    return t==0
        ? 0
        : 0.5*exp(a*x/(2.*D))*(exp(-x/(2.*D)*tmp)*erfc((x-tmp*t)/(2*sqrt(D*t)))
                               + exp(x/(2*D)*tmp)*erfc((x+tmp*t)/(2*sqrt(D*t))));
}

int main(int argc, char* argv[]) {

    int n=64;
    double L=6;
    double dx=L/(n-1);
    double tf=1;
    double dt=dx/a;

    Rn x(n), y(n);
    GracePlot p;

    for (int i=0; i<n; i++) {
        x(i) = i*dx;
    }

    for (double t=0; t<tf+dt; t+=dt) {
        for (int i=0; i<n; i++) {
            y(i) = solex1(x(i), t);
        }
        std::ostringstream str;
        str << "Time= " << t;
        p.Title(str.str());
        p.Xlabel("x"); p.Ylabel("f(x)");
        p.Legend("courbe de y=solex(x)");
        p.Plot(x, y);
        p.Sleep(0.5);
    }

    std::cout << "Hit return to close plot" << std::endl;
    char ch;
    cin.get(ch);
}


#endif /* TEST */

#endif /* HAVE_GRACE_NP_H */

} /* end of namespace LifeV */
