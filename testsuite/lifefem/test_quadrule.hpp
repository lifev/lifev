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

/*!
    @file
    @brief Quadrature Rule test

	@author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @author Umberto Villa <uvilla@emory.edu>
    @contributor
    @maintainer Umberto Villa <uvilla@emory.edu>

    @date 02-03-2010

The program tests the degree of exactness (DoE) and the convergence rate (CR)
of all the quadrature rule in 3D (Tetrahedra) or in 2D (Triangles).

The code produce the following output:

quadRuleTetra.txt ==> Show the Degree of Exactness of all the quadrature rules on Tetrahedral elements
quadRuleTetra.plt ==> Show the Convergence Rate of all the quadrature rules on Tetrahedral elements
                               using gnuplot
 */

#include <life/lifecore/life.hpp>
#include <life/lifefem/quadRule.hpp>
#include <string>
#include <fstream>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/numeric/ublas/vector.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include "SetOfFun.hpp"
namespace LifeV
{
typedef boost::numeric::ublas::vector<Real> Vector;
typedef std::vector<QuadRule const *> container_Type;
typedef container_Type::const_iterator constIterator_Type;

// This function checks the DEGREE of Exactness (DoE) of the quadrature rules
template<typename Mesh>
bool quad_check_doe(const RefFE &refFE, const GeoMap & geoMap, const container_Type &allQuad, std::string output_file)
{

    Mesh aMesh;
    UInt nEl(1);
    regularMesh3D( aMesh, 1, nEl, nEl, nEl);

    SetofFun fct;
    int nquadrule = allQuad.size();
    bool check(false);
    output_file = output_file+".txt";
    std::ofstream ofile(output_file.c_str());


    for (UInt nqr(0); nqr<nquadrule; ++nqr)
    {

        CurrentFE fe(refFE,geoMap, *allQuad[nqr]);
        ofile<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_"<<std::endl;
        ofile<<allQuad[nqr]->name()<<std::endl;
        ofile<<"Degree of Exactness: "<<allQuad[nqr]->degreeOfExactness()<<std::endl;
        ofile<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_"<<std::endl;

        ofile<<"polynomial \t\t test \t quadrature error"<<std::endl;

        for (UInt fun(0); fun<fct.nfun(); ++fun)
        {
            Real integral = 0.;
            for (UInt i=1; i<=aMesh.numElements(); ++i)
            {
                fe.updateJacQuadPt(aMesh.element(i));
                Real s = 0., x, y, z;
                for ( UInt ig(0); ig < fe.nbQuadPt(); ++ig )
                {
                    fe.coorQuadPt( x, y, z, ig );
                    s += fct.val( fun, x, y, z) * fe.weightDet( ig );
                }
                integral+=s;
            }
            Real err;
            err = integral-fct.ex_int(fun);
            // Check for Quadrature Rule Errors
            if (fct.degree(fun) <= allQuad[nqr]->degreeOfExactness())
            {
                if (fabs(err)<1e-9)
                {
                    check = true;
                    ofile<<fct.name(fun)<<" \t passed \t ("<<err<<")"<<std::endl;
                }
                else
                {
                    check = false;
                    ofile<<fct.name(fun)<<" \t FAILED \t ("<<err<<")"<<std::endl;
                }
            }

        }//end for on fun
        ofile<<std::endl<<std::endl;
    }//end for on qr
    return check;
}

// This function checks the convergence rate (CR) of the quadrature rules
template<typename Mesh>
bool quad_check_cr(	const RefFE &refFE, const GeoMap & geoMap, const container_Type &allQuad, std::string output_name)
{
    SetofFun fct;
    int fun(fct.nfun());

    int nrefine(4);
    Vector h(nrefine);

    for (int i=0; i<nrefine; ++i)
        h(i) = pow(.5,i);

    int nquadrule = allQuad.size();
    boost::numeric::ublas::matrix<double> err(nrefine,nquadrule+1);

    for (int iref = 0; iref<nrefine; ++iref)
    {

        Mesh aMesh;
        UInt nEl( pow(2.0,static_cast<double>(iref) ) );
        regularMesh3D( aMesh, 1, nEl, nEl, nEl);

        err(iref,0) = h(iref);
        for (int iqr = 0; iqr<nquadrule; ++iqr)
        {

            // ===================================================
            // Current FE classes for the problem under study with
            // mapping and quadrature rules
            // ===================================================

            CurrentFE fe(refFE,geoMap, *allQuad[iqr]);
            Real integral = 0.;
            for (UInt i=1; i<=aMesh.numElements(); ++i)
            {
                fe.updateJacQuadPt(aMesh.element(i));
                Real s = 0., x, y, z;
                for ( UInt ig(0); ig < fe.nbQuadPt(); ++ig )
                {
                    fe.coorQuadPt( x, y, z, ig );
                    s += fct.val( fun, x, y, z) * fe.weightDet( ig );
                }
                integral+=s;
            }
            err(iref,iqr+1) = fabs(integral-fct.ex_int(fun));


        }//end for on qr
    }//end for on iref

    std::string fname;
    {
        fname = output_name+".dat";
        std::ofstream ofile(fname.c_str());
        for (int i = 0; i<nrefine; ++i)
        {
            for (int j=0; j<nquadrule+1; ++j)
                ofile<<err(i,j)<<" \t";
            ofile<<std::endl;
        }
    }
    {
        fname = output_name+".plt";
        std::ofstream ofile(fname.c_str());
        ofile<<"set timestamp"<<std::endl;
        ofile<<"set logscale"<<std::endl;
        ofile<<"plot ";
        for (int i = 0; i<nquadrule; ++i)
            ofile<<"\""<<output_name<<".dat\" using 1:"<<i+2<<" with linespoints title '"<<allQuad[i]->name()<<"' ,\\"<<std::endl;
        ofile<<"\""<<output_name<<".dat\" using 1:($1*$1) with lines title 'order 2', \\"<<std::endl;
        ofile<<"\""<<output_name<<".dat\" using 1:($1**4) with lines title 'order 4', \\"<<std::endl;
        ofile<<"\""<<output_name<<".dat\" using 1:($1**6) with lines title 'order 6'";
    }

    return true;

}

}//end namespace
