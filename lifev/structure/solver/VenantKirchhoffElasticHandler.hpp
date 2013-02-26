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
 *  @file
 *  @brief This file contains an abstract class for Elastic Structures.
 *
 *  @version 1.0
 *  @date 01-06-2003
 *  @author Miguel Angel Fernandez
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _ELASTICSTRUCTUREHANDLER_H_
#define _ELASTICSTRUCTUREHANDLER_H_

#include <cmath>
#include <sstream>
#include <cstdlib>

#include <life/lifesolver/dataElasticStructure.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/core/fem/DOF.hpp>
#include <lifev/core/LifeV.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#include <lifev/core/fem/BCHandler.hpp>

namespace LifeV
{
/*!
  \class VenantKirchhoffElasticHandler
*/

template <typename Mesh>
class VenantKirchhoffElasticHandler:
    public VenantKirchhoffElasticData<Mesh>
{

public:
    //! @name Type definitions
    //@{

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > source_Type;
    typedef BCHandler                                bchandlerRaw_Type;
    typedef boost::shared_ptr<bchandlerRaw_type>     bchandler_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    VenantKirchhoffElasticHandler();

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature
      \param BCh boundary conditions for the displacement
    */
    setUp ( const VenantKirchhoffElasticData<Mesh>& data,
            const RefernceFE&                      refFE,
            const QuadratureRule&                   Qr,
            const QuadratureRule&                   bdQr,
            BCHandler&                        BCh );

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature
    */
    setUp ( const VenantKirchhoffElasticData<Mesh>& data,
            const ReferenceFE&                      refFE,
            const QuadratureRule&                   Qr,
            const QuadratureRule&                   bdQr);

    //! Destructor
    virtual ~VenantKirchhoffElasticHandler()
    {}

    //@}

    //! @name Methods
    //@{

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic force
      \param time present time
    */
    virtual void timeAdvance ( source_Type const& , const Real& time ) = 0;

    //! Solve the non-linear problem
    virtual void iterate() = 0;

    //! Postprocessing
    void postProcess();

    //! Sets initial condition for the displacment en velocity
    /*!
      \param d0 space function defining the initial displacement
      \param w0 space function defining the initial velocity
    */
    void initialize ( const Function& d0, const Function& w0 );

    /*!
      \param depName file containing the initial displacement vector
      \param velName file containing the initial velocity vector
      \param starT starting time when the displacement and the velocity vectors
             are set
    */
    void initialize ( const std::string& depName,
                      const std::string& velName,
                      Real             startT = 0.);

    //@}

    //! @name Set methods
    //@{

    //! checking if BC are set
    const bool setSolidBC() const
    {
        return M_setBC;
    }

    //! set the fluid BCs
    /*!
      \param BCh_solid BCHandler object containing the boundary conditions
    */
    void setSolidBC (BCHandler& BCh_solid)
    {
        M_BCh_solid = &BCh_solid;
        M_setBC = true;
    }

    //@}

    //! @name Get methods
    //@{

    //! Returns the displacement vector
    PhysVectUnknown<Vector>& getDisplacement();

    //! Returns the velocity vector
    PhysVectUnknown<Vector>& getVelocity();

    //! Returns the number of unknowns
    const UInt getDimension() const
    {
        return M_dim;
    }

    //! Returns the reference FE object
    const ReferenceFE&          getRefFE() const
    {
        return M_refFE;
    }

    //! Returns the current FE object
    CurrentFE&            getFe()
    {
        return M_fe;
    }
    //! Returns the current FE object
    const CurrentFE&      getFe() const
    {
        return M_fe;
    }
    //! Returns the current boundary FE object
    CurrentBoundaryFE&    getFeBd()
    {
        return M_feBd;
    }
    //! Returns the DOF object
    const DOF&            getdDof() const
    {
        return M_dof;
    }

    //! Returns the BCHandler object
    BCHandler& getBCh_solid()
    {
        return *M_BCh_solid;
    }
    //@}

private:

    //! @name Provate Methods
    /*!@{

    //!private copy constructor:this class should not be copied
    /*!
      if you need a copy you should implement it, so that it copies the shared pointer one by one, without copying the content.
    */
    VenantKirchhoffElasticHandler (VenantKirchhoffElasticHandler& T);
    {}

    /*!
      \param name string containing the name of the file
      \param unknown vector that must be initialized to the vector contained in the file
    */
    void readUnknown ( const std::string&       name,
                       PhysVectUnknown<Vector>& unknown);

    //! Reference FE
    const ReferenceFE&                           M_refFE;

    //! The DOF object
    DOF                                    M_dof;

    //! The number of total displacement dofs
    UInt                                   M_dim;

    //! Quadrature rule for volumic elementary computations
    const QuadratureRule&                        M_Qr;

    //! Quadrature rule for elementary computations
    const QuadratureRule&                        M_bdQr;

    //! Current FE
    CurrentFE                              M_fe;

    //! Current boundary FE
    CurrentBoundaryFE                            M_feBd;

    //! The displacement
    PhysVectUnknown<Vector>                M_d;
    PhysVectUnknown<Vector>                M_dRhs;

    //! The velocity
    PhysVectUnknown<Vector>                M_w;

    //! The current time
    Real                                   M_time;

    //! Aux. var. for PostProc
    UInt                                   M_count;

    //! The BC handler
    BCHandler*                             M_BCh_solid;

    bool                                   M_setBC;

};


//=========================================
// Constructor
//=========================================

template <typename Mesh>
VenantKirchhoffElasticHandler<Mesh>::
VenantKirchhoffElasticHandler( ) :
    VenantKirchhoffElasticData<Mesh> (  ),
    M_refFE                    (  ),
    M_dof                      (  ),
    M_dim                      (  ),
    M_Qr                       (  ),
    M_bdQr                     (  ),
    M_fe                       (  ),
    M_feBd                     (  ),
    M_d                        (  ),
    M_dRhs                     (  ),
    M_w                        (  ),
    M_time                     (  ),
    M_count                    (  ),
    M_BCh_solid                (  )
{}

template <typename Mesh>
VenantKirchhoffElasticHandler<Mesh>::
setUp ( VenantKirchhoffElasticData<Mesh>& data,
        const ReferenceFE&                      refFE,
        const QuadratureRule&                   Qr,
        const QuadratureRule&                   bdQr,
        BCHandler&                        BCh )
{
    VenantKirchhoffElasticData<Mesh> = data;
    M_refFE                    = refFE;
    M_dof                      ( this->mesh(), M_refFE ),
                               M_dim                      = M_dof.numTotalDof();
    M_Qr                       = Qr;
    M_bdQr                     = bdQr;
    M_fe                       ( M_refFE, getGeometricMap ( this->mesh() ), M_Qr )
    M_feBd                     ( M_refFE.boundaryFE(), getGeometricMap ( this->mesh() ).boundaryMap(), M_bdQr );
    M_d                        = M_dim;
    M_dRhs                     = M_dim;
    M_w                        = M_dim;
    M_time                     = 0;
    M_count                    = 0;
    M_BCh_solid                = &BCh;
}

template <typename Mesh>
VenantKirchhoffElasticHandler<Mesh>::
setUp ( const VenantKirchhoffElasticData<Mesh>& data,
        const ReferenceFE&                      refFE,
        const QuadratureRule&                   Qr,
        const QuadratureRule&                   bdQr)
{
    VenantKirchhoffElasticData<Mesh> = data;
    M_refFE                    = refFE;
    M_dof                      ( this->mesh(), M_refFE ),
                               M_dim                      = M_dof.numTotalDof();
    M_Qr                       = Qr;
    M_bdQr                     = bdQr;
    M_fe                       ( M_refFE, getGeometricMap ( this->mesh() ), M_Qr );
    M_feBd                     ( M_refFE.boundaryFE(), getGeometricMap ( this->mesh() ).boundaryMap(), M_bdQr );
    M_d                        = M_dim;
    M_dRhs                   = M_dim;
    M_w                        = M_dim ;
    //    _BCh( new BCHandler(0)),
    M_time                     = 0;
    M_count                    = 0;
    M_BCh_solid                = 0;
}

// Returns the displacement vector
template <typename Mesh>
PhysVectUnknown<Vector>&
VenantKirchhoffElasticHandler<Mesh>::disp()
{
    return M_d;
}


// Returns the velocity vector
template <typename Mesh>
PhysVectUnknown<Vector>&
VenantKirchhoffElasticHandler<Mesh>::w()
{
    return M_w;
}

// Postprocessing
template <typename Mesh>
void
VenantKirchhoffElasticHandler<Mesh>::postProcess()
{
    std::ostringstream index;
    std::string name, namedef;

    ++M_count;

    if ( fmod ( float ( M_count ), float ( this->_verbose ) ) == 0.0 )
    {
        std::cout << "  S-  Post-processing \n";
        index << ( M_count / this->_verbose );

        switch ( index.str().size() )
        {
            case 1:
                name = "00" + index.str();
                break;
            case 2:
                name = "0" + index.str();
                break;
            case 3:
                name = index.str();
                break;
        }

        namedef = "defor." + name + ".mesh";
        wr_medit_ascii_scalar ( "dep_x." + name + ".bb", M_d.giveVec(), this->mesh().numVertices() );
        wr_medit_ascii_scalar ( "dep_y." + name + ".bb", M_d.giveVec() + M_dim, this->mesh().numVertices() );
        wr_medit_ascii_scalar ( "dep_z." + name + ".bb", M_d.giveVec() + 2 * M_dim, this->mesh().numVertices() );
        wr_medit_ascii2 ( namedef, this->mesh(), M_d, this->_factor );

        // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),this->mesh().numVertices(),_dim_u);

        system ( ( "ln -sf " + namedef + " dep_x." + name + ".mesh" ).data() );
        system ( ( "ln -sf " + namedef + " dep_y." + name + ".mesh" ).data() );
        system ( ( "ln -sf " + namedef + " dep_z." + name + ".mesh" ).data() );

        // system(("ln -s "+this->mesh()_file+" veloc."+name+".mesh").data());

        wr_medit_ascii_scalar ( "veld_x." + name + ".bb", M_w.giveVec(), this->mesh().numVertices() );
        wr_medit_ascii_scalar ( "veld_y." + name + ".bb", M_w.giveVec() + M_dim, this->mesh().numVertices() );
        wr_medit_ascii_scalar ( "veld_z." + name + ".bb", M_w.giveVec() + 2 * M_dim, this->mesh().numVertices() );

        // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),this->mesh().numVertices(),_dim_u);

        system ( ( "ln -sf " + namedef + " veld_x." + name + ".mesh" ).data() );
        system ( ( "ln -sf " + namedef + " veld_y." + name + ".mesh" ).data() );
        system ( ( "ln -sf " + namedef + " veld_z." + name + ".mesh" ).data() );

        // system(("ln -s "+this->mesh()_file+" veloc."+name+".mesh").data());
    }
}

// Sets the initial condition
template <typename Mesh>
void
VenantKirchhoffElasticHandler<Mesh>::initialize ( const Function& d0, const Function& w0 )
{

    // Initialize velocity

    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV = M_refFE.nbDofPerVertex; // number of DOF per vertex
    UInt nDofpE = M_refFE.nbDofPerEdge;   // number of DOF per edge
    UInt nDofpF = M_refFE.nbDofPerFace;   // number of DOF per face
    UInt nDofpEl = M_refFE.nbDofPerVolume; // number of DOF per Volume

    UInt nElemV = GeoShape::S_numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::S_numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::S_numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's DOF on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's DOF on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's DOF on a Element

    ID nbComp = M_d.nbcomp(); // Number of components of the mesh velocity

    Real x, y, z;

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 0; iElem < this->mesh().numVolumes(); ++iElem )
    {

        M_fe.updateJac ( this->mesh().volume ( iElem ) );

        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on element vertices
            for ( ID iVe = 0; iVe < nElemV; ++iVe )
            {

                // Loop number of DOF per vertex
                for ( ID l = 0; l < nDofpV; ++l )
                {
                    lDof = iVe * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    M_fe.coorMap ( x, y, z, M_fe.refFE.xi ( lDof  ), M_fe.refFE.eta ( lDof ), M_fe.refFE.zeta ( lDof  ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        M_d ( icmp * M_dim + M_dof.localToGlobalMap ( iElem, lDof ) ) = d0 ( 0.0, x, y, z, icmp + 1 );
                        M_w ( icmp * M_dim + M_dof.localToGlobalMap ( iElem, lDof ) ) = w0 ( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( ID iEd = 0; iEd < nElemE; ++iEd )
            {

                // Loop number of DOF per edge
                for ( ID l = 0; l < nDofpE; ++l )
                {
                    lDof = nDofElemV + iEd * nDofpE + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    M_fe.coorMap ( x, y, z, M_fe.refFE.xi ( lDof ), M_fe.refFE.eta ( lDof ), M_fe.refFE.zeta ( lDof ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        M_d ( icmp * M_dim + M_dof.localToGlobalMap ( iElem, lDof ) ) = d0 ( 0.0, x, y, z, icmp + 1 );
                        M_w ( icmp * M_dim + M_dof.localToGlobalMap ( iElem, lDof )  ) = w0 ( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( ID iFa = 0; iFa < nElemF; ++iFa )
            {

                // Loop on number of DOF per face
                for ( ID l = 0; l < nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + iFa * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    M_fe.coorMap ( x, y, z, M_fe.refFE.xi ( lDof ), M_fe.refFE.eta ( lDof ), M_fe.refFE.zeta ( lDof ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        M_d ( icmp * M_dim + M_dof.localToGlobalMap ( iElem, lDof ) ) = d0 ( 0.0, x, y, z, icmp + 1 );
                        M_w ( icmp * M_dim + M_dof.localToGlobalMap ( iElem, lDof ) ) = w0 ( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }
        // Element based Dof
        // Loop on number of DOF per Element
        for ( ID l = 0; l < nDofpEl; ++l )
        {
            lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

            // Nodal coordinates
            M_fe.coorMap ( x, y, z, M_fe.refFE.xi ( lDof ), M_fe.refFE.eta ( lDof ), M_fe.refFE.zeta ( lDof ) );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            {
                M_d ( icmp * M_dim + M_dof.localToGlobalMap ( iElem, lDof ) ) = d0 ( 0.0, x, y, z, icmp + 1 );
                M_w ( icmp * M_dim + M_dof.localToGlobalMap ( iElem, lDof ) ) = w0 ( 0.0, x, y, z, icmp + 1 );
            }
        }
    }
}


// Sets the initial condition
template <typename Mesh>
void
VenantKirchhoffElasticHandler<Mesh>::initialize ( const std::string& depName,
                                                  const std::string& velName,
                                                  Real             startT)
{
    std::cout << "  S- restarting at time = " << startT << std::endl;

    M_count = (Int) (startT / this->timestep() - 0.5);

    // Loop on elements of the mesh
    for ( ID iElem = 0; iElem < this->mesh().numVolumes(); ++iElem )
    {
        M_fe.updateJac ( this->mesh().volume ( iElem ) );
    }
    readUnknown (depName, M_d);
    readUnknown (velName, M_w);
}


template <typename Mesh>
void
VenantKirchhoffElasticHandler<Mesh>::readUnknown ( const std::string&       name,
                                                   PhysVectUnknown<Vector>& unknown)
{
    std::string sdummy;
    std::string ext;
    UInt nsx, nsy, nsz;
    Int ndim;

    Int nDof = M_dim;

    std::string filenamex = name;
    ext = "_x.bb";
    filenamex.insert (filenamex.end(), ext.begin(), ext.end() );

    std::ifstream filex (filenamex.c_str(), std::ios::in);

    if (!filex)
    {
        std::cout << "Reading file " << filenamex
                  << " impossible" << std::endl;
        exit (1);
    }

    filex >> ndim;
    filex >> sdummy;
    filex >> nsx;

    if (nsx != M_dim)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamex << std::endl;
        exit (1);
    }

    filex >> sdummy;

    for (UInt ix = 0; ix < nsx; ++ix)
    {
        filex >> unknown[ix + 0 * nDof];
        //            unknown[ix + 0*nDof] /= factor;
    }

    filex.close();

    std::string filenamey = name;
    ext = "_y.bb";
    filenamey.insert (filenamey.end(), ext.begin(), ext.end() );

    std::ifstream filey (filenamey.c_str(), std::ios::in);

    if (!filey)
    {
        std::cout << "Reading file " << filenamey
                  << " impossible" << std::endl;
        exit (1);
    }

    filey >> ndim;
    filey >> sdummy;
    filey >> nsy;

    if (nsy != M_dim)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamex << std::endl;
        exit (1);
    }

    filey >> sdummy;

    for (UInt iy = 0; iy < nsy; ++iy)
    {
        filey >> unknown[iy + 1 * nDof];
        //        unknown[iy + 1*nDof] /= factor;
    }

    filey.close();

    std::string filenamez = name;
    ext = "_z.bb";
    filenamez.insert (filenamez.end(), ext.begin(), ext.end() );

    //     std::cout << "Reading INRIA solid file   (" << filenamez << ")"
    //               << ":" << std::endl;

    std::ifstream filez (filenamez.c_str(), std::ios::in);

    if (!filez)
    {
        std::cout << "Reading mesh file " << filenamez
                  << " impossible" << std::endl;
        exit (1);
    }

    filez >> ndim;
    filez >> sdummy;
    filez >> nsz;

    if (nsz != M_dim)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamex << std::endl;
        exit (1);
    }

    filez >> sdummy;

    for (UInt iz = 0; iz < nsz; ++iz)
    {
        filez >> unknown[iz + 2 * nDof];
        //        unknown[iz + 2*nDof] /= factor;
    }

    filez.close();
}


}

#endif
