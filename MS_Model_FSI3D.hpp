//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief A short description of the file content

    @author Paolo Crosetto <paolo.crosetto@epfl.ch>
    @date 19 Apr 2010

    A more detailed description of the file (if necessary)
 */

#ifndef MS_MODEL_FSI3D_H
#define MS_MODEL_FSI3D_H 1

#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

#include <lifemc/lifesolver/BCInterface.hpp>

#include <life/lifesolver/FSISolver.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/noexport.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/hdf5exporter.hpp>
#endif


namespace LifeV {

//! MS_Model_FSI3D -
/*!
    @author Paolo Crosetto
    @see Reference to papers (if available)

    Here write a long and detailed description of the class.

    For this purpose you can use a lot of standard HTML code.
    Here there is a list with some useful examples:

    For bold text use: <b>BOLD TEXT</b>

    For empatyze a word type @e followed by the word

    For verbatim a word type @c followed by the word

    For vertical space (empty lines) use: <br>

    For creating list type:
    <ol>
        <li> First element of the enumerated list
        <ul>
             <li> First element of the dotted sublist.
             <li> Second element of the dotted sublist
        </ul>
        <li> Second element of the enumerated list
        <li> Third element of the enumerated list
    </ol>

    For writing a warning type: @warning followed by the description
    of the warning

    It is possible to use a lot of other standard HTML commands.
    Visit http://www.stack.nl/~dimitri/doxygen/htmlcmds.html for
    a detailed list.

    For any other kind of information visit www.doxygen.org.
 */
class MS_Model_FSI3D:
       public virtual MS_PhysicalModel
{
public:

    //! @name Public Types
    //@{

    typedef MS_PhysicalModel                                                               super;
    typedef FSISolver                                                                      FSISolverType;
    typedef boost::shared_ptr< FSISolverType >                                             FSISolverPtrType;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh3D<LifeV::LinearTetra> > > filter_ptrtype;
    typedef LifeV::Ensight<LifeV::FSIOperator::mesh_type>                                  ensightfilter_type;
    typedef boost::shared_ptr<ensightfilter_type>                                          ensightfilter_ptrtype;
#ifdef HAVE_HDF5
    typedef LifeV::Hdf5exporter<LifeV::FSIOperator::mesh_type>                             hdf5filter_type;
    typedef boost::shared_ptr<hdf5filter_type>                                             hdf5filter_ptrtype;
#endif
    typedef RegionMesh3D<LinearTetra>                                                      mesh_type;
    typedef OseenShapeDerivative< mesh_type >                                              fluid_type;
    typedef VenantKirchhofSolver< mesh_type >                                              solid_type;
    typedef FSISolverType::vector_type                                                     vector_type;
    typedef FSISolverType::vector_ptrtype                                                  vector_ptrtype;
    typedef BCInterface< FSIOperator >                                                     fluid_bc_type;
    typedef BCInterface< FSIOperator >                                                     solid_bc_type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor
    MS_Model_FSI3D(FSISolverPtrType solver);

    MS_Model_FSI3D();

private:
    MS_Model_FSI3D( const MS_Model_FSI3D& FSI3D ){}
public:
    //! Destructor
    ~MS_Model_FSI3D();

    //@}



    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    void SetupData( const std::string& FileName );

    //! Setup the model.
    void SetupModel();

    //! Build the initial system (matrix and vectors).
    void BuildSystem();

    //! Update the system for (matrix and vectors).
    void UpdateSystem();

    //! Solve the problem.
    void SolveSystem();

    //! Save the solution
    void SaveSolution();

    //! Display some information about the model.
    void ShowMe();

    //@}


private:

    void
    SetUpExporter( filter_ptrtype& exporter, const std::string name, const GetPot& dataFile );
    void
    SetExporterFluid(const filter_ptrtype& exporter);
    void
    SetExporterSolid(const filter_ptrtype& exporter);
    void updateBCSegregated(const std::string& dataFileName);
    void updateBCMonolithic(const std::string& dataFileName);

    filter_ptrtype                        M_exporterSolid;
    filter_ptrtype                        M_exporterFluid;
    FSISolverPtrType                      M_solver;
    vector_ptrtype                        M_solidDisp;
    vector_ptrtype                        M_solidVel;
    vector_ptrtype                        M_velAndPressure;
    vector_ptrtype                        M_fluidDisp;
    Real                                  M_time;
    boost::shared_ptr< fluid_bc_type >       M_fluidBC;
    boost::shared_ptr< solid_bc_type >       M_solidBC;
    boost::shared_ptr< solid_bc_type >       M_harmonicExtensionBC;
    boost::shared_ptr< fluid_bc_type >       M_linearizedFluidBC;
    boost::shared_ptr< solid_bc_type >       M_linearizedSolidBC;
};

inline MS_PhysicalModel* createModelFSI3D()
{
    return new MS_Model_FSI3D();
}

} // Namespace LifeV

#endif /* MS_MODEL_FSI3D_H */
