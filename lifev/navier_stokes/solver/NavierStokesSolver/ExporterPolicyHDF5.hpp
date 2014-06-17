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
    @file ExporterPolicyHDF5 class
    @brief This class contains all the informations necessary to export to HDF5

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 04-12-2012
 */

#ifndef EXPORTERPOLICYHDF5_HPP
#define EXPORTERPOLICYHDF5_HPP

#include <iostream>
#include <string>
#include <boost/shared_ptr.hpp>


#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/util/LifeChrono.hpp>


namespace LifeV
{

template< class mesh_Type >
struct ExporterPolicyHDF5
{

    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef boost::shared_ptr<mesh_Type>             meshPtr_Type;
    typedef MapEpetra                                map_Type;
    typedef boost::shared_ptr<map_Type>              mapPtr_Type;
    typedef FESpace< mesh_Type, map_Type >           fespace_Type;
    typedef boost::shared_ptr< fespace_Type >        fespacePtr_Type;
    typedef ExporterHDF5<mesh_Type>                  exporter_Type;
    typedef boost::shared_ptr< exporter_Type >       exporterPtr_Type;

    void initExporter ( Teuchos::ParameterList& list,
                        vectorPtr_Type solution );
    void exportSolution ();
    void finalizeExporter();
    exporterPtr_Type M_exporter;

    // Host state variables getters
    virtual Displayer displayer() = 0;
    virtual meshPtr_Type mesh() const = 0;
    virtual fespacePtr_Type uFESpace() const = 0;
    virtual fespacePtr_Type pFESpace() const = 0;
    virtual Real currentTime() const = 0;

};

template< class mesh_Type >
void
ExporterPolicyHDF5< mesh_Type >::initExporter ( Teuchos::ParameterList& list,
                                                vectorPtr_Type solution )
{
    // Loading the parameters
    std::string outputPath = list.get ( "Output path", "." );
    outputPath.append ("/");
    std::string outputFilename = list.get ( "Output filename", "solution" );

    GetPot datafile;
    bool multipleMesh = list.get ( "Multiple mesh", false );
    if ( multipleMesh )
    {
        datafile.set ( "exporter/multimesh", "true" );
    }
    else
    {
        datafile.set ( "exporter/multimesh", "false" );
    }
    int start = list.get ( "Start", 0 );
    datafile.set ( "exporter/start", start );
    int save = list.get ( "Save", 1 );
    datafile.set ( "exporter/save", save );

    LifeChrono exporterSetupChrono;
    exporterSetupChrono.start();

    displayer().leaderPrint ( "Defining the exporter... " );
    M_exporter.reset ( new exporter_Type ( datafile, outputFilename ) );
    M_exporter->setPostDir ( outputPath ); // This is a test to see if M_post_dir is working
    M_exporter->setMeshProcId ( mesh(), mesh()->comm()->MyPID() );
    displayer().leaderPrint ( "done\n" );

    displayer().leaderPrint ( "Updating the exporter... " );

    // Pressure offset in the vector
    UInt pressureOffset = uFESpace()->fieldDim() * uFESpace()->dof().numTotalDof();

    M_exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace(),
                              solution, UInt ( 0 ) );
    M_exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace(),
                              solution, pressureOffset );
    displayer().leaderPrint ( "done\n" );

    exporterSetupChrono.stop();
    displayer().leaderPrintMax ("Exporter setup time: ", exporterSetupChrono.diff(), " s.\n");
}

template< class mesh_Type >
void
ExporterPolicyHDF5< mesh_Type >::exportSolution ()
{
    displayer().leaderPrint ( "Exporting solution at time t =  ", currentTime(), "... \n" );
    M_exporter->postProcess ( currentTime() );
}

template< class mesh_Type >
void
ExporterPolicyHDF5< mesh_Type >::finalizeExporter()
{
    LifeChrono finalizeChrono;
    finalizeChrono.start();
    M_exporter->closeFile();
    finalizeChrono.stop();
    displayer().leaderPrintMax ("Exporter finalization time: ", finalizeChrono.diff(), " s.\n");
}

} // namespace LifeV

#endif /* EXPORTERPOLICYHDF5_HPP */
