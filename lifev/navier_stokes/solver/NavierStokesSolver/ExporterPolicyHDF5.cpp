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

#include <lifev/navier_stokes/solver/NavierStokesSolver/ExporterPolicyHDF5.hpp>

#include <string>
#include <lifev/core/util/LifeChrono.hpp>

namespace LifeV
{

void
ExporterPolicyHDF5::initExporter ( Teuchos::ParameterList& list,
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

void
ExporterPolicyHDF5::exportSolution ()
{
    displayer().leaderPrint ( "Exporting solution at time t =  ", currentTime(), "... \n" );
    M_exporter->postProcess ( currentTime() );
}

void
ExporterPolicyHDF5::finalizeExporter()
{
    LifeChrono finalizeChrono;
    finalizeChrono.start();
    M_exporter->closeFile();
    finalizeChrono.stop();
    displayer().leaderPrintMax ("Exporter finalization time: ", finalizeChrono.diff(), " s.\n");
}

} // namespace LifeV
