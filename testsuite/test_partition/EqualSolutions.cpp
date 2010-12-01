/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Radu Popescu <radu.popescu@epfl.ch>
       Date: 2010-07-20

  Copyright (C) 2010 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file equalSolutions.cpp
   \author Radu Popescu <radu.popescu@epfl.ch>
   \date 2010-07-20
 */

#include "EqualSolutions.hpp"
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <EpetraExt_HDF5.h>
#include <iostream>
#include <sstream>
#include <string>

bool equalSolutions(char* fileA, char* fileB, int timesteps, double tolerance)
{
    Epetra_SerialComm comm;
    EpetraExt::HDF5 HDF5input(comm);

    bool correct(true);

    int velocity_length, dummy;
    int pressure_length;

    std::stringstream velocity_group_name, pressure_group_name;

    for (int i = 0; i < timesteps; ++i)
    {
        HDF5input.Open(fileA, H5F_ACC_RDONLY);

        velocity_group_name << "velocity.0000" << i;

        HDF5input.ReadMultiVectorProperties(velocity_group_name.str(),
                                            velocity_length,
                                            dummy);

        pressure_group_name << "pressure.0000" << i;
        HDF5input.ReadMultiVectorProperties(pressure_group_name.str(),
                                            pressure_length,
                                            dummy);

        Epetra_Map velocity_map(velocity_length, 0, comm);
        Epetra_Map pressure_map(pressure_length, 0, comm);

        Epetra_MultiVector* velocity_ref = new Epetra_MultiVector(velocity_map, 3);
        Epetra_MultiVector* pressure_ref = new Epetra_MultiVector(pressure_map, 1);


        HDF5input.Read(velocity_group_name.str(), velocity_ref, true, 0);
        HDF5input.Read(pressure_group_name.str(), pressure_ref, true, 0);

        HDF5input.Close();
        HDF5input.Open(fileB, H5F_ACC_RDONLY);

        Epetra_MultiVector* velocity = new Epetra_MultiVector(velocity_map, 3);
        Epetra_MultiVector* pressure = new Epetra_MultiVector(pressure_map, 1);

        HDF5input.Read(velocity_group_name.str(), velocity, true, 0);
        HDF5input.Read(pressure_group_name.str(), pressure, true, 0);

        HDF5input.Close();

        velocity_ref->Update(-1.0, *velocity, 1.0);
        pressure_ref->Update(-1.0, *pressure, 1.0);

        double error_norm_velocity[3];
        double error_norm_pressure;

        velocity_ref->NormInf(error_norm_velocity);
        pressure_ref->NormInf(&error_norm_pressure);

        std::cout << "Velocity error norms: " << error_norm_velocity[0] << " "
                  << error_norm_velocity[1] << " " << error_norm_velocity[2] << std::endl;
        std::cout << "Pressure error norm: " << error_norm_pressure << std::endl;

        for (int k = 0; k < 3; ++k)
        {
            if (error_norm_velocity[k] > tolerance)
            {
                correct = false;
            }
        }
        if (error_norm_pressure > tolerance)
        {
            correct = false;
        }

        velocity_group_name.str(std::string());
        velocity_group_name.clear();
        pressure_group_name.str(std::string());
        pressure_group_name.clear();

        delete velocity;
        delete velocity_ref;
        delete pressure;
        delete pressure_ref;
    }

    return correct;
}
