/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2004-09-22

  Copyright (C) 2004 EPFL

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
   \file SolverAztec.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-09-22
*/

#include <SolverAztec.hpp>

class GetPot;

namespace LifeV {

UInt SolverAztec::_solverNumber = 100;

SolverAztec::SolverAztec()
    : _matrix(0), _precond(0), _tempPattern(0), _tempMatrix(0)
{
    _data_org[AZ_N_internal] = 0;
    _data_org[AZ_N_border]= 0;
    _data_org[AZ_N_external]= 0;
    _data_org[AZ_N_neigh]= 0;
    _data_org[AZ_name]= _solverNumber++;
    AZ_set_proc_config(_proc_config, AZ_NOT_MPI);

    // let dataAztec set the defaults
    GetPot dataFile;
    DataAztec dataAztec(dataFile, "aztec");
    dataAztec.aztecOptionsFromDataFile(_options, _params);
}

SolverAztec::~SolverAztec() {
    if (_matrix) {
        AZ_matrix_destroy(&_matrix);
    }
    if (_precond) {
        AZ_precond_destroy(&_precond);
    }
}

SolverAztec* SolverAztec::New() {
    return new SolverAztec;
}

/*!
  \brief Gets the last residual norm that has been computed.

  \return last residual norm
*/
double SolverAztec::residualNorm() const {
    return _status[AZ_r];
}

void SolverAztec::setMatrix(MSRMatr<value_type> const& m) {
    _tempPattern.reset(0);
    _tempMatrix.reset(0);
    _setMatrix(m);
}

void SolverAztec::setMatrix(CSRMatr<CSRPatt, value_type> const& m) {
    _tempPattern.reset(new MSRPatt(*(m.Patt())));
    _tempMatrix.reset(new MSRMatr<value_type>(*_tempPattern, m));
    _setMatrix(*_tempMatrix);
}

void SolverAztec::
setMatrixFree(int nEq, void* data,
              void (*matvec)(double*,double*, AZ_MATRIX_STRUCT*, int*)) {
    _data_org[AZ_N_internal] = nEq;
    if (_matrix) {
        AZ_matrix_destroy(&_matrix);
    }
    if (_precond) {
        AZ_precond_destroy(&_precond);
    }
    _matrix = AZ_matrix_create(nEq);
    AZ_set_MATFREE(_matrix, data, matvec);
}

void SolverAztec::_setMatrix(MSRMatr<value_type> const& m) {
    int nEq = m.Patt()->nRows();
    _data_org[AZ_N_internal] = nEq;
    if (_matrix) {
        AZ_matrix_destroy(&_matrix);
    }
    if (_precond) {
        AZ_precond_destroy(&_precond);
    }
    _matrix = AZ_matrix_create(nEq);
    _precond = AZ_precond_create(_matrix, AZ_precondition, NULL);
    AZ_set_MSR(_matrix, (int*)(m.Patt()->giveRaw_bindx()),
               (double*)(m.giveRaw_value()), _data_org, 0, NULL, AZ_LOCAL);
}

void SolverAztec::solve( array_type& x, array_type const& b) {
    if (_matrix) {
        std::cerr << "[SolverAztec::solve]  Solving primal\n";
        AZ_iterate(&x[0], &b[0], _options, _params, _status, _proc_config,
                   _matrix, _precond, NULL);
    } else {
        std::cerr << "[SolverAztec::solve]  ERROR: Matrix not set\n";
    }
}

void SolverAztec::setOptionsFromGetPot(GetPot const& dataFile,
                                       std::string section) {
    DataAztec dataAztec(dataFile, section);
    dataAztec.aztecOptionsFromDataFile(_options, _params);
}

} // namespace LifeV
