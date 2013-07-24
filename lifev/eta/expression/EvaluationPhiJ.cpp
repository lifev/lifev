//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains the definition of the EvaluationPhiJ class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <lifev/eta/expression/EvaluationPhiJ.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

const flag_Type EvaluationPhiJ<1>::S_globalUpdateFlag = ET_UPDATE_NONE;

const flag_Type EvaluationPhiJ<1>::S_testUpdateFlag = ET_UPDATE_NONE;

const flag_Type EvaluationPhiJ<1>::S_solutionUpdateFlag = ET_UPDATE_PHI;

} // Namespace ExpressionAssembly

} // Namespace LifeV
