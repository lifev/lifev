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
     @brief This file contains the definition of the flags for the ETCurrentFE class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef ETCURRENTFEFLAG_HPP
#define ETCURRENTFEFLAG_HPP

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

/*! \page Update_flag Flags to update the ETCurrentFE

  \section flag_goal Why do we need flags?

  Flags are needed in order to save computation. Indeed, when the LifeV::ETCurrentFE structure has
  to be updated for a given element, it is possible that all the informations (e.g. second derivatives
  of the basis functions) are not needed. The flag carry the information of the values to be computed
  because they are required.

  \section flag_definition What is a flag?

  At first sight, we might think that a flag is a complicated class implemented to do exaclty what we want,
  with overloaded operators... Actually, it is much simpler: a LifeV::flag_Type is just an unsigned integer.

  \section flag_primitive How to define a flag?

  The flags use the binary representation of the integers to work. This enables a very fast definition and
  use of the flags. To understand it, let us make a simple example. Suppose that we can update three
  quantities A,B and C.

  The first step is to define a "primitive" flag for each of these quantities. These flags are defined as
  powers of 2. Here, we will define

  \code
  flag_Type UPDATE_A(1);
  flag_Type UPDATE_B(2);
  flag_Type UPDATE_C(4);
  \endcode

  We need here powers of 2 because this makes the binary representation of the "primitive" flags simple:
  UPDATE_A is 001, UPDATE_B is 010 and UPDATE_C is 100 (if the integers are coded on 3 bits, otherwise there
  are many zeros before). The fact that we want A to be update is then represented by a 1 in the third position,
  for B it is the second position and for C the first position.

  So now, if we want to build a flag that updates A and C, we want it to be in binary 101, that is 5. So,
  the flag UPDATE_AC will be defined as:

  \code
  flag_Type UPDATE_AC(5);
  \endcode

  \section flag_combine How to combine flags?

  With the last example, one could think that is it just a matter of addition, but it is not. Suppose that we
  want to combine the flags UPDATE_A and UPDATE_AC. If we add the corresponding values, we will get 1+5=6 that
  is represented in binary by 110. This would mean update B and C, what is not what we want!

  To combine flags, we use the binary operator | (bitwise "OR" operation) that do exactly the job that we want: 1|5=5.

  \section flag_detect How to detect flags?

  Now that we can build every flag that we want, we want to detect quickly the different flags. This is achievied
  by using the binary operator & (bitwise "AND" operation) and the "primitive" flags.

  Suppose that we want to know if A has to be updated. Then, we perform the operation "& UPDATE_A" on the incoming
  flag. If the result is zero, then we do not need to update it: for example, UPDATE_B & UPDATE_A = 0 . Otherwise,
  we have to update it: for example, UPDATE_AC & UPDATE_A = 1.


  \section flag_list Possible flags

  There are several flags defined for you. Here is a list of the possible flags that are usually used:

  <table>
   <tr>
    <th> Flag name </th> <th> Effect </th>
   </tr>
   <tr>
    <th> ET_UPDATE_QUAD_NODES </th> <th> Update everything needed to know the position of the quadrature nodes in the current cell </th>
   </tr>
   <tr>
    <th> ET_UPDATE_DPHI </th> <th>  Update everything needed to know the values of the derivatives of the basis functions in the quadrature nodes (in the current cell) </th>
   </tr>
   <tr>
    <th> ET_UPDATE_D2PHI </th> <th>  Update everything needed to know the values of the second derivatives of the basis functions in the quadrature nodes (in the current cell) </th>
   </tr>
   <tr>
    <th> ET_UPDATE_WDET </th> <th>  Update everything needed to know the determinant of the transformation multiplied by the weights of the quadrature. </th>
   </tr>

  </table>

  Note: in the old versions there was also the flag UPDATE_PHI. This flag has been removed, since the values of the basis functions in the quadrature nodes are always the same, so they do not need to be updated.
  Besides this usual flags, there are a couple of "primitive" flags, that update only a particular element in the currentFE structure. Be sure to know what you are doing before using them.

*/

// For more informations about the flags, please visite the documentation \ref Update_flag

/*! Typedef for the flag_Type */
typedef unsigned int flag_Type;

// PRIMITIVE FLAGS
// The flags containing "ONLY" are not intended to be directly used.
// They are rather meant to be composed into other flags that take
// care for the completness of the procedure.

// Do nothing flag
const flag_Type ET_UPDATE_NONE (0);

// Update cell coordinates
const flag_Type ET_UPDATE_ONLY_CELL_NODE (1);

// Update quadrature points coordinate in the current element
const flag_Type ET_UPDATE_ONLY_QUAD_NODE (2);

//Update the jacobian of the transformation
const flag_Type ET_UPDATE_ONLY_JACOBIAN (4);

// Update the determinant of the jacobian
const flag_Type ET_UPDATE_ONLY_DET_JACOBIAN (8);

// Update the inverse of the jacobian
const flag_Type ET_UPDATE_ONLY_T_INVERSE_JACOBIAN (16);

// Update the weighted determinant only
const flag_Type ET_UPDATE_ONLY_W_DET_JACOBIAN (32);

// Update the derivative of the basis functions
const flag_Type ET_UPDATE_ONLY_DPHI (64);

// Update the second derivative of the basis functions
const flag_Type ET_UPDATE_ONLY_D2PHI (128);

// Update the divergence of the basis functions
const flag_Type ET_UPDATE_ONLY_DIVERGENCE (256);

// Update the laplacian of the basis functions
const flag_Type ET_UPDATE_ONLY_LAPLACIAN (512);

// Update the diameter of the triangle
const flag_Type ET_UPDATE_ONLY_DIAMETER (1024);

// Update the measure of the triangle
const flag_Type ET_UPDATE_ONLY_MEASURE (2048);

// Update everything
const flag_Type ET_UPDATE_ALL (4096 - 1);


// COMPOSITE FLAGS
// These flags are composed of elementary flags. They take care
// that all the intermediate quantities required to compute
// the final (required) quantity are updated as well.

// Flag for the quadrature nodes in the current cell
const flag_Type ET_UPDATE_QUAD_NODE (ET_UPDATE_ONLY_CELL_NODE
                                     | ET_UPDATE_ONLY_QUAD_NODE);

// Flag for the gradient of the basis functions
const flag_Type ET_UPDATE_DPHI (ET_UPDATE_ONLY_CELL_NODE
                                | ET_UPDATE_ONLY_JACOBIAN
                                | ET_UPDATE_ONLY_DET_JACOBIAN
                                | ET_UPDATE_ONLY_T_INVERSE_JACOBIAN
                                | ET_UPDATE_ONLY_DPHI);

// Flag for the second derivate of the basis functions
const flag_Type ET_UPDATE_D2PHI (ET_UPDATE_ONLY_CELL_NODE
                                 | ET_UPDATE_ONLY_JACOBIAN
                                 | ET_UPDATE_ONLY_DET_JACOBIAN
                                 | ET_UPDATE_ONLY_T_INVERSE_JACOBIAN
                                 | ET_UPDATE_ONLY_D2PHI);

// Flag for the weighted determinant
const flag_Type ET_UPDATE_WDET (ET_UPDATE_ONLY_CELL_NODE
                                | ET_UPDATE_ONLY_JACOBIAN
                                | ET_UPDATE_ONLY_DET_JACOBIAN
                                | ET_UPDATE_ONLY_W_DET_JACOBIAN);

// Flag for the gradient of the basis functions
const flag_Type ET_UPDATE_DIVERGENCE (ET_UPDATE_ONLY_CELL_NODE
                                      | ET_UPDATE_ONLY_JACOBIAN
                                      | ET_UPDATE_ONLY_DET_JACOBIAN
                                      | ET_UPDATE_ONLY_T_INVERSE_JACOBIAN
                                      | ET_UPDATE_ONLY_DPHI
                                      | ET_UPDATE_ONLY_DIVERGENCE);

// Flag for the laplacian of the basis functions
const flag_Type ET_UPDATE_LAPLACIAN  (ET_UPDATE_ONLY_CELL_NODE
                                      | ET_UPDATE_ONLY_JACOBIAN
                                      | ET_UPDATE_ONLY_DET_JACOBIAN
                                      | ET_UPDATE_ONLY_T_INVERSE_JACOBIAN
                                      | ET_UPDATE_ONLY_D2PHI
                                      | ET_UPDATE_ONLY_LAPLACIAN);

// Flag for the diameter of the cell
const flag_Type ET_UPDATE_DIAMETER (ET_UPDATE_ONLY_CELL_NODE
                                    | ET_UPDATE_ONLY_DIAMETER);

// Flag for the diameter of the cell
const flag_Type ET_UPDATE_MEASURE (ET_UPDATE_WDET
                                   | ET_UPDATE_ONLY_MEASURE);



} // Namespace LifeV

#endif // ETCURRENTFEFLAG_HPP
