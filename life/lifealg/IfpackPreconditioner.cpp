/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2006-11-09

  Copyright (C) 2006 EPFL

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
   \file EpetraPreconditioner.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2006-11-09
 */

#include "IfpackPreconditioner.hpp"



namespace LifeV
{


IfpackPreconditioner::IfpackPreconditioner():
        super (),
        M_Prec()
{
}

IfpackPreconditioner::~IfpackPreconditioner()
{}

// IfpackPreconditioner::IfpackPreconditioner(operator_type& oper):
//         super(oper),
//         M_Prec()
// {
//     buildPreconditioner(oper);
// }



void
IfpackPreconditioner::setDataFromGetPot( const GetPot& dataFile,
                                         const std::string& section )
{

    //! See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    //! for more informations on the parameters

    M_overlapLevel = dataFile((section + "/ifpack/overlap").data(),     4);
    M_precType     = dataFile((section + "/ifpack/prectype").data(),"Amesos");

    //Teuchos::ParameterList list;

    createIfpackList(dataFile, section, M_List);

    //this->setList(M_List);


}

int
IfpackPreconditioner::buildPreconditioner(operator_type& oper)
{
    M_Oper = oper;

    //List.set("schwarz: combine mode", "Zero"); //
    //M_List.set("schwarz: filter singletons", true);

//    List.set("amesos: solver type", "Amesos_Lapack");

//     M_List.set("PrintTiming", true);
//     M_List.set("PrintStatus", true);

    Ifpack factory;

    //    M_Prec

    M_Prec.reset(factory.Create(M_precType, &M_Oper->getEpetraMatrix(), M_overlapLevel));
    //    M_Prec.reset(new prec_type(&A.getEpetraMatrix(), OverlapLevel));
    if ( !M_Prec.get() )
        { //! if not filled, I do not know how to diagonalize.
            ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
        }
    IFPACK_CHK_ERR(M_Prec->SetParameters(M_List));
    IFPACK_CHK_ERR(M_Prec->Initialize());
    IFPACK_CHK_ERR(M_Prec->Compute());
    return EXIT_SUCCESS;
}


// void
// IfpackPreconditioner::createList( const GetPot&              dataFile,
//                                   const std::string&         section,
//                                   Teuchos::ParameterList&    list)
// {
//     //! See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
//     //! for more informations on the parameters

//     bool displayList = dataFile((section + "/displayList").data(),     false);

//     std::string relaxationType              = dataFile((section + "/ifpack/relaxation/type").data(), "Jacobi");
//     int         relaxationSweeps            = dataFile((section + "/ifpack/relaxation/sweeps").data(), 1);
//     double      relaxationDampingFactor     = dataFile((section + "/ifpack/relaxation/damping_factor").data(), 1.0);
//     double      relaxationMinDiagValue      = dataFile((section + "/ifpack/relaxation/min_diagonal_value").data(), 0.0);
//     bool        relaxationZeroStartSolution = dataFile((section + "/ifpack/relaxation/zero_starting_solution").data(), true);

//     list.set("relaxation: type",                    relaxationType);
//     list.set("relaxation: sweeps",                  relaxationSweeps);
//     list.set("relaxation: damping factor",          relaxationDampingFactor);
//     list.set("relaxation: min diagonal value",      relaxationMinDiagValue);
//     list.set("relaxation: zero starting solution",  relaxationZeroStartSolution);


//     std::string partitionerType             = dataFile((section + "/ifpack/partitioner/type").data(), "metis");
//     int         partitionerOverlap          = dataFile((section + "/ifpack/partitioner/overlap").data(), 0);
//     int         partitionerLocalParts       = dataFile((section + "/ifpack/partitioner/local_parts").data(), 1);
//     int         partitionerRootNode         = dataFile((section + "/ifpack/partitioner/root_node").data(), 0);
//     bool        partitionerUseSymmGraph     = dataFile((section + "/ifpack/partitioner/use_symmetric_graph").data(), true);

//     list.set("partitioner: type",                partitionerType);
//     list.set("partitioner: overlap",             partitionerOverlap);
//     list.set("partitioner: local parts",         partitionerLocalParts);
//     list.set("partitioner: root node",           partitionerRootNode);
//     list.set("partitioner: use symmetric graph", partitionerUseSymmGraph);


//     std::string amesosSolverType            = dataFile((section + "/ifpack/amesos/solvertype").data(), "Amesos_KLU");

//     list.set("amesos: solver type", amesosSolverType);

//     double levelOfFill     = dataFile((section + "/ifpack/fact/level-of-fill").data(),      4.);
//     double ILUTlevelOfFill = dataFile((section + "/ifpack/fact/ilut_level-of-fill").data(), 4.);
//     double athr            = dataFile((section + "/ifpack/fact/absolute_threshold").data(), 0.);
//     double rthr            = dataFile((section + "/ifpack/fact/relative_threshold").data(), 1.);
//     double relaxValue      = dataFile((section + "/ifpack/fact/relax_value").data(),        0.);
//     double dropTolerance   = dataFile((section + "/ifpack/fact/drop_tolerance").data(),     1e-5);


//     list.set("fact: drop tolerance",     dropTolerance);
//     list.set("fact: level-of-fill",      levelOfFill);
//     list.set("fact: ilut level-of-fill", ILUTlevelOfFill);
//     list.set("fact: absolute threshold", athr);
//     list.set("fact: relative threshold", rthr);
//     list.set("fact: relax value",        relaxValue);

//     int combineMode              = dataFile((section + "/ifpack/schwarz/combine_mode").data(),    0);
//     Epetra_CombineMode schwarzCombineMode;

//      switch(combineMode)
//          {
//          case 0 :
//              schwarzCombineMode = Add;
//              break;
//          case 1 :
//              schwarzCombineMode = Zero;
//              break;
//          case 2 :
//              schwarzCombineMode = Insert;
//              break;
//          case 3 :
//              schwarzCombineMode = Average;
//              break;
//          case 4 :
//              schwarzCombineMode = AbsMax;
//              break;
//          default:
//              schwarzCombineMode = Zero;
//          }


//     bool schwarzComputeCondest          = dataFile((section + "/ifpack/schwarz/compute_condest").data(),   true);
//     std::string schwarzReorderingType   = dataFile((section + "/ifpack/schwarz/reordering_type").data(),      "none");
//     bool schwarzFilterSingletons        = dataFile((section + "/ifpack/schwarz/filter_singletons").data(), false);

//     list.set("schwarz: combine mode",       schwarzCombineMode);
//     list.set("schwarz: compute condest",    schwarzComputeCondest);
//     list.set("schwarz: reordering type",    schwarzReorderingType);
//     list.set("schwarz: filter singletons",  schwarzFilterSingletons);

//     if (displayList) list.print(std::cout);
// }


double
IfpackPreconditioner::Condest()
{
    return M_Prec->Condest();
}


EpetraPreconditioner::prec_raw_type*
IfpackPreconditioner::getPrec()
{
    return M_Prec.get();
}

void
IfpackPreconditioner::precReset()
{
    M_Oper.reset();
    M_Prec.reset();
}



void
createIfpackList( const GetPot&              dataFile,
                  const std::string&         section,
                  Teuchos::ParameterList&    list)
{
    //! See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    //! for more informations on the parameters

    bool displayList = dataFile((section + "/displayList").data(),     false);


    std::string relaxationType              = dataFile((section + "/ifpack/relaxation/type").data(), "Jacobi");
    int         relaxationSweeps            = dataFile((section + "/ifpack/relaxation/sweeps").data(), 1);
    double      relaxationDampingFactor     = dataFile((section + "/ifpack/relaxation/damping_factor").data(), 1.0);
    double      relaxationMinDiagValue      = dataFile((section + "/ifpack/relaxation/min_diagonal_value").data(), 0.0);
    bool        relaxationZeroStartSolution = dataFile((section + "/ifpack/relaxation/zero_starting_solution").data(), true);

    list.set("relaxation: type",                    relaxationType);
    list.set("relaxation: sweeps",                  relaxationSweeps);
    list.set("relaxation: damping factor",          relaxationDampingFactor);
    list.set("relaxation: min diagonal value",      relaxationMinDiagValue);
    list.set("relaxation: zero starting solution",  relaxationZeroStartSolution);


    std::string partitionerType             = dataFile((section + "/ifpack/partitioner/type").data(), "metis");
    int         partitionerOverlap          = dataFile((section + "/ifpack/partitioner/overlap").data(), 0);
    int         partitionerLocalParts       = dataFile((section + "/ifpack/partitioner/local_parts").data(), 1);
    int         partitionerRootNode         = dataFile((section + "/ifpack/partitioner/root_node").data(), 0);
    bool        partitionerUseSymmGraph     = dataFile((section + "/ifpack/partitioner/use_symmetric_graph").data(), true);

    list.set("partitioner: type",                partitionerType);
    list.set("partitioner: overlap",             partitionerOverlap);
    list.set("partitioner: local parts",         partitionerLocalParts);
    list.set("partitioner: root node",           partitionerRootNode);
    list.set("partitioner: use symmetric graph", partitionerUseSymmGraph);


    std::string amesosSolverType            = dataFile((section + "/ifpack/amesos/solvertype").data(), "Amesos_KLU");

    list.set("amesos: solver type", amesosSolverType);

    double levelOfFill     = dataFile((section + "/ifpack/fact/level-of-fill").data(),      4.);
    double ILUTlevelOfFill = dataFile((section + "/ifpack/fact/ilut_level-of-fill").data(), 4.);
    double athr            = dataFile((section + "/ifpack/fact/absolute_threshold").data(), 0.);
    double rthr            = dataFile((section + "/ifpack/fact/relative_threshold").data(), 1.);
    double relaxValue      = dataFile((section + "/ifpack/fact/relax_value").data(),        0.);
    double dropTolerance   = dataFile((section + "/ifpack/fact/drop_tolerance").data(),     1e-5);


    list.set("fact: drop tolerance",     dropTolerance);
    list.set("fact: level-of-fill",      levelOfFill);
    list.set("fact: ilut level-of-fill", ILUTlevelOfFill);
    list.set("fact: absolute threshold", athr);
    list.set("fact: relative threshold", rthr);
    list.set("fact: relax value",        relaxValue);

    int combineMode              = dataFile((section + "/ifpack/schwarz/combine_mode").data(),    0);
    Epetra_CombineMode schwarzCombineMode;

     switch(combineMode)
         {
         case 0 :
             schwarzCombineMode = Add;
             break;
         case 1 :
             schwarzCombineMode = Zero;
             break;
         case 2 :
             schwarzCombineMode = Insert;
             break;
         case 3 :
             schwarzCombineMode = Average;
             break;
         case 4 :
             schwarzCombineMode = AbsMax;
             break;
         default:
             schwarzCombineMode = Zero;
         }


    bool schwarzComputeCondest          = dataFile((section + "/ifpack/schwarz/compute_condest").data(),   true);
    std::string schwarzReorderingType   = dataFile((section + "/ifpack/schwarz/reordering_type").data(),      "none");
    bool schwarzFilterSingletons        = dataFile((section + "/ifpack/schwarz/filter_singletons").data(), true);

    list.set("schwarz: combine mode",       schwarzCombineMode);
    list.set("schwarz: compute condest",    schwarzComputeCondest);
    list.set("schwarz: reordering type",    schwarzReorderingType);
    list.set("schwarz: filter singletons",  schwarzFilterSingletons);

    if (displayList) list.print(std::cout);
}









// void
// IfpackPreconditioner::createList( const GetPot& /*dataFile*/ )
// {
// }


// namespace
// {
// EpetraPreconditioner* createIfpack(){ std::cout << "*******************"<< std::endl;return new IfpackPreconditioner(); }
// static bool reg = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack )
//                    && (std::cout << "*************** SIMONE" << std::endl) ) ;
// }

//} // namespace Epetra

} // namespace LifeV
