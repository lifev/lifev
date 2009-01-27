/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

   This file is part of the LifeV library

   Author(s): Gilles Fourestey <gilles.fourestey@imag.fr>
   Simone Deparis   <simone.deparis@epfl.ch>
   Date: 2007-06-01

   Copyright (C) 2008 EPFL

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

   Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file MLTester.hpp
 */


#ifndef _MLTESTER_H_
#define _MLTESTER_H_

#include <life/lifesolver/Oseen.hpp>


namespace LifeV
{
/*
 */

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class MLTester : public Oseen<Mesh, SolverType>
{

public :

    typedef Oseen<Mesh, SolverType> super;

    typedef DataNavierStokes<Mesh> data_type;

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> source_type;

    typedef Mesh mesh_type;

    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;

//    typedef SolverType                    solver_type;
    typedef typename super::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>        matrix_ptrtype;
    typedef typename super::vector_type      vector_type;

    typedef typename super::prec_raw_type    prec_raw_type;
    typedef typename super::prec_type        prec_type;


    // constructor
    MLTester( const data_type&          dataType,
              FESpace<Mesh, EpetraMap>& uFESpace,
              FESpace<Mesh, EpetraMap>& pFESpace,
              Epetra_Comm&              comm );

    MLTester( const data_type&          dataType,
              FESpace<Mesh, EpetraMap>& uFESpace,
              FESpace<Mesh, EpetraMap>& pFESpace,
              Epetra_Comm&              comm,
              const EpetraMap           monolithicMap,
              const UInt                offset=0);



    // destructor
    ~MLTester(){};


    //

    void createDefaultList( const GetPot&       dataFile,
                            const std::string&  section );

    void testML(  bchandler_raw_type& bch );

private :

    Teuchos::ParameterList          M_MLDefaultList;


};


template< typename Mesh, typename SolverType >
MLTester<Mesh, SolverType>::
MLTester( const data_type&          dataType,
          FESpace<Mesh, EpetraMap>& uFESpace,
          FESpace<Mesh, EpetraMap>& pFESpace,
          Epetra_Comm&              comm ):
    super(dataType, uFESpace, pFESpace, comm)
{};



template< typename Mesh, typename SolverType >
MLTester<Mesh, SolverType>::
MLTester( const data_type&          dataType,
          FESpace<Mesh, EpetraMap>& uFESpace,
          FESpace<Mesh, EpetraMap>& pFESpace,
          Epetra_Comm&              comm,
          const EpetraMap           monolithicMap,
          const UInt                offset):
    super(dataType, uFESpace, pFESpace, comm, monolithicMap, offset)
{};




template < typename Mesh, typename SolverType >
void
MLTester<Mesh, SolverType>::createDefaultList( const GetPot&       dataFile,
                                               const std::string&  section )
{
    std::cout << "Section = " << section << std::endl;
    createMLList(dataFile, section, M_MLDefaultList);
};


template < typename Mesh, typename SolverType >
void MLTester<Mesh, SolverType>::
testML( bchandler_raw_type& bch )
{


    Chrono chrono;


    // matrix and vector assembling communication

    this->leaderPrint("  f-  Finalizing the matrix and vectors ...    ");

    chrono.start();


    this->M_matrNoBC->GlobalAssemble();

    if (this->M_stab)
        this->M_matrStab->GlobalAssemble();

    matrix_ptrtype matrFull( new matrix_type( this->M_localMap, this->M_matrNoBC->getMeanNumEntries()));

    updateStab(*matrFull);
    getFluidMatrix(*matrFull);

    vector_type rhsFull (this->M_rhsNoBC);

//     matrFull.reset(new matrix_type(*M_matrNoBC));
//     M_rhsFull = M_rhsNoBC;

    chrono.stop();

    this->leaderPrintMax("done in ", chrono.diff() );

    // boundary conditions update
    this->M_comm->Barrier();

    this->leaderPrint("  f-  Applying boundary conditions ...         ");

    chrono.start();
    applyBoundaryConditions( *matrFull, rhsFull, bch);

    chrono.stop();

    this->leaderPrintMax("done in " , chrono.diff());

    // solving the system

    //    Epetra_RowMatrix* A = this->M_linearSolver.GetUserMatrix();

    //this->SetParameters(true);

    bool TestPassed = true;


    // parameter list definition

    vector<string> TestList;

    TestList.push_back("IFPACK");
//     TestList.push_back("Aztec");
//     TestList.push_back("Jacobi");
//     TestList.push_back("Chebyshev");

    vector<std::string> DefaultParamList;


    DefaultParamList.push_back("ML");
    //    DefaultParamList.push_back("SA");
//     DefaultParamList.push_back("DD");
//     DefaultParamList.push_back("DD-ML");
//     DefaultParamList.push_back("NSSA");


    vector<std::string> CoarseParamList;

    CoarseParamList.push_back("Aztec");
    CoarseParamList.push_back("IFPACK");
    CoarseParamList.push_back("Jacobi");
//     CoarseParamList.push_back("ML symmetric Gauss-Seidel");
    CoarseParamList.push_back("symmetric Gauss-Seidel");
    //    CoarseParamList.push_back("ML Gauss-Seidel");
    CoarseParamList.push_back("Gauss-Seidel");
    CoarseParamList.push_back("Chebyshev");
    CoarseParamList.push_back("MLS");
    //    CoarseParamList.push_back("Hiptmair");
    CoarseParamList.push_back("Amesos-KLU");
    CoarseParamList.push_back("Amesos-Superlu");
    CoarseParamList.push_back("Amesos-UMFPACK");
    CoarseParamList.push_back("Amesos-Superludist");
    CoarseParamList.push_back("Amesos-MUMPS");
    //    CoarseParamList.push_back("user-defined");
    CoarseParamList.push_back("SuperLU");
//     CoarseParamList.push_back("IFPACK-Chebyshev");
    CoarseParamList.push_back("self");
    CoarseParamList.push_back("do-nothing");
    CoarseParamList.push_back("IC");
    CoarseParamList.push_back("ICT");
    CoarseParamList.push_back("ILU");
    //CoarseParamList.push_back("ILUT");


    vector<string> PreOrPost;
//     PreOrPost.push_back("pre");
//     PreOrPost.push_back("post");
    PreOrPost.push_back("both");

    vector<double> Damping;
    Damping.push_back(1.00);

    std::string AZstatus = "";

    int    bestIter = 100;
    double bestTime = 100.;

    std::string bestTimePoP      = "";
    std::string bestIterPoP      = "";


    std::string bestTimeSmoother = "";
    std::string bestIterSmoother = "";

    std::string bestTimeDefParamList = "";
    std::string bestIterDefParamList = "";


    double bestTimeDamping  = 0.;
    double bestIterDamping  = 0.;

    double bestTimeSweeps  = 0.;
    double bestIterSweeps  = 0.;

    std::ofstream ofile( "results.txt" );

    for (int smi = 0; smi < DefaultParamList.size(); smi++)
    for (int cl = 0; cl < CoarseParamList.size(); cl++)
        for (int sweeps = 1 ; sweeps < 2 ; sweeps += 2)
            for (unsigned int i = 0; i < TestList.size(); ++i)
                for (unsigned int j = 0; j < PreOrPost.size(); ++j)
                    for (unsigned int k = 0; k < Damping.size(); ++k)
                        {
                            if (this->M_comm->MyPID() == 0)
                                {
                                    std::cout << std::endl;
                                    //                                    std::cout << "### Testing       " << DefaultParamList[smi];
                                    std::cout << "### Testing       " << CoarseParamList[cl];
                                    std::cout << " with " << TestList[i] << " ";
                                    std::cout << "### sweeps      = " << sweeps << " ";
                                    std::cout << "### pre or post = " << PreOrPost[j] << " ";
                                    std::cout << "### damping     = " << Damping[k] << std::endl;

                                    ofile << std::endl;
                                    //                                    ofile << "### Testing       " << DefaultParamList[smi];
                                    ofile << "### Testing       " << CoarseParamList[cl];
                                    ofile << " with " << TestList[i] << " ";
                                    ofile << "### sweeps      = " << sweeps << " ";
                                    ofile << "### pre or post = " << PreOrPost[j] << " ";
                                    ofile << "### damping     = " << Damping[k] << std::endl;
                                }

                            // ========================= //
                            // build ML with ML smoother //
                            // ========================= //

                            Teuchos::ParameterList MLList = M_MLDefaultList;

                            //ML_Epetra::SetDefaults(DefaultParamList[smi],MLList);

                            string mlSmootherType;
                            if (TestList[i] == "Gauss-Seidel")
                                mlSmootherType = "ML Gauss-Seidel";
                            else if (TestList[i] == "symmetric Gauss-Seidel")
                                mlSmootherType = "ML symmetric Gauss-Seidel";
                            else if (TestList[i] == "Chebyshev")
                                {
                                    mlSmootherType = TestList[i];
                                    MLList.set("smoother: polynomial order", sweeps);
                                }
                            else
                                mlSmootherType = TestList[i];

                            MLList.set("smoother: type", mlSmootherType);
                            MLList.set("smoother: pre or post", PreOrPost[j]);
                            //                            MLList.set("smoother: sweeps", sweeps);
                            MLList.set("smoother: damping factor", Damping[k]);

//                             MLList.set("smoother: type (level 0)","IFPACK");
//                             MLList.set("smoother: type (level 1)","IFPACK");
//                             MLList.set("smoother: type (level 2)","IFPACK");

                            MLList.set("coarse: type", CoarseParamList[cl]);

                            MLList.set("ML output", 0);

                            //MLList.print(std::cout);

                            this->M_prec.reset( PRECFactory::instance().createObject( "ML" ) );
                            this->M_prec->setList(MLList);

                            Epetra_Time Time(*this->M_comm);

                            this->M_sol.getEpetraVector().Scale(0.);

                            this->solveSystem(matrFull, rhsFull, this->M_sol, this->M_linearSolver, this->M_prec);


                            this->resetPrec();

                            const double* status =  this->M_linearSolver.getAztecStatus();

                            if( status[AZ_why] == AZ_normal         ) AZstatus = "N";
                            else if( status[AZ_why] == AZ_maxits    ) AZstatus = "M";
                            else if( status[AZ_why] == AZ_loss      ) AZstatus = "L";
                            else if( status[AZ_why] == AZ_ill_cond  ) AZstatus = "I";
                            else if( status[AZ_why] == AZ_breakdown ) AZstatus = "B";

                            int MLIters   = this->M_linearSolver.NumIters();
                            double MLTime = Time.ElapsedTime();

                            if (AZstatus == "N")
                                {
                                    if (MLIters < bestIter)
                                        {
                                            bestIter          = MLIters;
                                            bestIterPoP       = PreOrPost[j];
                                            bestIterSmoother  = TestList[i];
                                            bestIterDefParamList  = DefaultParamList[smi];
                                            bestIterDamping   = Damping[k];
                                            bestIterSweeps    = sweeps;
                                        }

                                    if (MLTime < bestTime)
                                        {
                                            bestTime          = MLTime;
                                            bestTimePoP       = PreOrPost[j];
                                            bestTimeSmoother  = TestList[i];
                                            bestTimeDefParamList  = DefaultParamList[smi];
                                            bestTimeDamping   = Damping[k];
                                            bestTimeSweeps    = sweeps;
                                        }
                                }


                            if (this->M_comm->MyPID() == 0)
                                {
                                    std::cout << "ML     : iters     = " << MLIters;
                                    std::cout << " time      = " << MLTime << " (s) with "
                                              << AZstatus << " status, residual = "
                                              << status[AZ_scaled_r] << std::endl;
                                    ofile << "ML     : iters     = " << MLIters;
                                    ofile << " time      = " << MLTime << " (s) with "
                                          << AZstatus << " status, residual = "
                                          << status[AZ_scaled_r] << std::endl;


                                }



                        }


    if (this->M_comm->MyPID() == 0)
        {
            std::cout << "Best time : " << bestTime << std::endl;
            std::cout << "  Def. Param = " << bestTimeDefParamList;

            std::cout << "  Smoother   = " << bestTimeSmoother << std::endl;
            std::cout << "  PoP        = " << bestTimePoP << std::endl;
            std::cout << "  Damping    = " << bestTimeDamping << std::endl;
            std::cout << "  Sweeps     = " << bestTimeSweeps << std::endl;

            std::cout << "Best iter : " << bestIter << std::endl;
            std::cout << "  Def. Param = " << bestIterDefParamList;
            std::cout << "  Smoother   = " << bestIterSmoother << std::endl;
            std::cout << "  PoP        = " << bestIterPoP << std::endl;
            std::cout << "  Damping    = " << bestIterDamping << std::endl;
            std::cout << "  Sweeps     = " << bestIterSweeps << std::endl;


            ofile << "Best time : " << bestTime << std::endl;
            ofile << "  Def. Param = " << bestTimeDefParamList;

            ofile << "  Smoother   = " << bestTimeSmoother << std::endl;
            ofile << "  PoP        = " << bestTimePoP << std::endl;
            ofile << "  Damping    = " << bestTimeDamping << std::endl;
            ofile << "  Sweeps     = " << bestTimeSweeps << std::endl;

            ofile << "Best iter : " << bestIter << std::endl;
            ofile << "  Def. Param = " << bestIterDefParamList;
            ofile << "  Smoother   = " << bestIterSmoother << std::endl;
            ofile << "  PoP        = " << bestIterPoP << std::endl;
            ofile << "  Damping    = " << bestIterDamping << std::endl;
            ofile << "  Sweeps     = " << bestIterSweeps << std::endl;
        }




    ofile.close();

}




}

#endif
