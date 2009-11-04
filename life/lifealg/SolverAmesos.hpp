/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file SolverAmesos.hpp

*/



#ifndef _SolverAmesos_H
#define _SolverAmesos_H

#include <Amesos.h>
#include <Amesos_BaseSolver.h>


#include <Teuchos_ParameterList.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifearray/EpetraMatrix.hpp>

#include <life/lifecore/chrono.hpp>
#include <life/lifecore/displayer.hpp>




class GetPot;

namespace LifeV
{
class SolverAmesos
{
public:


        /** @name Typedefs
     */
    //@{

    typedef double                           value_type;

    typedef SolverAmesos                     solver_type;

    typedef EpetraMatrix<double>             matrix_type;
    typedef EpetraVector                     vector_type;

    typedef void                             prec_raw_type;
    typedef void                             prec_type;
    typedef boost::shared_ptr<matrix_type>   matrix_ptrtype;
    typedef boost::shared_ptr<EpetraVector>  vector_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    //! @param filename GetPot data file containing options for solver in
    //!        section "aztec"
    SolverAmesos(Epetra_Comm&);

    //SolverAmesos();
    //@}

    /** @name Accessors
     */
    //@{

    //! Total number of iterations
    int    NumIters();

    //! True Residual
    double TrueResidual();

    //@}

    /** @name  Mutators
     */
    //@{

    //! set matrix from EpetraMatrix
    void setMatrix(matrix_type& m);

    void setOperator(Epetra_Operator& op);

    //! set Epetra_Operator preconditioner
    //void setPreconditioner( prec_type _prec );

    void setDataFromGetPot( const GetPot& dfile, const std::string& section );
    void SetParameters(bool cerr_warning_if_unused = false);

    void setTolMaxiter(const double tol, const int maxiter=-1);


    //! describes the verobisty level
    enum VerboseLevel { NONE, SUMMARY, LAST };
    //! sets verbosity level
    void SetVerbose(const VerboseLevel verb=NONE);


    //@}


    /** @name  Methods
     */
    //@{

    /*
      solve the problem \f$ A x = b \f$

      \c A has been entered via \c setMatrix .

      return the number of iterations, M_maxIter+1 if solve failed

    */
    //int solve( vector_type& x, vector_type& b );

    double computeResidual( vector_type& __X, vector_type& __B );

    bool precSet() const {return true;}
    void precReset() { return; }


    // return the Aztec status


    void printStatus();

//     void getStatus( double status[AZ_STATUS_SIZE])
//     { M_solver.GetAllAztecStatus( status );}

    //@}

    /** Solves the system and returns the number of iterations.
        @param  matrFull,
        @param  rhsFull,
        @param  sol,
        @param  prec,
        @param  reuse,
        @param  retry = true

        returns number of iterations. If negative, the solver did not converge,
        the preconditionar has been recomputed, and a second solution is tried
    */
    int solveSystem(  vector_type&     rhsFull,
                      vector_type&     sol,
                      matrix_ptrtype&  prec,
                      bool const       reuse,
                      bool const       retry=true);

    void setUpPrec    (const GetPot& dataFile,  const std::string& section);
//     void setPrec      (prec_raw_type* prec);

    //prec_type& getPrec(){return 0;}

private:

    matrix_type::matrix_ptrtype  M_matrix;

    Amesos                 M_factory;

    Epetra_LinearProblem   M_problem;

    std::string            M_solverType;
    Teuchos::ParameterList M_TrilinosParameterList;

    // Teuchos::ParameterList M_timingsList;
    // timing variable

    double                 M_sfact_time;
    double                 M_nfact_time;
    double                 M_solve_time;

    //

    double                 M_mtx_conv_time;
    double                 M_mtx_redist_time;
    double                 M_vec_redist_time;

    //

    bool                   M_redistribute;
    bool                   M_printTiming;
    bool                   M_printStatus;

    //

    int                    M_maxIter;
    double                 M_tol;
    int                    M_maxIterSolver;
    int                    M_maxIterForReuse;

    Displayer              M_Displayer;
    Epetra_Comm&           M_comm;



};

} // namespace LifeV

#endif /* __SolverAmesos_H */
