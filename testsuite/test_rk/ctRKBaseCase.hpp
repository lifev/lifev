#ifndef __CTRK_BASE_CASE_HH
#define __CTRK_BASE_CASE_HH 1

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <ChorinTemamRK.hpp>
#include <iostream>

using namespace LifeV;

/*!
 * \struct CTRKcaseBase
 * \brief Interface for Chorin-Temam with RK2 time stepping  case study
 */

struct CTRKcaseBase
{
  public:
    // typedefs
    typedef boost::function<Real (const Real&, const Real&, const Real&, const Real&, const ID&)> func_type;
    typedef boost::shared_ptr<BCHandler> bc_ptrtype;

    // no ctor
    virtual ~CTRKcaseBase() {}
    
    // methods called by CTRK upstream class
    GetPot &get_data_hdl() {return C_data;}
    bc_ptrtype get_bcHu() {return C_bcHu;}
    bc_ptrtype get_bcHp() {return C_bcHp;}
    void set_base_data(const GetPot& _data, Epetra_Comm *_comm) 
    {
        C_comm = _comm;
	C_data = _data;
	C_hasdata = true;
    }

    // handy die
    void die(std::string msg)
    { 
    	if (!C_comm->MyPID()) {
	  std::cout << "  --  ERROR: " << msg << std::endl;
	  std::cout << "  --  Exiting ..." << std::endl;
	}
	exit(1);
    }

    // Use data interface
    virtual void set_user_data()
    {
        die("Cannot set data without user input.");
    }

    // BC allocation w.r.t. their number
    void create_bcs()
    {
        C_bcHu.reset(new BCHandler(C_num_bcs, BCHandler::HINT_BC_NONE));
	C_bcHp.reset(new BCHandler(C_num_bcs, BCHandler::HINT_BC_NONE));
    }

    // User bcs interface
    virtual void set_bcs()
    {
        die("Cannot set BCs without user input.");
    }
 
  protected:
    // mandatory members
    GetPot C_data;
    bool C_hasdata;
    int C_num_bcs;
    bc_ptrtype C_bcHu;
    bc_ptrtype C_bcHp;
    Epetra_Comm *C_comm;
   
}; // struct CTRKcaseBase

#endif /* __CTRK_BASE_CASE_HH */
