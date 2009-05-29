#ifndef __CT_BASE_CASE_HH
#define __CT_BASE_CASE_HH 1

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <ChorinTemam.hpp>
#include <iostream>

using namespace LifeV;

/*!
 * \struct CTcaseBase
 * \brief Interface for Chorin-Temam / projection case study
 */

struct CTcaseBase
{
  public:
    // typedefs
    typedef boost::function<Real (const Real&, const Real&, const Real&, const Real&, const ID&)> func_type;
    typedef boost::shared_ptr<BCHandler> bc_ptrtype;
  
    // no ctor
    virtual ~CTcaseBase() {}

    // methods called by CT upstream class
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
    
    // User data interface
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
	die("Cannot set BCs witout user input.");
    }

  protected:
    GetPot C_data;
    bool C_hasdata;
    int C_num_bcs;
    bc_ptrtype C_bcHu;
    bc_ptrtype C_bcHp;
    Epetra_Comm *C_comm;

}; // struct CTcaseBase

#endif /* __CT_BASE_CASE_HH */
