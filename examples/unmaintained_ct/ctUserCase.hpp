#ifndef __CT_USER_CASE_HH
#define __CT_USER_CASE_HH 1

#include <ctBaseCase.hpp>

using namespace LifeV;

/*!
 * \struct CTcaseUser
 * \brief User defined case study. This is the (a priori only) piece of code to modify
 * to account for possible test studies.
 */

struct CTcaseUser
            : public CTcaseBase
{
public:
    // no ctor
    ~CTcaseUser() {}

private:
    // user defined members (those read in the data file)
    Real C_density;
    Real C_viscosity;
    Real C_dimX;
    Real C_dimY;
    Real C_dimZ;
    Real C_diameter;
    Real C_ampl;
    Real C_amplstep;

    // user defined members (eg. reference used in the mesh)
    int INLET;
    int WALL;
    int SLIPWALL;
    int OUTLET;
    int CYLINDER;
    std::vector<ID> zComp;

    // user defined (boost) functions
    func_type u3D_zero;
    func_type u3D_in;

    // user defined methods to be used like functions in BC, etc.

    // zero
    Real u3DZero(const Real&, const Real&, const Real&, const Real&, const ID&)
    {return(0);}

    // time independant parabolic profile
    Real u3Dcyl(const Real&, const Real&, const Real& y, const Real&, const ID& i)
    {
        if (i == 1)
            return ( C_ampl / (C_dimY * C_dimY)*(y + C_dimY)*(C_dimY-y) );
        else
            return(0.);
    }

    // time dependant parabolic profile
    Real u3Dcyl_dyn(const Real& t, const Real&, const Real& y, const Real&, const ID& i)
    {
        if (i==1)
        {
            if (t < C_amplstep)
                return ( (t/C_amplstep) * C_ampl / (C_dimY * C_dimY)*(y + C_dimY)*(C_dimY-y) );
            else
                return ( C_ampl / (C_dimY * C_dimY)*(y + C_dimY)*(C_dimY-y) );
        }
        else
            return(0.);
    }

public:

    // user defined data settings
    void set_user_data()
    {
        std::cout << "  c- Setting user defined data." << std::endl;

        if (!C_hasdata) die("Data case has not be loaded");
        // now read whatever data we need
        C_density   = C_data ("fluid/physics/density", 1.);
        C_viscosity = C_data ("fluid/physics/vicosity", 1.);
        C_dimX = C_data ("fluid/problem/dimX", 40.);
        C_dimY = C_data ("fluid/problem/dimY", 20.);
        C_dimZ = C_data ("fluid/problem/dimZ", 4.);
        C_diameter = C_data ("fluid/problem/diameter", 1.);
        C_ampl = C_data ("fluid/problem/amplitude", 1.);
        C_amplstep = C_data ("fluid/problem/amplstep", 2.);

        // set the number of bcs
        C_num_bcs = 5;

        // specify boundary vertices references
        INLET    = 40;
        WALL     = 60;
        SLIPWALL = 61;
        OUTLET   = 50;
        CYLINDER = 70;
        // bind some locally defined methods to (boost) functions for use in BC, etc.
        u3D_zero = boost::bind(&CTcaseUser::u3DZero, this, _1, _2, _3, _4, _5);
        u3D_in = boost::bind(&CTcaseUser::u3Dcyl_dyn, this, _1, _2, _3, _4, _5);
    }

    // set boundary conditions
    void set_bcs()
    {
        std::cout << "  c- Setting user defined boundary conditions." << std::endl;

        zComp.resize(1);
        zComp[0] = 3;
        BCFunctionBase uIn(u3D_in);
        BCFunctionBase uZero(u3D_zero);

        // bc for the velocity
        C_bcHu->addBC("Inlet",    INLET, 	  Essential, Full,      uIn,   3);
        C_bcHu->addBC("Outlet",   OUTLET,	  Natural,   Full,      uZero, 3);
        C_bcHu->addBC("Wall",     WALL, 	  Essential, Full,      uZero, 3);
        C_bcHu->addBC("SlipWall", SLIPWALL, Essential, Component, uZero, zComp);
        C_bcHu->addBC("Cylinder", CYLINDER, Essential, Full,      uZero, 3);

        // bc for the pressure
        C_bcHp->addBC("Inlet",    INLET,    Natural,   Scalar,    uZero);
        C_bcHp->addBC("Outlet",   OUTLET,   Essential, Scalar,    uZero);
        C_bcHp->addBC("Wall",     WALL,     Natural,   Scalar,    uZero);
        C_bcHp->addBC("SlipWall", SLIPWALL, Natural,   Scalar,    uZero);
        C_bcHp->addBC("Cylinder", CYLINDER, Natural,   Scalar,    uZero);

    }

}; // struct CTcaseUser

#endif /* __CT_USER_CASE_HH */
