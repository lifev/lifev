
/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
 
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

#ifndef _ELEM_OPER_EXT
#define _ELEM_OPER_EXT

/******************************************************************
 Definition of ElementarOperators and all the utilities
 for the assembly via Expression Templates (See Veldhuizen, Haney, Furnish)
 
 Release 0.1 - A. Veneziani, January 2001
 Adaptation to the new f.e. JFG 04/2002 (remark : several class
 implementation should now be put in a .cc)
*****************************************************************/
#include "lifeV.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "currentFE.hpp"

namespace LifeV
{
typedef double Real;

class Function
{
public:
    virtual Real operator() ( Real x, Real y, Real z )
    {
        return 0;
    }
};


template <typename A>
class Tensor
{
public:
    Tensor( int n = NDIM )
    {
        m.resize( n );
        for ( int i = 0;i < n;i++ )
        {
            m[ i ].resize( n );
        }
    }

    //Operatore di assegnamento
    void operator() ( int i,int j, A a )
    {
        m[ i ][ j ] = a;
    }

    //Operatori di restituzione
    A operator() ( int i,int j )
    {
        return ( m[ i ][ j ] );
    }

    Real operator() ( int i,int j, Real x, Real y, Real z )
    {
        return ( *m[ i ][ j ] ) ( x, y, z );
    }

private:
    std::vector<std::vector<A> > m;
};

typedef Tensor<Function*> FTens;
typedef Tensor<Real> RTens;

//////////////////////
// Available elementar operators
/////////////////////
class Mass
{
public:
    Mass( CurrentFE* fe ) : _fe( fe )
    {}
    ;

    CurrentFE* fe_ptr()
    {
        return _fe;
    };

    Real operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic = 0, int jc = 0 )
    {
        return _fe->phi( i, iq ) * _fe->phi( j, iq );
    }

private:
    CurrentFE* _fe; // MUST be a pointer
};


// December 2001: Vector formulation
class VMass
{
public:
    VMass( CurrentFE* fe, int ncomp ) : _fe( fe ), _ncomp( ncomp )
    {}
    ;

    CurrentFE* fe_ptr()
    {
        return _fe;
    };

    // Versione semplice: DOVREBBE ANDARE BENE PER MOLTE COSE
    Real operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic, int jc )
    {
        if ( ic == jc )
            return _fe->phi( i, iq ) * _fe->phi( j, iq );
        else
            return 0;
    }


    // CHIEDRE AD ALAIN per quanto segue
    // For the VBR Pattern (Aztec): operator returns the single component of the local matrix
    // (for each quadrature node) which is a _ncomp*_ncomp vector organized as in the inner comment
    std::vector<Real> operator() ( int i, int j, int iq, Real x, Real y, Real z )
    {
        std::vector<Real> v( _ncomp * _ncomp );
        for ( int jc = 1;jc <= _ncomp;++jc )
        {
            for ( int ic = 1;ic <= _ncomp;++ic )
            {
                // [[a,b,c][d,e,f][g,h,i]] => [a,d,g,b,e,h,c,f,i];
                ic == jc ? v[ ( jc - 1 ) * _ncomp + ic ] = _fe->phi( i, iq ) * _fe->phi( j, iq ) : v[ ( jc - 1 ) * _ncomp + ic ] = 0.;
            }
        };
        return v;
    }

private:
    CurrentFE* _fe; // MUST be a pointer
    int _ncomp; // CHIEDERE A LUCA: IO LO SPECIFICHEREI STATIC CONST: che ne pensi ????
};

/////
class Stiff
{
public:
    Stiff( CurrentFE* fe ) : _fe( fe )
    {}
    ;
    // Stiff()
    //   { ASSERT_PRE(FE::getDerivSwitch(FirstDeriv),
    //      "stiffness matrix needs first derivatives");
    //   }

    CurrentFE* fe_ptr()
    {
        return _fe;
    };

    Real operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic = 0, int jc = 0 )
    {
        Real s = 0;
        for ( int icoor = 0; icoor < _fe->nbCoor;icoor++ )
        {
            s += _fe->phiDer( i, icoor, iq ) * _fe->phiDer( j, icoor, iq );
        }
        return s;
    }

    Real operator() ( int i, int j, int iq, int icoor, int jcoor, int ic = 0, int jc = 0 )
    {
        return _fe->phiDer( i, icoor, iq ) * _fe->phiDer( j, jcoor, iq );
    }


private:
    CurrentFE* _fe; // MUST be a pointer
};


// December 2001: Vector formulation
class VStiff
{
public:
    VStiff( CurrentFE* fe, int ncomp ) : _fe( fe ), _ncomp( ncomp )
    {}
    ;

    CurrentFE* fe_ptr()
    {
        return _fe;
    };

    Real by_comp( int i, int j, int iq, Real x, Real y, Real z )
    {
        Real s = 0;
        for ( int icoor = 0; icoor < _fe->nbCoor;icoor++ )
            s += _fe->phiDer( i, icoor, iq ) * _fe->phiDer( j, icoor, iq );

        return s;
    }

    // Versione semplice: DOVREBBE ANDARE BENE PER MOLTE COSE
    Real operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic, int jc )
    {
        if ( ic == jc )
        {
            Real s = 0;
            for ( int icoor = 0; icoor < _fe->nbCoor;icoor++ )
            {
                s += _fe->phiDer( i, icoor, iq ) * _fe->phiDer( j, icoor, iq );
            }
            return s;
        }
        else
            return 0;
    }


    // CHIEDERE AD ALAIN per quanto segue
    // For the VBR Pattern (Aztec): operator returns the single component of the local matrix
    // (for each quadtarure node) which is a _ncomp*_ncom vector organized as in the inner comment
    std::vector<Real> operator() ( int i, int j, int iq, Real x, Real y, Real z )
    {
        std::vector<Real> v( _ncomp * _ncomp );
        for ( int jc = 1;jc <= _ncomp;++jc )
        {
            for ( int ic = 1;ic <= _ncomp;++ic )
            {
                // [[a,b,c][d,e,f][g,h,i]] => [a,d,g,b,e,h,c,f,i];
                ic == jc ? v[ ( jc - 1 ) * _ncomp + ic ] = by_comp( i, j, iq, x, y, z ) : v[ ( jc - 1 ) * _ncomp + ic ] = 0.;
            }
        };
        return v;
    }

private:
    CurrentFE* _fe; // MUST be a pointer
    int _ncomp; // CHIEDERE A LUCA: IO LO SPECIFICHEREI STATIC CONST: che ne pensi ????
};


// A glance to the mixed (and NS) problems.....
template <int coor>
class Grad
{
public:
    Grad<coor>( CurrentFE* fe ) : _fe1( fe ), _fe2( fe )
    {}
    ;

    CurrentFE* fe1_ptr()
    {
        return _fe1;
    };
    CurrentFE* fe2_ptr()
    {
        return _fe2;
    };

    Grad<coor>( CurrentFE* fe1, CurrentFE* fe2 ) : _fe1( fe1 ), _fe2( fe2 ), _ncomp( 1 )
    {}
    ;
    // Grad()
    //   { ASSERT_PRE(FE::getDerivSwitch(FirstDeriv),
    //      "stiffness matrix needs first derivatives");
    //   }

    // Correction Alain, 25/01/02, ic and jc were missing.
    Real operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic = 0, int jc = 0 )
    {
        return _fe2->phi( i, iq ) * _fe1->phiDer( j, coor, iq );
    }

private:
    CurrentFE* _fe1; // MUST be pointers
    CurrentFE* _fe2; // MUST be pointers
    int _ncomp;
};

// A vectorial mixed operator:
// psi(j)*div(v_phi(i)) where psi is FE1 and vectorial v_phi is FE2
class Vdiv
{
public:
    Vdiv( CurrentFE * fe, int ncomp ) : _fe1( fe ), _fe2( fe ), _ncomp( ncomp )
    {}

    CurrentFE* fe1_ptr()
    {
        return _fe1;
    }
    CurrentFE* fe2_ptr()
    {
        return _fe2;
    }

    Vdiv( CurrentFE * fe1, CurrentFE * fe2, int ncomp ) : _fe1( fe1 ), _fe2( fe2 ), _ncomp( ncomp )
    {}

    Real operator() ( int i, int j, int iq, Real x, Real y,
                      Real z, int ic, int jc )
    {
        return -_fe1->phi( i, iq ) * _fe2->phiDer( j, jc, iq );
    }

private:
    CurrentFE* _fe1; // MUST be pointers
    CurrentFE* _fe2; // MUST be pointers
    int _ncomp;
};

// trVdiv Transpose operator of Vdiv:
// psi(i)*div(v_phi(j)) where psi is FE1 and vectorial v_phi is FE2
class trVdiv
{
public:

    trVdiv( CurrentFE * fe, int ncomp ) : _fe1( fe ), _fe2( fe ), _ncomp( ncomp )
    {}

    CurrentFE* fe1_ptr()
    {
        return _fe1;
    }
    CurrentFE* fe2_ptr()
    {
        return _fe2;
    }

    trVdiv( CurrentFE * fe1, CurrentFE * fe2, int ncomp ) : _fe1( fe1 ), _fe2( fe2 ), _ncomp( ncomp )
    {}

    Real operator() ( int i, int j, int iq, Real x, Real y,
                      Real z, int ic, int jc )
    {
        return -_fe2->phi( j, iq ) * _fe1->phiDer( i, ic, iq );
    }

private:
    CurrentFE* _fe1; // MUST be pointers
    CurrentFE* _fe2; // MUST be pointers
    int _ncomp;
};


/////////////////////////////////////////////
// WRAPPER Classes
//
// P = Real for a scalar problem
// P = ElemMat for a vector problem
//
template <typename P, typename A>
class EOExpr
{
private:
    A _a;
public:
    EOExpr( const A& eo ) : _a( eo )
    {}

    P operator() ( int i, int iq, Real x, Real y, Real z, int ic )
    {
        return _a( i, iq, x, y, z, ic );
    }

    P operator() ( int i, int j, int iq, int icoor, int jcoor )
    {
        return _a( i, j, iq, icoor, jcoor );
    }

    P operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic, int jc )
    {
        return _a( i, j, iq, x, y, z, ic, jc );
    }

    P operator() ( int i, int j, int iq, int icoor, int jcoor, int ic, int jc )
    {
        return _a( i, j, iq, icoor, jcoor, ic, jc );
    }

    //   P operator()(int i, int iq, Real x, Real y, Real z,int ic=0,int jc=0)
    //     { return _a(i,iq,x,y,z);}



    //  int nc() {return _a.nc();}

};
//////////////////////////////////////////////
//
/*template<typename P, typename A>
class SymEOExpr: public EOExpr<P,A>{
 public:
   SymEOExpr(const A& eo): EOExpr<P,A>(eo)
     {}
     };*/ 
//////////////////////////////////////////////
//
// Example of use: -mu*Laplace(u) + sigma*(u)
//
// Mass<fe> mass();
// Stiff<fe> stiff();
// // Grad<fe,0> gradx();
// // Grad<fe,1> grady();
// // Grad<fe,2> gradz();
// Real mu=1.;
// Real sigma = 1.;
// Real beta1=0.1, beta2=0.2 , beta3=0.3;
// // typedef EOExpr<Real,EORBinOp<Real,Mass,EORMult<Real>>> OMass;
// // typedef EOExpr<Real,EORBinOp<Real,Stiff,EORMult<Real>>> OStiff;
// // typedef EOExpr<Real,EOBinOp<Real,OMass,OStiff,EOAdd<Real>>> Oper;
// // Oper oper=mu*stiff+sigma*mass
// assemble(oper,mesh,geomap,fe,dof,rhs,A,b);
//
/////////////////////////////////
//
// Binary operations
//
template <typename P, typename A, typename B, typename Op>
class EOBinOp
{
private:
    A _a;
    B _b;
    // Prova Gennaio 2002
    //  int _ncomp;
public:
    EOBinOp( const A& a, const B& b ) : _a( a ), _b( b )
    {
        //     if (a.nc()!=b.nc()) {std::cerr << "Error: you can't compose differential operators with different components"<<std::endl; exit(1);}
        //     _ncomp=a.nc();
    };

    P operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic = 0, int jc = 0 )
    {
        return Op::apply( _a( i, j, iq, x, y, z, ic, jc ), _b( i, j, iq, x, y, z, ic, jc ) );
    };

    P operator() ( int i, int iq, Real x, Real y, Real z, int ic = 0 )
    {
        return Op::apply( _a( i, iq, x, y, z, ic ), _b( i, iq, x, y, z, ic ) );
    };

};
//
template <typename P, typename A, typename Op>
class EORBinOp
{
private:
    A _a;
    Real _b;
    //  int _ncomp;
public:
    EORBinOp( const A& a, const Real b ) : _a( a ), _b( b )
    {}
    ;
    EORBinOp( const Real b, const A& a ) : _a( a ), _b( b )
    {}
    ;

    P operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic = 0, int jc = 0 )
    {
        return Op::apply( _a( i, j, iq, x, y, z, ic, jc ), _b );
    };
};

template <typename P, typename A, typename Op>
class EOFTBinOp
{
private:
    A _a;
    FTens _t;
    int _ncomp;
public:
    EOFTBinOp( const FTens& t, const A& a ) : _a( a ), _t( t )
    {}
    ;
    EOFTBinOp( const A& a, const FTens& t ) : _a( a ), _t( t )
    {}
    ;

    P operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic = 0, int jc = 0 ) const
    {
        P s = 0;
        for ( int icoor = 0; icoor < NDIM; icoor++ )
        {
            for ( int jcoor = 0; jcoor < NDIM; jcoor++ )
            {
                s += Op::apply( _a( i, j, iq, icoor, jcoor, ic, jc ), _t( icoor, jcoor, x, y, z ) );
            }
        }
        return s;
    }

};

template <typename P, typename A, typename Op>
class EORTBinOp
{
private:
    A _a;
    RTens _t;
    int _ncomp;
public:
    EORTBinOp( const RTens& t, const A& a ) : _a( a ), _t( t )
    {}
    ;
    EORTBinOp( const A& a, const RTens& t ) : _a( a ), _t( t )
    {}
    ;

    P operator() ( int i, int j, int iq, Real x, Real y, Real z, int ic = 0, int jc = 0 ) const
    {
        P s = 0;
        for ( int icoor = 0; icoor < NDIM; icoor++ )
        {
            for ( int jcoor = 0; jcoor < NDIM; jcoor++ )
            {
                s += Op::apply( _a( i, j, iq, icoor, jcoor, ic, jc ), _t( icoor, jcoor ) );
            }
        }
        return s;
    }

};

template <typename P, typename A, typename Op>
class EOFBinOp
{
private:
    A _a;
    Function* _f;
    int _ncomp;
public:
    EOFBinOp( const A& a, const Function& f ) : _a( a ), _f( &f )
    {}
    ;
    EOFBinOp( const Function& f, const A& a ) : _a( a ), _f( &f )
    {}
    ;


    P operator() ( int i, int j, int iq, Real x, Real y, Real z ) const
    {
        return Op::apply( _a( i, j, iq, x, y, z ), ( *_f ) ( x, y, z ) );
    };

};



//
// Low-level operations
//
////////////////////////////////
//
// Addition
//
template <typename P>
class EOAdd
{
public:
    EOAdd()
    {}
    ;

    static inline P apply( P a, P b )
    {
        return a + b;
    }
};
//////////////////////
//
// Scalar multiplication
template <typename P>
class EORMult
{
public:
    EORMult()
    {}
    ;
    static inline P apply( P a, Real b )
    {
        return a * b;
    }
};
///////////////////////////
//
// Binary Operators
//
template <typename P, typename A, typename B>
EOExpr<P, EOBinOp<P, EOExpr<P, A>, EOExpr<P, B>, EOAdd<P> > >
operator+( const EOExpr<P, A>& a, const EOExpr<P, B>& b )
{
    typedef EOBinOp<P, EOExpr<P, A>, EOExpr<P, B>, EOAdd<P> > ExprT;
    return EOExpr<P, ExprT>( ExprT( a, b ) );
}

template <typename P, typename A>
EOExpr<P, EORBinOp<P, EOExpr<P, A>, EORMult<P> > >
operator*( const EOExpr<P, A>& a, const Real b )
{
    typedef EORBinOp<P, EOExpr<P, A>, EORMult<P> > ExprT;
    return EOExpr<P, ExprT>( ExprT( a, b ) );
}

template <typename P, typename A>
EOExpr<P, EORBinOp<P, EOExpr<P, A>, EORMult<P> > >
operator*( const Real b, const EOExpr<P, A>& a )
{
    typedef EORBinOp<P, EOExpr<P, A>, EORMult<P> > ExprT;
    return EOExpr<P, ExprT>( ExprT( a, b ) );
}

template <typename P, typename A>
EOExpr<P, EOFTBinOp<P, EOExpr<P, A>, EORMult<P> > >
operator*( const FTens& t, const EOExpr<P, A>& a )
{
    typedef EOFTBinOp<P, EOExpr<P, A>, EORMult<P> > ExprT;
    return EOExpr<P, ExprT>( ExprT( t, a ) );
}

template <typename P, typename A>
EOExpr<P, EORTBinOp<P, EOExpr<P, A>, EORMult<P> > >
operator*( const RTens& t, const EOExpr<P, A>& a )
{
    typedef EORTBinOp<P, EOExpr<P, A>, EORMult<P> > ExprT;
    return EOExpr<P, ExprT>( ExprT( t, a ) );
}

template <typename P, typename A>
EOExpr<P, EOFBinOp<P, EOExpr<P, A>, EORMult<P> > >
operator*( const Function& f, const EOExpr<P, A>& a )
{
    typedef EOFBinOp<P, EOExpr<P, A>, EORMult<P> > ExprT;
    return EOExpr<P, ExprT>( ExprT( a, f ) );
}

template <typename P, typename A>
EOExpr<P, EOFBinOp<P, EOExpr<P, A>, EORMult<P> > >
operator*( const EOExpr<P, A>& a, const Function& f )
{
    typedef EOFBinOp<P, EOExpr<P, A>, EORMult<P> > ExprT;
    return EOExpr<P, ExprT>( ExprT( a, f ) );
}

// STABILIZATION STUFF

/****************AJOUT DE PATRICK************************/
// Test d'AJOUT de PHIFCT sur l'operateur ASSEMBLE
//===========================================
// Function TSquare: da 3 a 5 volte + veloce di pow(InValue,2)
//===========================================
template < typename T >
inline T TSquare( T InValue )
{
    return ( InValue * InValue );
}

//===========================================
// Class PhiFct: Serve per l`assemblaggio del vettore dei termini noti del sistema lineare
//===========================================
template <typename FCT>
class PhiFct
{
public:

    PhiFct<FCT>( CurrentFE* fe, FCT& fcT ) : _fe( fe ), _fcT( fcT )
    {}
    CurrentFE* fe_ptr()
    {
        return _fe;
    }

    Real operator() ( int j,int iq, Real x, Real y, Real z )
    {
        return _fe->phi( j, iq ) * _fcT( x, y, z );
    }

    Real operator() ( int j,int iq, Real x, Real y, Real z, int ic )
    {
        return _fe->phi( j, iq ) * _fcT( x, y, z, ic );
    }


    Real operator() ( int j,int iq, Real x, Real y, Real z, Real t, int ic )
    {
        return _fe->phi( j, iq ) * _fcT( x, y, z, t, ic );
    }

    Real operator() ( int j,int iq, Real x, Real y, Real z, Real t )
    {
        return _fe->phi( j, iq ) * _fcT( x, y, z, t );
    }

private:

    CurrentFE* _fe;
    FCT _fcT;
};

//===========================================
// Class PurePhi
//===========================================
class PurePhi
{
public:

    PurePhi( CurrentFE* fe ) : _fe( fe )
    {}

    CurrentFE* fe_ptr()
    {
        return _fe;
    }

    Real operator() ( UInt k,UInt iq )
    {
        return _fe->phi( k, iq );
    }

private:

    CurrentFE* _fe;
};

//===========================================
// Class PureGrad
//===========================================
class PureGrad
{
public:

    PureGrad( CurrentFE* fe ) : _fe( fe )
    {}

    CurrentFE* fe_ptr()
    {
        return _fe;
    }

    Real operator() ( int k,int icoor, int iq )
    {
        return _fe->phiDer( k, icoor, iq );
    }

private:

    CurrentFE* _fe;
};

//===========================================
// Class PureDer2
//===========================================
class PureDer2
{
public:

    PureDer2( CurrentFE* fe ) : _fe( fe )
    {}

    CurrentFE* fe_ptr()
    {
        return _fe;
    }

    Real operator() ( UInt k ,UInt icoor, UInt jcoor, UInt iq )
    {
        return _fe->phiDer2( ( int ) k, ( int ) icoor, ( int ) jcoor, ( int ) iq );
    }

private:

    CurrentFE* _fe;
};

// Operatori da usare nell`implementazione dei metodi di stabilizzazione:
//===========================================
// Class Laplacian
//===========================================
template < typename PHD >
class Laplacian
{
public:

    Laplacian< PHD > ( PHD& phD ) : _phD( phD )
    {}

    Real operator() ( UInt k ,UInt iq )
    {
        Real s = 0;

        for ( int icoor = 0;icoor < _phD.fe_ptr() ->nbCoor;icoor++ )
        {
            s += _phD( k, icoor, icoor, iq );
        }

        return s;
    }

private:

    PHD _phD;
};

//===========================================
// Class LS
//===========================================
template < typename L>
class LS
{
public:

    LS<L>( L &Lap, CurrentFE* fe, Real sigma, Real mu ) : _Lap( Lap ), _fe( fe ), _sigma( sigma ), _mu( mu )
    {}
    CurrentFE* fe_ptr()
    {
        return _fe;
    }
    Real operator() ( UInt i,UInt iq )
    {
        return _mu * _Lap( i, iq ) + _sigma * _fe->phi( i, iq );

    }

private:

    L _Lap;
    CurrentFE* _fe;
    Real _sigma;
    Real _mu;

};

//===========================================
// Class LSS
//===========================================
template <typename PHD, typename VT>
class LSS
{
public:

    LSS<PHD, VT>( PHD& phD, VT& Beta ) : _phD( phD ), _Beta( Beta )
    {}

    Real operator() ( UInt k,UInt iq )
    {
        Real s = 0;

        for ( int icoor = 0;icoor < _phD.fe_ptr() ->nbCoor;icoor++ )
        {
            s += _Beta[ icoor ] * _phD( k, icoor, iq );
        }

        return s;
    }

private:

    PHD _phD;
    VT _Beta;
};

//===========================================
// Class Stab
//===========================================
// The Stab class has to be used for problems which have a regular mesh.
// In fact, the parameter delta which serves for the stabilization is fixed at the beginning of the program. You have to calculate it before (depending on the mesh).
template <typename LS, typename LSS, typename FCT>
class Stab
{
public:

    Stab<LS, LSS, FCT>( LS& ols, LSS& olss, FCT& fct, Real Delta, Real Ro ) : _ols( ols ), _olss( olss ), _fct( fct ), _Delta( Delta ), _Ro( Ro )
    {}

    Real operator() ( UInt i,UInt j, UInt iq, Real x , Real y , Real z, UInt ic = 0, UInt jc = 0 )
    {

        return _Delta * ( ( _ols( i, iq ) + _olss( i, iq ) ) * ( _olss( j, iq ) + _Ro * ( _ols( j, iq ) ) ) );
    }

    Real operator() ( UInt j,UInt iq, Real x, Real y, Real z, UInt ic = 0, UInt jc = 0 )
    {
        return _Delta * ( ( _fct( x, y, z ) ) * ( _olss( j, iq ) + _Ro * ( _ols( j, iq ) ) ) );
    }

private:

    LS _ols;
    LSS _olss;
    FCT _fct;
    Real _Delta;
    Real _Ro;
};


// Class VStab
//==============================================
template <typename LS, typename LSS, typename FCT>
class VStab
{
public:
    VStab<LS, LSS, FCT>( LS& ols, LSS& olss, FCT& fct, Real Delta, Real Ro ) : _ols( ols ), _olss( olss ), _fct( fct ), _Delta( Delta ), _Ro( Ro )
    {}

    Real operator() ( UInt i,UInt j, UInt iq, Real x, Real y, Real z, UInt ic, UInt jc )
    {
        if ( ic == jc )
            return _Delta * ( ( _ols( i, iq ) + _olss( i, iq ) ) * ( _olss( j, iq ) + _Ro * ( _ols( j, iq ) ) ) );
        else
            return 0.;
    }

    Real operator() ( UInt j,UInt iq, Real x, Real y, Real z, UInt ic = 0, UInt jc = 0 )
    {
        return _Delta * ( ( _fct( x, y, z ) ) * ( _olss( j, iq ) + _Ro * ( _ols( j, iq ) ) ) );
    }

private:

    LS _ols;
    LSS _olss;
    FCT _fct;
    Real _Delta;
    Real _Ro;
};



//=====================================================================================================
// Operator HMStab
//=====================================================================================================
// The HmStab class  has the same goal as the Stab class but it allows the user to stabilize the problem even if the mesh is non-regular.
// In fact, in this case the parameter of stabilization delta is calculated for each tetraedre.
template <typename OLS, typename OLSS, typename FCT, typename MESH>
class HMStab
{
public:

    HMStab<OLS, OLSS, FCT, MESH>( OLS& ols, OLSS& olss, FCT& fct, MESH& Me, CurrentFE* fe, Real& Mu, int& Ro, std::vector<double>& Tv );

    void updateLocalElement();

    Real operator() ( UInt i,UInt j, UInt ig, Real x, Real y, Real z, UInt ic = 0, UInt jc = 0 )
    {
        if ( ic == jc )  // da verificare !!!!!!!!!! A. Veneziani - Novembre 2002
        {
            if ( i == 0 && j == 0 && ig == 0 )
            {
                updateLocalElement();

            }
            return ( _Delta * ( ( _ols( i, ig ) + _olss( i, ig ) ) * ( _olss( j, ig ) + _Ro * ( _ols( j, ig ) ) ) ) );
        }
        else
            return 0.;
    }

    Real operator() ( int j,int ig, Real x, Real y, Real z, int ic )
    {
        return _Delta * _fct( x, y, z, ic ) * ( _olss( j, ig ) + _Ro * ( _ols( j, ig ) ) );
    }

private:

    OLS _ols;    // Parte simmetrica dell'operatore
    OLSS _olss;           // Parte antisimmetrica dell'operatore
    FCT _fct;
    // Termine noto dell'equazione differenziale
    Real _Delta;           // Parametro di stabilizzazione da calcolare
    int _Ro;    // Seleziono il metodo di stabilizzazione (costruttore)
    Real _V;    // Volume Tetraedro (calcolare con geomap)
    Real _T;    // Modulo del trasporto (calcolare dal vettore di trasporto: costruttore)
    Real _S;    // Sfericita' (calcolare per ogni tetraedro)
    Real _At;    // Area totale del tetraedro in esame
    Real _Pe;    // Peclet (viene calcolato per ogni tetredro)
    Real _Mu;    // Diffusione (inizializzata dal costruttore)
    int _Cn;    // Contatore tetredri (inizializzata dal costruttore a 1)

    std::vector<double> _Tv;         // Vettore componenti trasporto (costruttore)

    int _C[ 4 ][ 3 ];   // Matrice connessione del tetraedro (inizializzata dal costruttore)
    Real _Vcoor[ 4 ][ 3 ];  // Coordinate dei vertici del tetraedro
    Real _L[ 6 ];   // Lunghezza dei lati
    Real _P[ 4 ];   // Semiperimetri
    Real _A[ 4 ];   // Aree delle facce

    MESH const & _Me;  // Riferimento alla RegionMesh3D (costruttore): estraggo le coordinate dei vertici
    CurrentFE* _fe;   // Riferimento al CurrentFE (costruttore): calcolo i volumi con measure
};

//---------------------------------------------------------------------------------------------------------------------------
// inline HMStab<OLS,OLSS,FCT>::HMStab(OLS& ols,OLSS& olss,FCT& fct,int Ro)
//---------------------------------------------------------------------------------------------------------------------------
template <typename OLS, typename OLSS, typename FCT, typename MESH>
inline HMStab<OLS, OLSS, FCT, MESH>::HMStab( OLS& ols, OLSS& olss, FCT& fct, MESH& Me, CurrentFE* fe, Real& Mu, int& Ro, std::vector<double>& Tv ) : _ols( ols ), _olss( olss ), _fct( fct ), _Me( Me ), _fe( fe )  //,_Tv(Tv), _Mu(Mu),_Ro(Ro),_V(0.0) ,_T(0.),_S(0.0) ,_Delta(Delta),_At(0.0) ,_Pe(0.0) ,_Cn(1)
{
    _Cn = 1;
    _Pe = 0.0;
    _At = 0.0;
    _Delta = 0.0;
    _S = 0.0;
    _T = 0.0;
    _V = 0.0;
    _Ro = Ro;
    _Mu = Mu;
    _Tv = Tv;
    int i = 0;
    int j = 0;
    Real somma = 0;

    // Inizializzo la matrice di connessione
    _C[ 0 ][ 0 ] = 0 ;
    _C[ 1 ][ 0 ] = 3 ;
    _C[ 2 ][ 0 ] = 0 ;
    _C[ 3 ][ 0 ] = 1 ;
    _C[ 0 ][ 1 ] = 2 ;
    _C[ 1 ][ 1 ] = 4 ;
    _C[ 2 ][ 1 ] = 1 ;
    _C[ 3 ][ 1 ] = 2 ;
    _C[ 0 ][ 2 ] = 4 ;
    _C[ 1 ][ 2 ] = 5 ;
    _C[ 2 ][ 2 ] = 3 ;
    _C[ 3 ][ 2 ] = 5 ;

    // Inizializzo le coordinate dei vertici
    for ( i = 0;i < 4;i++ )
    {
        for ( j = 0;j < 3;j++ )
        {
            _Vcoor[ i ][ j ] = 0.;
        }
    }

    // Inizializzo le lunghezze dei lati
    for ( i = 0;i < 6;i++ )
    {
        _L[ i ] = 0.;
    }

    // Inizializzo i semiperimetri e le Aree
    for ( i = 0;i < 4;i++ )
    {
        _P[ i ] = 0.;
        _A[ i ] = 0.;
    }

    // Inizializzo il modulo del trasporto
    for ( i = 0;i < 3;i++ )
    {
        somma += TSquare( _Tv[ i ] );               //pow(_Tv[i],2);
    }
    _T = pow( somma, .5 );
}

//---------------------------------------------------------------------------------------------------------------------------
// inline void HMStab<OLS,OLSS,FCT>::updateLocalElement()
//---------------------------------------------------------------------------------------------------------------------------
template <typename OLS, typename OLSS, typename FCT, typename MESH>
inline void HMStab<OLS, OLSS, FCT, MESH>::updateLocalElement()
{
    int i = 0;
    int j = 0;
    int k = 0;
    Real somma = 0;

    _fe->updateJac( _Me.volumeList( _Cn ) );
    _V = _fe->measure();
    //std::cout << "Volume " << _V << std::endl;

    for ( j = 1;j < 5;j++ )
    {
        _Vcoor[ j - 1 ][ 0 ] = _Me.volumeList( _Cn ).point( j ).x();
        _Vcoor[ j - 1 ][ 1 ] = _Me.volumeList( _Cn ).point( j ).y();
        _Vcoor[ j - 1 ][ 2 ] = _Me.volumeList( _Cn ).point( j ).z();
    }

    // Calcolo i lati del tetredro in esame:

    for ( j = 0;j < 3;j++ )
    {
        somma = 0;
        for ( i = 0;i < 3;i++ )
            somma += TSquare( _Vcoor[ 0 ][ i ] - _Vcoor[ j + 1 ][ i ] );    // pow(_Vcoor[0][i]-_Vcoor[ j + 1 ][i],2) ;
        _L[ j ] = sqrt( somma );                               // pow(somma,.5) ;
    }

    for ( j = 1;j < 3;j++ )
    {
        somma = 0;
        for ( i = 0;i < 3;i++ )
            somma += TSquare( _Vcoor[ 1 ][ i ] - _Vcoor[ j + 1 ][ i ] );    // pow(_Vcoor[1][i]-_Vcoor[ j + 1 ][i],2) ;
        _L[ j + 2 ] = sqrt( somma );                            // pow(somma,.5);
    }

    for ( j = 1;j < 2;j++ )
    {
        somma = 0;
        for ( i = 0;i < 3;i++ )
            somma += TSquare( _Vcoor[ 2 ][ i ] - _Vcoor[ j + 2 ][ i ] );    // pow(_Vcoor[2][i]-_Vcoor[ j + 2 ][i],2) ;
        _L[ 5 ] = sqrt( somma ) ;                                // pow(somma,.5);
    }

    // Calcolo i semiperimetri del tetredro in esame:

    for ( i = 0;i < 4;i++ )
    {
        somma = 0;
        for ( j = 0;j < 3;j++ )
        {
            k = _C[ i ][ j ];
            somma += _L[ k ];
        }
        _P[ i ] = somma / 2;
    }

    // Calcolo le aree delle facce del tetredro in esame:

    for ( i = 0;i < 4;i++ )
    {
        somma = 1;
        for ( j = 0;j < 3;j++ )
        {
            k = _C[ i ][ j ];
            somma *= _P[ i ] - _L[ k ];
        }
        _A[ i ] = sqrt( _P[ i ] * somma );                          // pow(_P[i] * somma,.5)  ;
    }

    // Calcolo l'area totale:

    _At = 0;
    for ( i = 0;i < 4;i++ )
        _At += _A[ i ];

    // Calcolo la sfericita':
    _S = 6 * _V / _At;
    //_S = 0.042265;

    // Calcolo il numero di Peclet:
    _Pe = ( _T * _S ) / ( 2 * _Mu );
    // std::cout <<  " valeur du numero de Peclet :  " <<  _Pe << std::endl;
    //std::cout << std::endl << "Peclet " << _Pe ;

    // Assegno il valore del parametro di stabilizzazione:
    if ( _Pe >= 1 )
        _Delta = 0.5 * ( _S * TSquare( _T ) ) / ( sqrt( TSquare( _T ) * _V ) );
    else
        _Delta = 0;
    _Cn++;
    // std::cout <<  " valeur de delta :  " <<  _Delta << std::endl;
}



//===========================================================================================================
// Operator StabUP
//===========================================================================================================
template <typename STIFF, typename MESH>
class StabUP
{
public:

    StabUP<STIFF, MESH>( STIFF& stiff, MESH& Me, CurrentFE* fe, std::vector<double>& Tv );

    void updateLocalElementUP();

    Real operator() ( UInt i ,UInt j, UInt ig, Real x, Real y, Real z, UInt ic = 0, UInt jc = 0 )
    {
        if ( i == 0 && j == 0 && ig == 0 )
        {
            updateLocalElementUP();
        }
        //return (0.0);
        // std::cout << _S << " " << _T << std::endl;
        return ( _S * _T * _stiff( i, j, ig, x, y, z ) );
        //return ( 0.12*0.042265*_stiff(i,j,ig,x,y,z));
    }


private:

    STIFF _stiff;    // Parte simmetrica dell'operatore

    Real _V;    // Volume Tetraedro (calcolare con geomap)
    Real _T;    // Modulo del trasporto (calcolare dal vettore di trasporto: costruttore)
    Real _S;    // Sfericita' (calcolare per ogni tetraedro)
    Real _At;    // Area totale del tetraedro in esame

    int _Cn;    // Contatore tetredri (inizializzata dal costruttore a 1)

    std::vector<double> _Tv;         // Vettore componenti trasporto (costruttore)

    int _C[ 4 ][ 3 ];   // Matrice connessione del tetraedro (inizializzata dal costruttore)
    Real _Vcoor[ 4 ][ 3 ];  // Coordinate dei vertici del tetraedro
    Real _L[ 6 ];   // Lunghezza dei lati
    Real _P[ 4 ];   // Semiperimetri
    Real _A[ 4 ];   // Aree delle facce

    MESH const & _Me;  // Riferimento alla RegionMesh3D (costruttore): estraggo le coordinate dei vertici
    CurrentFE* _fe;  // Riferimento al CurrentFE: calcolo i volumi con measure

};

//-----------------------------------------------------------------------------------------------------------------
// inline StabUP<OLS,OLSS,FCT>::StabUP(OLS& ols,OLSS& olss,FCT& fct,int Ro)
//-----------------------------------------------------------------------------------------------------------------
template <typename STIFF, typename MESH>
inline StabUP<STIFF, MESH>::StabUP ( STIFF& stiff, MESH& Me, CurrentFE* fe , std::vector<double>& Tv ) : _stiff( stiff ), _Me( Me ), _fe( fe )
{
    _Cn = 1;
    _At = 0.0;
    _S = 0.0;
    _T = 0.0;
    _V = 0.0;
    _Tv = Tv;
    int i = 0;
    int j = 0;
    Real somma = 0;

    // Inizializzo la matrice di connessione
    _C[ 0 ][ 0 ] = 0;
    _C[ 1 ][ 0 ] = 3;
    _C[ 2 ][ 0 ] = 0;
    _C[ 3 ][ 0 ] = 1;
    _C[ 0 ][ 1 ] = 2;
    _C[ 1 ][ 1 ] = 4;
    _C[ 2 ][ 1 ] = 1;
    _C[ 3 ][ 1 ] = 2;
    _C[ 0 ][ 2 ] = 4;
    _C[ 1 ][ 2 ] = 5;
    _C[ 2 ][ 2 ] = 3;
    _C[ 3 ][ 2 ] = 5;

    // Inizializzo le coordinate dei vertici
    for ( i = 0;i < 4;i++ )
    {
        for ( j = 0;j < 3;j++ )
        {
            _Vcoor[ i ][ j ] = 0.;
        }
    }
    // Inizializzo le lunghezze dei lati
    for ( i = 0;i < 6;i++ )
    {
        _L[ i ] = 0.;
    }
    // Inizializzo i semiperimetri e le Aree
    for ( i = 0;i < 4;i++ )
    {
        _P[ i ] = 0.;
        _A[ i ] = 0.;
    }
    // Inizializzo il modulo del trasporto
    for ( i = 0;i < 3;i++ )
    {
        somma += TSquare( _Tv[ i ] );               //pow(_Tv[i],2);
    }
    _T = pow( somma, .5 );

}

//---------------------------------------------------------------------------------------------------------------------------
// inline void StabUP<OLS,OLSS,FCT>::updateLocalElementUP()
//---------------------------------------------------------------------------------------------------------------------------
template <typename STIFF, typename MESH>
inline void StabUP<STIFF, MESH>::updateLocalElementUP()
{
    int i = 0;
    int j = 0;
    int k = 0;
    Real somma = 0;

    _fe->updateJac( _Me.volumeList( _Cn ) );
    _V = _fe->measure();

    for ( j = 1;j <= 4;j++ )
    {
        _Vcoor[ j - 1 ][ 0 ] = _Me.volumeList( _Cn ).point( j ).x();
        _Vcoor[ j - 1 ][ 1 ] = _Me.volumeList( _Cn ).point( j ).y();
        _Vcoor[ j - 1 ][ 2 ] = _Me.volumeList( _Cn ).point( j ).z();
    }

    // Calcolo i lati del tetredro in esame:

    for ( j = 0;j < 3;j++ )
    {
        somma = 0;
        for ( i = 0;i < 3;i++ )
            somma += TSquare( _Vcoor[ 0 ][ i ] - _Vcoor[ j + 1 ][ i ] );    // pow(_Vcoor[0][i]-_Vcoor[ j + 1 ][i],2) ;
        _L[ j ] = sqrt( somma );                               // pow(somma,.5) ;
    }

    for ( j = 1;j < 3;j++ )
    {
        somma = 0;
        for ( i = 0;i < 3;i++ )
            somma += TSquare( _Vcoor[ 1 ][ i ] - _Vcoor[ j + 1 ][ i ] );    // pow(_Vcoor[1][i]-_Vcoor[ j + 1 ][i],2) ;
        _L[ j + 2 ] = sqrt( somma );                            // pow(somma,.5);
    }

    for ( j = 1;j < 2;j++ )
    {
        somma = 0;
        for ( i = 0;i < 3;i++ )
            somma += TSquare( _Vcoor[ 2 ][ i ] - _Vcoor[ j + 2 ][ i ] );    // pow(_Vcoor[2][i]-_Vcoor[ j + 2 ][i],2) ;
        _L[ 5 ] = sqrt( somma );                                // pow(somma,.5);
    }

    // Calcolo i semiperimetri del tetredro in esame:
    for ( i = 0;i < 4;i++ )
    {
        somma = 0;
        for ( j = 0;j < 3;j++ )
        {
            k = _C[ i ][ j ];
            somma += _L[ k ];
        }
        _P[ i ] = somma / 2;
    }

    // Calcolo le aree delle facce del tetredro in esame:
    for ( i = 0 ; i < 4 ; i ++ )
    {
        somma = 1;
        for ( j = 0;j < 3;j++ )
        {
            k = _C[ i ][ j ];
            somma *= _P[ i ] - _L[ k ];
        }
        _A[ i ] = sqrt( _P[ i ] * somma );                          // pow(_P[i] * somma,.5)  ;
    }

    // Calcolo l'area totale:
    _At = 0;
    for ( i = 0;i < 4;i++ )
        _At += _A[ i ];

    // Calcolo la sfericita':
    _S = 6 * _V / _At;
    _Cn++;
}
}
#endif

