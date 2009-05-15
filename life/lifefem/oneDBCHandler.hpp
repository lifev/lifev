/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file oneDBC.hpp
  \author Lucia Mirabella <lucia@mathcs.emory.edu>
  \author Tiziano Passerini <tiziano@mathcs.emory.edu>
  \date 08/2006

  \brief This file contains a class for the boundary conditions for 1D models.

*/

#ifndef _ONED_BCHANDLER_HPP_
#define _ONED_BCHANDLER_HPP_

// #include <life/lifesolver/oneDNonLinModelParam.hpp>
#include <life/lifefem/oneDBCFunctions.hpp>

namespace ublas = boost::numeric::ublas;

namespace LifeV
{

/*!
 *\class OneDBC
 * Class featuring methods to handle boundary conditions.
 */
template< class FLUX >
class OneDBC
{
public:
    /** @name Constructor
     */
    //@{
    //! Constructor
    OneDBC( const OneDModelHandler::ScalVec_vector& U_thistime,
            const FLUX& fluxFun,
            const Real& dimDof,
            const std::string& side );
    //@}
    /** @name Destructor
     */
    //@{
    ~OneDBC(){}
    //@}
    /** @name Public Methods
     */
    //@{
    //! Compute [A,Q] at the boundary
    OneDModelHandler::Vec2D Uboundary(const OneDModelHandler::OneDModelHandler::ScalVec& U1,
                                      const OneDModelHandler::OneDModelHandler::ScalVec& U2) const;
    //! Return the boundary Dof
    inline UInt boundaryDof() const{return _M_boundaryDof;}
    //! Apply boundary conditions
    void applyBC(const Real& time_val, OneDModelHandler::Vec2D& BC_dir);

    // careful I need them to be setters!
    // do not remove the reference
    inline OneDBCFunctionPointer& rhs( const std::string& line )
    { return _M_rhs_at_line[line]; }

    inline std::string& variable( const std::string& line )
    { return _M_variable_at_line[line]; }

    inline OneDModelHandler::Vec2D& matrixrow( const std::string& line )
    { return _M_matrixrow_at_line[line]; }

    //    OneDBCFunctionPointer& f2(){ return _M_f2; }

    inline bool& isInternal() { return _M_isInternal; }

    //@}
protected:
    /** @name Protected Methods
     */
    //@{
    //! Impose the chosen boundary condition
    void compute_resBC( const Real& time_val );

    void compute_resBC_line( std::string line, OneDModelHandler::Vec2D left_eigvec,
                             OneDModelHandler::Vec2D U, OneDModelHandler::Vec2D W, Real& rhs );

    //! solve a 2x2 linear system by the Cramer method (for the boundary systems)
    OneDModelHandler::Vec2D _solveLinearSyst2x2(const OneDModelHandler::Vec2D& line1,
                                                const OneDModelHandler::Vec2D& line2,
                                                const OneDModelHandler::Vec2D& rhs2d) const;
    //@}

    //    OneDBCFunctionPointer _M_f1, _M_f2;

    bool _M_isInternal;

    std::map<std::string, std::string> _M_variable_at_line;

    std::map<std::string, OneDModelHandler::Vec2D> _M_matrixrow_at_line;

    std::map<std::string, OneDBCFunctionPointer> _M_rhs_at_line;

    std::map<std::string, OneDBCStringValue> _M_oneDBCMapStringValues;

    //! Result of the 2x2 linear system to be solved at each side
    OneDModelHandler::Vec2D _M_resBC;

    //! Reference to the solver current unknowns (U)
    const OneDModelHandler::ScalVec_vector& _M_U_thistime;
    //! Value of U at the boundary, at the neighboring internal node
    //    OneDModelHandler::Vec2D _M_line1, _M_line2;
    //! Boundary Dof (right or left)
    UInt _M_boundaryDof;
    //! Reference to the solver non linear flux functions
    const FLUX& _M_fluxFun;
};



template< class FLUX >
class OneDBCHandler
{
public:

    OneDBCHandler( const OneDModelHandler::ScalVec_vector& U_thistime,
                   const FLUX& fluxFun,
                   const Real& dimDof );

    void setBC( const OneDBCFunctionPointer& funptr,
                std::string const& border, std::string const& line, std::string const& var );

    void setBC( OneDBCFunctionPointer const& funptr,
                std::string const& border, std::string const& line, std::string const& var,
                OneDModelHandler::Vec2D const& matrixrow );

    template< class SOURCE >
    void setDefaultBC(const BasicOneDMesh& mesh,
                      const SOURCE& sourceFun,
                      const Real& dt);

    //! Apply boundary conditions
    void applyBC(const Real& time_val,
                 OneDModelHandler::Vec2D& left_BC_dir,
                 OneDModelHandler::Vec2D& right_BC_dir);

    // getters!!!
    inline OneDBC<FLUX>& BC( const std::string & bound ){ return *(_M_boundary[bound]); }

    OneDBCFunctionPointer& leftBCFunction( const std::string & line );
    OneDBCFunctionPointer& rightBCFunction( const std::string & line );

    inline bool& leftBCReady( const std::string & line ){ return _M_boundarybool["left"][line]; }
    inline bool& rightBCReady( const std::string & line ){ return _M_boundarybool["right"][line]; }

    inline void setBCLeft_internalnode(){ _M_boundary["left"]->isInternal()=true; }
    inline void setBCRight_internalnode(){ _M_boundary["right"]->isInternal()=true; }

private:
    //! Trick to use strings in C++ switch construct
    std::map<std::string, OneDBCStringValue> _M_oneDBCHandlerMapStringValues;
    //! Reference to the solver current unknowns (U)
    const OneDModelHandler::ScalVec_vector& _M_U_thistime;
    //! Reference to the solver non linear flux functions
    const FLUX& _M_fluxFun;

    std::map<std::string, boost::shared_ptr<OneDBC<FLUX> > > _M_boundary;

    std::map<std::string, std::map< std::string, bool > > _M_boundarybool;
};


template< class FLUX >
OneDBC<FLUX>::OneDBC( const OneDModelHandler::ScalVec_vector& U_thistime,
                      const FLUX& fluxFun,
                      const Real& dimDof,
                      const std::string& side )
    :
    _M_resBC(2),
    _M_U_thistime(U_thistime),
    _M_fluxFun(fluxFun)
{
    std::string prefix = "boundcond/" + side;
    Debug( 6311 ) << "[OneDBC::OneDBC] data file prefix = " << prefix << "\n";

    (side == "left") ? _M_boundaryDof=0 : _M_boundaryDof = dimDof - 1;

    _M_isInternal=false;

    _M_variable_at_line["first"]="not set";
    _M_variable_at_line["second"]="not set";

    _M_matrixrow_at_line["first"]=OneDModelHandler::Vec2D(2);
    _M_matrixrow_at_line["second"]=OneDModelHandler::Vec2D(2);

    _M_oneDBCMapStringValues["W1"] = OneDBCW1;
    _M_oneDBCMapStringValues["W2"] = OneDBCW2;
    _M_oneDBCMapStringValues["A"] = OneDBCA;
    _M_oneDBCMapStringValues["Q"] = OneDBCQ;
    _M_oneDBCMapStringValues["fun"] = OneDBCFUN;

}



template< class FLUX >
OneDModelHandler::Vec2D OneDBC<FLUX>::Uboundary(const OneDModelHandler::ScalVec& U1,
                                                const OneDModelHandler::ScalVec& U2) const
{
    OneDModelHandler::Vec2D Ubound(2);
    Ubound[0] = U1(_M_boundaryDof); Ubound[1] = U2(_M_boundaryDof);
    return Ubound;
}



template< class FLUX >
void OneDBC<FLUX>::compute_resBC(const Real& time_val)
{
    OneDModelHandler::Vec2D rhsBC(2);

    //! Eigen values of the jacobian diffFlux (= dF/dU = H)
    Real  eigval1, eigval2;
    //! Left eigen vectors for the eigen values eigval1 and eigval2
    OneDModelHandler::Vec2D left_eigvec1(2), left_eigvec2(2),
        left_eigvec_first(2), left_eigvec_second(2);

    OneDModelHandler::Vec2D U_boundary(2), W_boundary(2);

    for( UInt i=0; i<2; ++i ) {
        U_boundary[i] = _M_U_thistime[i](_M_boundaryDof);
        W_boundary[i] = _M_U_thistime[2+i](_M_boundaryDof);
    }

    _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary[0], U_boundary[1],
                                            eigval1, eigval2,
                                            left_eigvec1[0], left_eigvec1[1],
                                            left_eigvec2[0], left_eigvec2[1],
                                            _M_boundaryDof);

    Debug(6311) << "[OneDBC::compute_resBC_line] 1\n";
    Debug(6311) << "rhsBC.size() = " << rhsBC.size() << "\n";
    rhsBC[0]=_M_rhs_at_line["first"]->evaluate(time_val);
    rhsBC[1]=_M_rhs_at_line["second"]->evaluate(time_val);
    Debug(6311) << "rhsBC[0] = " << rhsBC[0] << "\n";;
    Debug(6311) << "rhsBC[1] = " << rhsBC[1] << "\n";;

    // this is not general
    // PRECONDITION (typical situation):
    //   on first line, left boundary, I impose W1
    //   on second line, left boundary, I impose W2
    //   on first line, right boundary, I impose W2
    //   on second line, right boundary, I impose W1
    // the code does not check for coherence (you cannot impose
    // the same variable on both lines!)
    Debug(6311) << "[OneDBC::compute_resBC_line] 2\n";
    _M_oneDBCMapStringValues[ _M_variable_at_line["first"] ] == OneDBCW1 ? //"W1"
        left_eigvec_first = left_eigvec1 :
        left_eigvec_first = left_eigvec2;

    Debug(6311) << "[OneDBC::compute_resBC_line] 3\n";
    _M_oneDBCMapStringValues[ _M_variable_at_line["second"] ] == OneDBCW1 ? //"W1"
        left_eigvec_second = left_eigvec1 :
        left_eigvec_second = left_eigvec2;

    compute_resBC_line("first", left_eigvec_first,
                       U_boundary, W_boundary, rhsBC[0]);
    compute_resBC_line("second", left_eigvec_second,
                       U_boundary, W_boundary, rhsBC[1]);

    Debug(6311) << "[OneDBC::compute_resBC] solving linear system with "
                << "\n\tfirst line = " << _M_matrixrow_at_line["first"][0]
                << ", " << _M_matrixrow_at_line["first"][1]
                << "\n\tsecond line = " << _M_matrixrow_at_line["second"][0]
                << ", " << _M_matrixrow_at_line["second"][1]
                << "\n\trhs = " << rhsBC[0]
                << ", " << rhsBC[1]
                << "\n";

    _M_resBC=_solveLinearSyst2x2(_M_matrixrow_at_line["first"],
                                 _M_matrixrow_at_line["second"], rhsBC);

}


template< class FLUX >
void OneDBC<FLUX>::compute_resBC_line( std::string line, OneDModelHandler::Vec2D left_eigvec,
                                       OneDModelHandler::Vec2D U, OneDModelHandler::Vec2D W, Real& rhs )
{
    ASSERT_PRE( left_eigvec.size() == 2 && U.size() == 2 && W.size() == 2,
                "compute_resBC_line works only for 2D vectors");

    Debug(6311) << "[OneDBC::compute_resBC_line] "
                << "_M_matrixrow_at_line.size() = "
                << _M_matrixrow_at_line.size() << ", "
                << "line = " << line << ", "
                << "_M_matrixrow_at_line[line].size() = "
                << _M_matrixrow_at_line[line].size() << ", "
                << "\n";

    Real LnUn, add;

    switch( _M_oneDBCMapStringValues[ _M_variable_at_line[line] ] )
        {
        case OneDBCW1: //"W1"
            //       ASSERT(eigval1<0. && eigval2<0.,
            // "The eigenvalues do no have the expected signs (lam1<0 and lam2<0).");
            _M_matrixrow_at_line[line] = left_eigvec;
            LnUn = dot( left_eigvec, U );
            add = LnUn - W[0];
            rhs += add;

            break;
        case OneDBCW2: //"W2"
            _M_matrixrow_at_line[line] = left_eigvec;
            LnUn = dot( left_eigvec, U );
            add = LnUn - W[1];
            rhs += add;
            break;
        case OneDBCA: //"A"
            _M_matrixrow_at_line[line][0] = 1.; _M_matrixrow_at_line[line][1] = 0.;
            break;
        case OneDBCQ: //"Q"
            _M_matrixrow_at_line[line][0] = 0.; _M_matrixrow_at_line[line][1] = 1.;
            break;
        case OneDBCFUN: //linear combination of A, Q
            break;
        default: std::cout << "\n[OneDBC::compute_resBC] Wrong boundary variable as " << line
                           << " condition at node " << _M_boundaryDof;
        }

    Debug(6311) << "[OneDBC::compute_resBC_line] to impose variable "
                << _M_variable_at_line[line]
                << ", " << line << " line = "
                << _M_matrixrow_at_line[line][0] << ", "
                << _M_matrixrow_at_line[line][1] << "\n";

}


template< class FLUX >
void OneDBC<FLUX>::applyBC(const Real& time_val,
                           OneDModelHandler::Vec2D& BC_dir)
{
    ASSERT_PRE( BC_dir.size() == 2,
                "applyBC works only for 2D vectors");

    if( _M_isInternal )
        {
            Debug(6311) << "[OneDBC::compute_resBC] found internal boundary\n";
        }
    else
        {
            compute_resBC(time_val);
            for( UInt i=0; i<2; ++i )
                BC_dir[i]=_M_resBC[i];
        }

    Debug(6311) << "[OneDBC::applyBC] at node " << _M_boundaryDof
                << " imposing [ A, Q ] = [ " << BC_dir[0]
                << ", " << BC_dir[1] << " ]\n";
}


/*! Matrix A is given by two pairs corresponding to the 2 _M_lines.
  A = [_M_matrixrow_at_line["first"];
  _M_matrixrow_at_line["second"] ]
  return A^{-1} * rhs2d
*/
template< class FLUX >
OneDModelHandler::Vec2D
OneDBC<FLUX>::_solveLinearSyst2x2(const OneDModelHandler::Vec2D& line1,
                                  const OneDModelHandler::Vec2D& line2,
                                  const OneDModelHandler::Vec2D& rhs2d) const
{
    ASSERT_PRE( line1.size() == 2 && line2.size() == 2 && rhs2d.size() == 2,
                "_solveLinearSyst2x2 works only for 2D vectors");

    long double aa11, aa12, aa21, aa22;

    aa11 = line1[0];  aa12 = line1[1];
    aa21 = line2[0];  aa22 = line2[1];

    long double determinant = aa11 * aa22 - aa12 * aa21;

    Debug(6311) << "[OneDBC::_solveLinearSyst2x2] solving linear system with "
                << "\n\tfirst line = " << aa11
                << ", " << aa12
                << "\n\tsecond line = " << aa21
                << ", " << aa22
                << "\n\trhs = " << rhs2d[0]
                << ", " << rhs2d[1]
                << "\n\tdet1 = " << aa11 * aa22
                << ", det2 = " << aa12 * aa21
                << "\n\tdet = " << determinant
                << "\n";

    ASSERT( determinant != 0,
            "Error: the 2x2 system on the boundary is not invertible."
            "\nCheck the boundary conditions.");

    OneDModelHandler::Vec2D res(2);
    res[0] = ( aa22 * rhs2d[0] - aa12 * rhs2d[1] ) / determinant;
    res[1] = ( - aa21 * rhs2d[0] + aa11 * rhs2d[1] ) / determinant;
    return res;
}



template< class FLUX >
OneDBCHandler<FLUX>::OneDBCHandler( const OneDModelHandler::ScalVec_vector& U_thistime,
                                    const FLUX& fluxFun,
                                    const Real& dimDof )
    :
    _M_U_thistime(U_thistime),
    _M_fluxFun(fluxFun)
{
    Debug( 6311 ) << "[OneDBCHandler::OneDBCHandler] Creating OneDBC classes.\n";

    _M_boundary["left"].reset( new OneDBC<FLUX>( U_thistime, /*W_thistime,*/ fluxFun, dimDof, "left" ) );
    _M_boundary["right"].reset( new OneDBC<FLUX>( U_thistime, /*W_thistime,*/ fluxFun, dimDof, "right" ) );

    _M_boundarybool["left"].insert( make_pair("first", false));
    _M_boundarybool["left"].insert( make_pair("second", false));
    _M_boundarybool["right"].insert( make_pair("first", false));
    _M_boundarybool["right"].insert( make_pair("second", false));

    _M_oneDBCHandlerMapStringValues["left"] = OneDBCLeftBoundary;
    _M_oneDBCHandlerMapStringValues["right"] = OneDBCRightBoundary;
    _M_oneDBCHandlerMapStringValues["first"] = OneDBCFirstRHS;
    _M_oneDBCHandlerMapStringValues["second"] = OneDBCSecondRHS;

}


template< class FLUX >
template< class SOURCE >
void OneDBCHandler<FLUX>::setDefaultBC( const BasicOneDMesh& mesh,
                                        const SOURCE& sourceFun,
                                        const Real& dt )
{
    if(!_M_boundarybool["left"].operator[]("first"))
        {
            OneDBCFunctionPointer point ( new Reimann<FLUX>( mesh, _M_fluxFun,
                                                             _M_U_thistime, //_M_W_thistime,
                                                             "left", "W1" ) );
            Debug( 6311 ) << "[OneDBCHandler::setDefaultBC] Invoking setBC.\n";
            setBC(point, "left", "first", "W1");
        }
    if(!_M_boundarybool["left"].operator[]("second"))
        {
            OneDBCFunctionPointer point ( new Compatibility<FLUX, SOURCE>( mesh,
                                                                           _M_fluxFun, sourceFun, _M_U_thistime, //_M_W_thistime,
                                                                           dt, "left", "W2" ) );
            Debug( 6311 ) << "[OneDBCHandler::setDefaultBC] Invoking setBC.\n";
            setBC(point, "left", "second", "W2");
        }
    if(!_M_boundarybool["right"]["first"])
        {
            OneDBCFunctionPointer point ( new Reimann<FLUX>( mesh, _M_fluxFun,
                                                             _M_U_thistime, //_M_W_thistime,
                                                             "right", "W2" ) );
            Debug( 6311 ) << "[OneDBCHandler::setDefaultBC] Invoking setBC.\n";
            setBC(point, "right", "first", "W2");
        }
    if(!_M_boundarybool["right"]["second"])
        {
            OneDBCFunctionPointer point ( new Compatibility<FLUX, SOURCE>( mesh,
                                                                           _M_fluxFun, sourceFun, _M_U_thistime, //_M_W_thistime,
                                                                           dt, "right", "W1" ) );
            Debug( 6311 ) << "[OneDBCHandler::setDefaultBC] Invoking setBC.\n";
            setBC(point, "right", "second", "W1");
        }

}

template< class FLUX >
void OneDBCHandler<FLUX>::applyBC(const Real& time_val,
                                  OneDModelHandler::Vec2D& left_BC_dir,
                                  OneDModelHandler::Vec2D& right_BC_dir)
{
    ASSERT_PRE( left_BC_dir.size() == 2 && right_BC_dir.size() == 2,
                "applyBC works only for 2D vectors");

    _M_boundary["left"]->applyBC( time_val, left_BC_dir );
    _M_boundary["right"]->applyBC( time_val, right_BC_dir );

    Debug(6311) << "[OneDBCHandler::applyBC] at left "
                << " imposing [ A, Q ] = [ " << left_BC_dir[0]
                << ", " << left_BC_dir[1] << " ]\n";
    Debug(6311) << "[OneDBCHandler::applyBC] at right "
                << " imposing [ A, Q ] = [ " << right_BC_dir[0]
                << ", " << right_BC_dir[1] << " ]\n";
};



template< class FLUX >
void OneDBCHandler<FLUX>::setBC( OneDBCFunctionPointer const& funptr,
                                 std::string const& border,
                                 std::string const& line,
                                 std::string const& var )
{
    _M_boundarybool[border][line] = true;
    _M_boundary[border]->variable(line)=var;
    _M_boundary[border]->rhs(line) = funptr;
    Debug( 6311 ) << "[OneDBCHandler::setBC] imposing function at "
                  << border << " boundary ("
                  << line << " line), variable "
                  << var << ".\n";
}



template< class FLUX >
void OneDBCHandler<FLUX>::setBC( OneDBCFunctionPointer const& funptr,
                                 std::string const& border,
                                 std::string const& line,
                                 std::string const& var,
                                 OneDModelHandler::Vec2D const& matrixrow )
{
    ASSERT_PRE( matrixrow.size() == 2,
                "setBC works only for 2D vectors");

    setBC(funptr, border, line, var);
    Debug( 6311 ) << "[OneDBCHandler::setBC] imposing matrix row as well.\n";
    _M_boundary[border]->matrixrow(line)=matrixrow;
}



template< class FLUX >
OneDBCFunctionPointer&
OneDBCHandler<FLUX>::leftBCFunction( const std::string& line )
{
    return _M_boundary["left"]->rhs(line);
}



template< class FLUX >
OneDBCFunctionPointer&
OneDBCHandler<FLUX>::rightBCFunction( const std::string& line )
{
    return _M_boundary["right"]->rhs(line);
}


}
#endif
