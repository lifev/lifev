/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
#include <life/lifefem/dof.hpp>

namespace LifeV
{
//! Constructor
Dof::Dof( const LocalDofPattern& _fe, UInt off ) : fe( _fe ), _offset( off ), _totalDof( 0 ),
        _nEl( 0 ), nlv( 0 ), nle( 0 ), nlf( 0 ), _ltg(),
        _numFaces(0),_ltgByFace(),_gtlByFace()
{
    //Getting the face
    switch( _fe.nbLocalDof() )
    {
        case 2:
            // No _fToP (it is 1D)
            _numLocalDofByFace = 1;
            break;
        case 4:
            _fToP = LinearTetra::fToP;
            _numLocalDofByFace = 3;
            break;
        case 5:
            _fToP = LinearTetraBubble::fToP;
            _numLocalDofByFace = 3;
            break;
        case 10:
            _fToP = QuadraticTetra::fToP;
            _numLocalDofByFace = 6;
            break;
        case 8:
            _fToP = LinearHexa::fToP;
            _numLocalDofByFace = 4;
            break;
        case 27:
            _fToP = QuadraticHexa::fToP;
            _numLocalDofByFace = 27;
            break;
        default:
            std::cout << "Warning: This refFE is not available for the dof by face." << std::endl;
            _numLocalDofByFace = 0;
            break;
    }

    for ( UInt i = 0; i < 5; ++i )
        _ncount[ i ] = 0;
}

//! Copy constructor
Dof::Dof( const Dof & dof2 ) : fe( dof2.fe ), _offset( dof2._offset ),
        _totalDof( dof2._totalDof ), _nEl( dof2._nEl ),
        nlv( dof2.nlv ), nle( dof2.nle ), nlf( dof2.nlf ),
        _ltg( dof2._ltg ),_numFaces(dof2._numFaces),
        _ltgByFace(dof2._ltgByFace),_gtlByFace(dof2._gtlByFace),
        _fToP(dof2._fToP),_numLocalDofByFace(dof2._numLocalDofByFace)
{
    if ( &dof2 == this )
        return ;

    //fe=dof2.fe;
    //  _offset=dof2._offset;
    //  _ltg=dof2._ltg;
    //  _totalDof=dof2._totalDof;
    for ( UInt i = 0; i < 5; ++i )
        _ncount[ i ] = dof2._ncount[ i ];
}

//! Ouput
void Dof::showMe( std::ostream & out, bool verbose ) const
{
    out << " Degree of Freedom (Dof) Object" << std::endl;
    out << " Total Dof Stored             " << _totalDof << std::endl;
    out << " With offset (min. Dof Id) =  " << _offset << std::endl;
    out << " Dof's on Vertices  from " << _ncount[ 0 ] << " , to:" << _ncount[ 1 ] - 1 << std::endl;
    out << " Dof's on Edges     from " << _ncount[ 1 ] << " , to:" << _ncount[ 2 ] - 1 << std::endl;
    out << " Dof's on Faces     from " << _ncount[ 2 ] << " , to:" << _ncount[ 3 ] - 1 << std::endl;
    out << " Dof's on Volumes   from " << _ncount[ 3 ] << " , to:" << _ncount[ 4 ] - 1 << std::endl;
    if ( verbose )
    {
        out << "************************************************************" << std::endl;
        out << "           Local to Global DOF table" << std::endl;
        out << "************************************************************" << std::endl;
        out << "Element Id   Loc. N.   Global N.  #  Element Id  Loc. N. Global N. " << std::endl;


        for ( UInt i = 0; i < _nEl;++i )
        {
            for ( UInt j = 0; j < numLocalDof();++j )
            {
                out.width( 10 );
                out << i + 1;
                out.width( 10 );
                out << j + 1;
                out.width( 10 );
                out << localToGlobal( i + 1, j + 1 );
                out << " # ";
                if ( (i*numLocalDof()+j) % 2 != 0 )
                    out << std::endl;
            }

        }out << std::endl;

    }

}

void Dof::showMeByFace(std::ostream& out, bool verbose) const{
  out << "--------------------------------------------------------------------------------" << std::endl;
  out << " Degree of freedom by face object " << std::endl;
  out << "--------------------------------------------------------------------------------" << std::endl;

  out << " Offset (min Dof Id) = " << _offset << std::endl;
  out << " Number of local dof per face = " << _numLocalDofByFace << std::endl;

  if(verbose){
    out << "*********************************************************************************" << std::endl;
    out << " Local-to-global DOF table (DOF grouped by internal face)" << std::endl;
    out << "*********************************************************************************" << std::endl;
    out << "=================================================================================" << std::endl;
    out << "Face ID     Local DOF   Global DOF  " << std::endl;
    out << "=================================================================================" << std::endl;

    for(UInt i = 0; i < _numFaces; ++i){
      for(UInt j = 0; j < _numLocalDofByFace; ++j){
    out.width(12);
    out << i + 1;
    out.width(12);
    out << j + 1;
    out.width(12);
    out << localToGlobal(i+1, j+1);
    out << " # ";
    if(j % 2 != 0) out << std::endl;
      } // for j
    } //for i
  } // if verbose
}

ID Dof::localToGlobalByFace(const ID& faceId, const ID& localDOF, bool& exist ) const
{
    ASSERT_PRE( (_numLocalDofByFace>0) , "This data are not available for this reference element");
    std::map<ID,ID>::const_iterator mapIt(_gtlByFace.find(faceId) );

    if(mapIt != _gtlByFace.end())
    {
        exist = true;
        //return _ltgByFace((*mapIt).second, localDOF);
        return _ltgByFace[(*mapIt).second][localDOF-1];
    }
    else
    {
        exist = false;
        return 0;
    }
}

//End of namespace LifeV
}
