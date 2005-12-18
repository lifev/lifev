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
#ifndef KNM_H_
#define KNM_H_

// une tentative qui ne marche pas
// de tableau constant
#include <cassert>
#include <iostream>
#include <iomanip>

#include <cmath>

//using namespace std;
#define const_R   R
#ifdef CHECK_KN
#include <cstdlib>
inline void Check_Kn( const char * str, const char * file, int line )
{
    std::cerr << "CHECK_KN: " << str << " in file: " << file << ", line " << line << std::endl;
    abort();
}
#define K_assert(i)  if (!(i)) Check_Kn(#i,__FILE__,__LINE__);

#else
#define K_assert(i)
#endif 
// version du 29 fev 2000
//  correction   for (... lj++,ui++  qui apelle le produit scalaire
//  petite correction  assert
// ajoute de operateur /= et *= sur des vecteurs
//   suppression de constructeur qui pose de probleme
//   correction oper +=  ...  dans RNM_op.h ligne 56  change = en oper
// version de 25 nov 99  sans const R
// ajoute de '.' pour extraire une colonne, ou ligne  , ...
//  version du 22 nov 1999   cast to KN_<const R>
//   version du 21 nov 1999  correction delete
//  version du 13 nov 1999
//  version du 18 mars 99
//  F. Hecht
// attention les indexations les indexations peuvent changer
//  puisque que l'on peut prendre la transposer d'une matrice
// tableau
// mais ils partent de 0
// version corrigee du 15/11/98
// version avec sous tableau  ---  mars 99
// -------
//  remarque du 8 mars 99 FH
// class pour prendre des sous-tableau
// attention aux PB de continute dans les tableaux
// on a supposer que les tableaux multi indices pouvait est vue comme
// un tableau continue ce qui est generalement faux quand l'on en
// prend un sous tableau
//   exemple: un tableau 3,5 est numerote comme:
//    0  3  6  9 12
//    1  4  7 10 13
//    2  5  8 11 14
//             step
//   indexi  n 1
//   indexj  m n
//   est le sous tableau  3,3  n'est pas numeroter consecutivement
//
//    Donc la fonction  IsVector1() nous dit si un tableau
//    ý un 2 ou 3 indices est ou non consecutif en memoire
//
//  ----------------------------------
//   version du 21 novembre 2000 FH
//   ---  initialisation --
namespace LifeV
{

template <class R>
class KNMKL_ ;
template <class R>
class KNMK_ ;
template <class R>
class KNM_ ;
template <class R>
class KN_ ;

template <class R>
class KNMKL ;
template <class R>
class KNMK ;
template <class R>
class KNM ;
template <class R>
class KN ;

template <class R>
class Add_KN_;
template <class R>
class Sub_KN_;
template <class R>
class Mulc_KN_;
template <class R>
class Add_Mulc_KN_;
template <class R>
class Mul_KNM_KN_;
template <class R>
class MatriceCreuseMulKN_;
template <class R>
class MatriceCreuseDivKN_;

class ShapeOfArray;

class FromTo
{
public:
    int from, to;
    FromTo( int i, int j ) : from( i ), to( j )
    {
        K_assert( i < j );
    }
};

class SubArray
{
public:
    const int n, step, start;
    //  SubArray(char  nn): n(-1),step(1),start(0) {}
    explicit SubArray( int nn, int sta = 0, int s = 1 ) : n( nn ), step( s ), start( sta )
    {}
    SubArray( const FromTo& ft ) : n( ft.to - ft.from + 1 ), step( 1 ), start( ft.from )
    {}
    SubArray( const ShapeOfArray & ); // all
    int end() const
    {
        return start + step * n;
    }
    int last() const
    {
        return start + step * ( n - 1 );
    }
    int len1() const
    {
        return step * ( n - 1 );
    }
};


class ShapeOfArray
{
protected:
public:

    const int n;     //   n  nb of item
    const int step;  //   step  nb of between 2 item
    const int next;  //  the   next array of same type in matrix for subarray
    // by default  no next
    ShapeOfArray( const ShapeOfArray & s, int nn ) : n( s.n ), step( s.n ), next( nn )
    {}
    ShapeOfArray( int nn ) : n( nn ), step( 1 ), next( -1 )
    {}

    ShapeOfArray( int nn, int s ) : n( nn ), step( s ), next( -1 )
    {}

    ShapeOfArray( int nn, int s, int nextt ) : n( nn ), step( s ), next( nextt )
    {}

    ShapeOfArray( const ShapeOfArray &old, const SubArray &sub )
            : n( sub.n ), step( old.step*sub.step ), next( old.next )
    {
        K_assert( ( sub.last() ) * old.step <= old.last() );
    } // a constructor

    ShapeOfArray( const ShapeOfArray &old, int stepo, int start )
            : n( old.n - start ), step( old.step*stepo ), next( old.next )
    {
        K_assert( n >= 0 );
    }

    int end() const
    {
        return n * step;
    }
    int last() const
    {
        return ( n -1 ) * step;
    }
    int constant() const
    {
        return step == 0;
    }
    int index( int k ) const
    {
        K_assert( ( k >= 0 ) && ( ( k < n ) || !step ) );
        return step*k;
    }
    ShapeOfArray operator*( int stepp ) const
    {
        return ShapeOfArray( n, step * stepp, next );
    }
    bool SameShape( const ShapeOfArray & a ) const
    {
        return !step || !a.step || a.n == n ;
    }
    int N( const ShapeOfArray & a )
    {
        return step ? n : a.n;
    } // size of 2 shape


    // protected:
    int operator[] ( int k ) const
    {
        K_assert( ( k>=0 ) && ( ( k < n ) || !step ) );
        return step*k;
    }

};

std::ostream & operator<<( std::ostream & f, const ShapeOfArray & s );

inline bool SameShape( const ShapeOfArray & a, const ShapeOfArray & b )
{
    return !a.step || !b.step || a.n == b.n ;
}

inline int N( const ShapeOfArray & a, const ShapeOfArray & b )
{
    K_assert( SameShape( a, b ) );
    return a.step ? a.n : b.n ;
}

inline SubArray::SubArray( const ShapeOfArray & s )
        : n( s.n ), step( s.step ), start( 0 )
{}




template <class R>
std::ostream & operator<<( std::ostream & f, const KN_<const_R> & v ) ;


template <class R>
class KN_: public ShapeOfArray
{
protected:
    R *v;
public:
    int N() const
    {
        return this->n;
    }
    int size() const
    {
        return this->step ? N() * this->step : N();
    }
    operator R *() const
    {
        return this->v;
    }
    KN_( const KN_<R> & u ) : ShapeOfArray( u ), v( u.v )
    {}
    KN_( const KN_<R> & U, const SubArray & sa ) : ShapeOfArray( U, sa ), v( U.v + U.index( sa.start ) )
    {}

    KN_ operator() ( const SubArray & sa ) const
    {
        return KN_( *this, sa );
    } // sub array
    R & operator[] ( int i ) const
    {
        return v[ index( i ) ];
    }
    R & operator() ( int i ) const
    {
        return v[ index( i ) ];
    }

    R operator,( const KN_<const_R> & v ) const; // dot  product

    const KN_& operator =( const KN_<const_R> & u ) ;
    const KN_& operator +=( const KN_<const_R> & u ) ;
    const KN_& operator -=( const KN_<const_R> & u ) ;

    const KN_& operator *=( const KN_<const_R> & u ) ;
    const KN_& operator /=( const KN_<const_R> & u ) ;


    const KN_& operator = ( const_R a ) ;
    const KN_& operator +=( const_R a ) ;
    const KN_& operator -=( const_R a ) ;
    const KN_& operator /=( const_R a ) ;
    const KN_& operator *=( const_R a ) ;

    R KNMmin() const ; // min -> KNMmin to avoid clash name with aztec, JFG
    R KNMmax() const ; // max -> KNMmax to avoid clash name with aztec, JFG
    R sum() const ;
    KN_ & map( R ( * ) ( R ) );

    const KN_& operator =( const Add_KN_<R> & u ) ;
    const KN_& operator+=( const Add_KN_<R> & u ) ;
    const KN_& operator-=( const Add_KN_<R> & u ) ;

    const KN_& operator =( const Sub_KN_<R> & u ) ;
    const KN_& operator-=( const Sub_KN_<R> & u ) ;
    const KN_& operator+=( const Sub_KN_<R> & u ) ;

    const KN_& operator =( const Mulc_KN_<R> & u ) ;
    const KN_& operator+=( const Mulc_KN_<R> & u ) ;
    const KN_& operator-=( const Mulc_KN_<R> & u ) ;

    const KN_& operator =( const Add_Mulc_KN_<R> & u ) ;
    const KN_& operator+=( const Add_Mulc_KN_<R> & u ) ;
    const KN_& operator-=( const Add_Mulc_KN_<R> & u ) ;

    const KN_& operator =( const Mul_KNM_KN_<R> & u ) ;
    const KN_& operator+=( const Mul_KNM_KN_<R> & u ) ;
    const KN_& operator-=( const Mul_KNM_KN_<R> & u ) ;

    const KN_& operator =( const MatriceCreuseMulKN_<R> & ) ;
    const KN_& operator =( const MatriceCreuseDivKN_<R> & ) ;

    friend std::ostream & operator<< <R>( std::ostream & f, const KN_<const_R> & v ) ;

private:

    KN_& operator++()
    {
        K_assert( next >= 0 );
        v += next;
        return *this;
    } //    ++U
    KN_& operator--()
    {
        K_assert( next >= 0 );
        v -= next;
        return *this;
    } //    --U
    KN_ operator++( int )
    {
        K_assert( next >= 0 );
        KN_ old = *this;
        v = v + next;
        return old;
    } // U++
    KN_ operator--( int )
    {
        K_assert( next >= 0 );
        KN_ old = *this;
        v = v - next;
        return old;
    } // U++

    KN_( R *u, const ShapeOfArray & s ) : ShapeOfArray( s ), v( u )
    {}
    KN_( R *u, int nn, int s ) : ShapeOfArray( nn, s ), v( u )
    {}
    KN_( R *u, int nn, int s, int nextt ) : ShapeOfArray( nn, s, nextt ), v( u )
    {}
    KN_( R *u, int nn ) : ShapeOfArray( nn ), v( u )
    {}
    KN_( const KN_<R> & u, int offset ) : ShapeOfArray( u ), v( &u[ offset ] )
    {}
    KN_( const KN_<R> & u, const ShapeOfArray &sh, int startv = 0 )
            : ShapeOfArray( sh*u.step ), v( &u[ startv ] )
    {}
    KN_( const KN_<R> & u, int nnext, const ShapeOfArray &sh, int startv = 0 )
            : ShapeOfArray( sh.n, sh.step*u.step, nnext ), v( &u[ startv ] )
    { }


    // friend class KN_<R>;
    friend class KNM_<R>;
    friend class KNMK_<R>;
    friend class KNMKL_<R>;
    friend class KN<R>;
    friend class KNM<R>;
    friend class KNMK<R>;
    friend class KNMKL<R>;

};

template <class R>
class KNM_: public KN_<R>
{
public:
    ShapeOfArray shapei;
    ShapeOfArray shapej;
public:
    int IsVector1() const
    {
        return ( shapei.n * shapej.n ) == this->n ;
    }
    int N() const
    {
        return shapei.n;
    }
    int M() const
    {
        return shapej.n;
    }
    int size() const
    {
        return shapei.n * shapej.n;
    }

    KNM_( R* u, const ShapeOfArray & s,
          const ShapeOfArray & si,
          const ShapeOfArray & sj )
            : KN_<R>( u, s ), shapei( si ), shapej( sj )
    {}
    KNM_( R* u, int n, int m )
            : KN_<R>( u, ShapeOfArray( n*m ) ), shapei( n, 1, n ), shapej( m, n, 1 )
    {}
    KNM_( R* u, int n, int m, int s )
            : KN_<R>( u, ShapeOfArray( n*m, s ) ), shapei( n, 1, n ), shapej( m, n, 1 )
    {}
    KNM_( KN_<R> u, int n, int m )
            : KN_<R>( u, ShapeOfArray( m*n ) ), shapei( n, 1, n ), shapej( m, n, 1 )
    { }

    KNM_( const KN_<R> &u, const ShapeOfArray & si, const ShapeOfArray & sj, int offset = 0 )
            : KN_<R>( &u[ offset ], si.last() + sj.last() + 1, u.step ), shapei( si ), shapej( sj )
    {
        K_assert( offset >= 0 && this->n + ( this->v - ( R* ) u ) <= u.n );
    }
    KNM_( const KN_<R> &u, const ShapeOfArray & si, const ShapeOfArray & sj, int offset, int nnext )
            : KN_<R>( &u[ offset ], si.last() + sj.last() + 1, u.step, nnext ), shapei( si ), shapej( sj )
    {
        K_assert( offset >= 0 && this->n + ( this->v - ( R* ) u ) <= u.n );
    }

    KNM_( KNM_<R> U, const SubArray & si, const SubArray & sj )
            : KN_<R>( U, SubArray( U.ij( si.len1(), sj.len1() ) + 1, U.ij( si.start, sj.start ) ) ),
            shapei( U.shapei, si ), shapej( U.shapej, sj )
    {}

    KNM_( KNM_<R> U, const SubArray & sa, const SubArray & si, const SubArray & sj )
            : KN_<R>( U, SubArray( sa ) ), shapei( U.shapei, si ), shapej( U.shapej, sj )
    {}

    KNM_( const KNM_<R> & u )
            : KN_<R>( u ), shapei( u.shapei ), shapej( u.shapej )
    {}

    KNM_ operator() ( const SubArray & sa, const SubArray & sb ) const
    {
        return KNM_( *this, sa, sb );
    } // sub array

    int ij( const int i, const int j ) const
    {
        return shapei.index( i ) + shapej.index( j );
    }
    int indexij( int i, int j ) const
    {
        return this->index( shapei.index( i ) + shapej.index( j ) );
    }
    R & operator() ( int i,int j ) const
    {
        return this->v[ indexij( i, j ) ];
    }
    //Alain (28/06/02): version for unsigned int.
    unsigned int indexij( unsigned int i, unsigned int j ) const
    {
        return this->index( shapei.index( i ) + shapej.index( j ) );
    }
    R & operator() ( unsigned int i,unsigned int j ) const
    {
        return this->v[ indexij( i, j ) ];
    }
    //Alain (18/10/02): version for long unsigned int.
    long unsigned int indexij( long unsigned int i, long unsigned int j ) const
    {
        return this->index( shapei.index( i ) + shapej.index( j ) );
    }
    R & operator() ( long unsigned int i,long unsigned int j ) const
    {
        return this->v[ indexij( i, j ) ];
    }
    //

    KN_<R> operator() ( const char,int j ) const   // une colonne j  ('.',j)
    {
        return KN_<R>( &this->v[ this->index( shapej.index( j ) ) ], shapei * this->step );
    }
    KN_<R> operator() ( int i ,const char ) const   // une ligne i  (i,'.')
    {
        return KN_<R>( &this->v[ this->index( shapei.index( i ) ) ], shapej * this->step );
    }
    KN_<R> operator() ( const char,const char ) const   // tous
    {
        return * this;
    }

    // Alain (21/11/01) : modification of t(), constructor
    // KNM_<R>(*this,shapej,shapei,v) was not defined.
    KNM_<R> t() const
    {
        return KNM_<R>( this->v, *this, shapej, shapei );
    }

    const KNM_& operator =( const KNM_<const_R> & u ) ;
    const KNM_& operator =( const_R a ) ;
    const KNM_& operator+=( const_R a ) ;
    const KNM_& operator-=( const_R a ) ;
    const KNM_& operator/=( const_R a ) ;
    const KNM_& operator*=( const_R a ) ;
    const KNM_& operator+=( const KNM_<const_R> & u ) ;
    const KNM_& operator-=( const KNM_<const_R> & u ) ;
    const KNM_& operator*=( const KNM_<const_R> & u ) ;
    const KNM_& operator/=( const KNM_<const_R> & u ) ;

private:
    KNM_& operator++()
    {
        this->v += this->next;
        return *this;
    } // ++U
    KNM_& operator--()
    {
        this->v -= this->next;
        return *this;
    } // ++U
    KNM_ operator++( int )
    {
        KNM_<R> old = *this;
        this->v = this->v + this->next;
        return old;
    } // U++
    KNM_ operator--( int )
    {
        KNM_<R> old = *this;
        this->v = this->v - this->next;
        return old;
    } // U--


    friend class KN_<R>;
    // friend class KNM_<R>;
    friend class KNMK_<R>;
    friend class KNMKL_<R>;
    friend class KN<R>;
    friend class KNM<R>;
    friend class KNMK<R>;
    friend class KNMKL<R>;
};




template <class R>
class KNMK_: public KN_<R>
{
    friend class KNMK<R>;
public:
    ShapeOfArray shapei;
    ShapeOfArray shapej;
    ShapeOfArray shapek;
public:
    int IsVector1() const
    {
        return ( shapei.n * shapej.n * shapek.n ) == this->n ;
    }
    int N() const
    {
        return shapei.n;
    }
    int M() const
    {
        return shapej.n;
    }
    int K() const
    {
        return shapek.n;
    }
    int size() const
    {
        return shapei.n * shapej.n * shapek.n;
    }
    KNMK_( const ShapeOfArray & s,
           const ShapeOfArray & si,
           const ShapeOfArray & sj,
           const ShapeOfArray & sk,
           R * u )
            : KN_<R>( u, s ), shapei( si ), shapej( sj ), shapek( sk )
    {}

    KNMK_( R* u, const int n, const int m, const int k )
            : KN_<R>( u, ShapeOfArray( n*m*k ) ), shapei( n, 1, n ), shapej( m, n, 1 ), shapek( k, n*m, n*m )
    {}
    ;

    //  KNMK_(const KN_<R> & u,int n,int m,int k)
    //   : KN_<R>(ShapeOfArray(n*m*k)),shapei(n,1,n),shapekj(m,n,1),u),
    //     shapek(k,n*m,n*m){};

    KNMK_( const KNMK_<R> &U, const SubArray & si, const SubArray & sj, const SubArray & sk ) :
            KN_<R>( U, SubArray( U.ijk( si.len1(), sj.len1(), sk.len1() ) + 1,
                                 U.ijk( si.start, sj.start, sk.start ) ) ),
            shapei( U.shapei, si ),
            shapej( U.shapej, sj ),
            shapek( U.shapek, sk )
    {}

    KNMK_( const KNMK_<R> & u ) : KN_<R>( u ), shapei( u.shapei ), shapej( u.shapej ), shapek( u.shapek )
    {}


    int ijk( const int i, const int j, const int k ) const
    {
        return shapei.index( i ) + shapej.index( j ) + shapek.index( k );
    }
    int indexijk( int i, int j, int k ) const
    {
        return this->index( shapei.index( i ) + shapej.index( j ) + shapek.index( k ) );
    }

    R & operator() ( int i,int j, int k ) const
    {
        return this->v[ indexijk( i, j, k ) ];
    }

    //  pas de tableau suivant
    KN_<R> operator() ( const char ,int j, int k ) const
    { // le tableau (.,j,k)
        return KN_<R>( *this, -1, shapei, shapej[ j ] + shapek[ k ] );
    }
    KN_<R> operator() ( int i,const char , int k ) const
    { // le tableau (i,.,k)
        return KN_<R>( *this, -1, shapej, shapei[ i ] + shapek[ k ] );
    }
    KN_<R> operator() ( int i,int j, const char ) const
    { // le tableau (i,j,.)
        return KN_<R>( *this, -1, shapek, shapei[ i ] + shapej[ j ] );
    }
    //
    KNM_<R> operator() ( const char ,const char , int k ) const
    { // le tableau (.,.,k)
        return KNM_<R>( *this, shapei, shapej, shapek[ k ], shapek.next );
    } // step = n*m
    //attention les suivants ne marche pas
    KNM_<R> operator() ( const char ,int j, const char ) const
    { // le tableau (.,j,.)
        return KNM_<R>( *this, shapei, shapek, shapej[ j ], -1 /*shapej.next*/ );
    } // step = n

    KNM_<R> operator() ( int i,const char , const char ) const
    { // le tableau (i,.,.)
        return KNM_<R>( *this, shapej, shapek, shapei[ i ], -1 /*shapei.next*/ );
    }  // step = 1

    const KNMK_& operator =( const KNMK_<const_R> & u ) ;
    const KNMK_& operator+=( const KNMK_<const_R> & u ) ;
    const KNMK_& operator-=( const KNMK_<const_R> & u ) ;
    const KNMK_& operator/=( const KNMK_<const_R> & u ) ;
    const KNMK_& operator*=( const KNMK_<const_R> & u ) ;
    const KNMK_& operator =( const_R a ) ;
    const KNMK_& operator+=( const_R a ) ;
    const KNMK_& operator-=( const_R a ) ;
    const KNMK_& operator/=( const_R a ) ;
    const KNMK_& operator*=( const_R a ) ;

    KNMK_ operator() ( SubArray si,SubArray sj, SubArray sk ) const
    {
        return KNMK_( *this, si, sj, sk );
    }

private:
    //  KNMK_&  operator++(){v += next;return *this;} // ++U
    //  KNMK_&  operator--(){v -= next;return *this;} // --U
    //  KNMK_  operator++(int ){KNMK_ old=*this;v = v +next;return old;} // U++
    //  KNMK_  operator--(int ){KNMK_ old=*this;v = v -next;return old;} // U--


    friend class KNM_<R>;
    friend class KN_<R>;

};

//===========================================================================
//===========================================================================
//===========================================================================
// A. Veneziani: KNMKL_ class for handling 4 indexes - October 2002

template <class R>
class KNMKL_: public KN_<R>
{
    friend class KNMKL<R>;
public:
    ShapeOfArray shapei;
    ShapeOfArray shapej;
    ShapeOfArray shapek;
    ShapeOfArray shapel;
public:
    int IsVector1() const
    {
        return ( shapei.n * shapej.n * shapek.n * shapel.n ) == this->n ;
    }
    int N() const
    {
        return shapei.n;
    }
    int M() const
    {
        return shapej.n;
    }
    int K() const
    {
        return shapek.n;
    }
    int L() const
    {
        return shapel.n;
    }
    int size() const
    {
        return shapei.n * shapej.n * shapek.n * shapel.n;
    }
    KNMKL_( const ShapeOfArray & s,
            const ShapeOfArray & si,
            const ShapeOfArray & sj,
            const ShapeOfArray & sk,
            const ShapeOfArray & sl,
            R * u )
            : KN_<R>( u, s ), shapei( si ), shapej( sj ), shapek( sk ), shapel( sl )
    {}

    KNMKL_( R* u, const int n, const int m, const int k, const int l )
            : KN_<R>( u, ShapeOfArray( n*m*k*l ) ), shapei( n, 1, n ), shapej( m, n, 1 ), shapek( k, n*m, n*m ), shapel( l, n*m*k, n*m*k )
    {}
    ;


    KNMKL_( const KNMK_<R> &U, const SubArray & si, const SubArray & sj, const SubArray & sk, const SubArray & sl ) :
            KN_<R>( U, SubArray( U.ijkl( si.len1(), sj.len1(), sk.len1(), sl.len1() ) + 1,
                                 U.ijkl( si.start, sj.start, sk.start, sl.start ) ) ),
            shapei( U.shapei, si ),
            shapej( U.shapej, sj ),
            shapek( U.shapek, sk ),
            shapel( U.shapel, sl )
    {}

    KNMKL_( const KNMK_<R> & u ) : KN_<R>( u ), shapei( u.shapei ), shapej( u.shapej ), shapek( u.shapek ), shapel( u.shapel )
    {}


    int ijkl( const int i, const int j, const int k, const int l ) const
    {
        return shapei.index( i ) + shapej.index( j ) + shapek.index( k ) + shapel.index( l );
    }
    int indexijk( int i, int j, int k, int l ) const
    {
        return this->index( shapei.index( i ) + shapej.index( j ) + shapek.index( k ) + shapel.index( l ) );
    }

    R & operator() ( int i,int j, int k, int l ) const
    {
        return this->v[ indexijk( i, j, k, l ) ];
    }

    //  pas de tableau suivant
    KN_<R> operator() ( const char ,int j, int k, int l ) const
    { // le tableau (.,j,k,l)
        return KN_<R>( *this, -1, shapei, shapej[ j ] + shapek[ k ] + shapel[ l ] );
    }

    KN_<R> operator() ( int i,const char , int k, int l ) const
    { // le tableau (i,.,k,l)
        return KN_<R>( *this, -1, shapej, shapei[ i ] + shapek[ k ] + shapel[ l ] );
    }

    KN_<R> operator() ( int i,int j, const char, int l ) const
    { // le tableau (i,j,.,l)
        return KN_<R>( *this, -1, shapek, shapei[ i ] + shapej[ j ] + shapel[ l ] );
    }

    KN_<R> operator() ( int i,int j, int k, const char ) const
    { // le tableau (i,j,k,.)
        return KN_<R>( *this, -1, shapel, shapei[ i ] + shapej[ j ] + shapek[ k ] );
    }


    // ATTENTION: This part has to be verified
    //  KNM_<R>  operator()(const char ,const char ,int k, int l)  const  { // le tableau (.,.,k,l)
    //    return KNM_<R>(*this,shapei,shapej,shapek[k]+shapel[l],-1/*shapel.next*/);} // step = k*n*m
    //  KNM_<R>  operator()(const char ,int j,const char , int l)  const  { // le tableau (.,j,.,l)
    //    return KNM_<R>(*this,shapei,shapek,shapej[j]+shapel[l],-1/*shapel.next*/);} // step = k*n*m
    //  KNM_<R>  operator()(const char ,int j, int k, const char)  const  { // le tableau (.,j,k,.)
    //    return KNM_<R>(*this,shapei,shapel,shapej[j]+shapek[k],-1/*shapek.next*/);} // step = n*m
    //  //attention les suivants ne marche pas
    //  KNM_<R>  operator()(int i, const char ,const char, int l)  const  { // le tableau (i,.,.,l)
    //         return KNM_<R>(*this,shapei,shapel,shapej[j]+shapek[k],-1/*shapej.next*/);} // step = n
    //  KNM_<R>  operator()(int i, const char ,int k,const char )  const  { // le tableau (i,.,k,.)
    //         return KNM_<R>(*this,shapei,shapek,shapej[j]+shapel[l],-1/*shapel.next*/);} // step = n
    //  KNM_<R>  operator()(int i, int j, const char ,const char )  const  { // le tableau (i,j,.,.)
    //         return KNM_<R>(*this,shapei,shapej,shapek[k]+shapel[l],-1/*1*/);}  // step = 1

    // //
    //  KNMK_<R>  operator()(const char ,const char ,int k)  const  { // le tableau (.,.,.,l)
    //         return KNMK_<R>(*this,shapei,shapej,shapek,shapel[l],shapel.next);} // step = k*n*m
    //  KNMK_<R>  operator()(const char ,const char ,int k)  const  { // le tableau (.,.,k,.)
    //         return KNMK_<R>(*this,shapei,shapej,shapel,shapek[k],shapek.next);} // step = n*m
    //  //attention les suivants ne marche pas
    //  KNMK_<R>  operator()(const char ,int j,const char )  const  { // le tableau (.,j,.,.)
    //         return KNMK_<R>(*this,shapei,shapek,shapel,shapej[j],-1/*shapej.next*/);} // step = n
    //  KNMK_<R>  operator()(int i,const char ,const char )  const  { // le tableau (i,.,.,.)
    //         return KNMK_<R>(*this,shapej,shapek,shapei,shapei[i],-1/*shapei.next*/);}  // step = 1

    const KNMKL_& operator =( const KNMKL_<const_R> & u ) ;
    const KNMKL_& operator+=( const KNMKL_<const_R> & u ) ;
    const KNMKL_& operator-=( const KNMKL_<const_R> & u ) ;
    const KNMKL_& operator/=( const KNMKL_<const_R> & u ) ;
    const KNMKL_& operator*=( const KNMKL_<const_R> & u ) ;
    const KNMKL_& operator =( const_R a ) ;
    const KNMKL_& operator+=( const_R a ) ;
    const KNMKL_& operator-=( const_R a ) ;
    const KNMKL_& operator/=( const_R a ) ;
    const KNMKL_& operator*=( const_R a ) ;

    KNMKL_ operator() ( SubArray si,SubArray sj, SubArray sk ) const
    {
        return KNMKL_( *this, si, sj, sk );
    }

private:
    //  KNMK_&  operator++(){v += next;return *this;} // ++U
    //  KNMK_&  operator--(){v -= next;return *this;} // --U
    //  KNMK_  operator++(int ){KNMK_ old=*this;v = v +next;return old;} // U++
    //  KNMK_  operator--(int ){KNMK_ old=*this;v = v -next;return old;} // U--

    friend class KNMK_<R>;
    friend class KNM_<R>;
    friend class KN_<R>;

};


//===========================================================================
//===========================================================================
//===========================================================================

template <class R>
class KN : public KN_<R>
{
public:

    // explicit  KN(const R & u):KN_<R>(new R(uu),1,0) {}
    KN( const int nn ) : KN_<R>( new R[ nn ], nn )
    {}
    KN( const int nn, R ( *f ) ( int i ) ) : KN_<R>( new R[ nn ], nn )
    {
        for ( int i = 0;i < this->n;i++ )
            this->v[ i ] = f( i );
    }
    KN( const int nn, const R & a ) : KN_<R>( new R[ nn ], nn )
    {
        KN_<R>::operator=( a );
    }
    KN( const int nn, int s, const R a ) : KN_<R>( new R[ nn ], nn, s )
    {
        KN_<R>::operator=( a );
    }
    template <class S>
    KN( const KN_<S> & s ) : KN_<R>( new R[ s.n ], s.n )
    {
        for ( int i = 0;i < this->n;i++ )
            this->v[ i ] = s[ i ];
    }
    template <class S>
    KN( const KN_<S> & s, R ( *f ) ( S ) ) : KN_<R>( new R[ s.n ], s.n )
    {
        for ( int i = 0;i < this->n;i++ )
            this->v[ i ] = f( s[ i ] );
    }
    explicit KN( const KN<R> & u ) : KN_<R>( new R[ u.n ], u.n )
    {
        KN_<R>::operator=( u );
    }
    explicit KN( const KN_<R> & u ) : KN_<R>( new R[ u.n ], u.n )
    {
        KN_<R>::operator=( u );
    }

    ~KN()
    {
        // should use boost::shared_ptr
        delete [] this->v;
    }

    const KN& operator =( const_R a )
    {
        KN_<R>::operator= ( a );
        return *this;
    }
    const KN& operator =( const KN_<R>& a )
    {
        KN_<R>::operator= ( a );
        return *this;
    }

    const KN& operator =( const KN<R>& a )
    {
        KN_<R>::operator= ( a );
        return *this;
    }

    const KN& operator =( const Add_KN_<R> & u )
    {
        KN_<R>::operator=( u );
        return *this;
    }
    const KN& operator =( const Sub_KN_<R> & u )
    {
        KN_<R>::operator=( u );
        return *this;
    }
    const KN& operator =( const Mulc_KN_<R> & u )
    {
        KN_<R>::operator=( u );
        return *this;
    }
    const KN& operator =( const Add_Mulc_KN_<R> & u )
    {
        KN_<R>::operator=( u );
        return *this;
    }
    const KN& operator =( const Mul_KNM_KN_<R> & u )
    {
        KN_<R>::operator=( u );
        return *this;
    }
    const KN& operator =( const MatriceCreuseMulKN_<R> & A )
    {
        KN_<R>::operator=( A );
        return *this;
    }
    const KN& operator =( const MatriceCreuseDivKN_<R> & A )
    {
        KN_<R>::operator=( A );
        return *this;
    }

    const KN& operator -=( const_R a )
    {
        KN_<R>::operator-=( a );
        return *this;
    }
    const KN& operator -=( const KN_<R>& a )
    {
        KN_<R>::operator-= ( a );
        return *this;
    }
    const KN& operator -=( const Add_KN_<R> & u )
    {
        KN_<R>::operator-=( u );
        return *this;
    }
    const KN& operator -=( const Sub_KN_<R> & u )
    {
        KN_<R>::operator-=( u );
        return *this;
    }
    const KN& operator -=( const Mulc_KN_<R> & u )
    {
        KN_<R>::operator-=( u );
        return *this;
    }
    const KN& operator -=( const Add_Mulc_KN_<R> & u )
    {
        KN_<R>::operator-=( u );
        return *this;
    }
    const KN& operator -=( const Mul_KNM_KN_<R> & u )
    {
        KN_<R>::operator-=( u );
        return *this;
    }

    const KN& operator +=( const_R a )
    {
        KN_<R>::operator += ( a );
        return *this;
    }
    const KN& operator += ( const KN_<R>& a )
    {
        KN_<R>::operator+= ( a );
        return *this;
    }
    const KN& operator +=( const Add_KN_<R> & u )
    {
        KN_<R>::operator+=( u );
        return *this;
    }
    const KN& operator +=( const Sub_KN_<R> & u )
    {
        KN_<R>::operator+=( u );
        return *this;
    }
    const KN& operator +=( const Mulc_KN_<R> & u )
    {
        KN_<R>::operator+=( u );
        return *this;
    }
    const KN& operator +=( const Add_Mulc_KN_<R> & u )
    {
        KN_<R>::operator+=( u );
        return *this;
    }
    const KN& operator +=( const Mul_KNM_KN_<R> & u )
    {
        KN_<R>::operator+=( u );
        return *this;
    }


    const KN& operator/=( const_R a )
    {
        KN_<R>::operator/=( a );
        return *this;
    }
    const KN& operator*=( const_R a )
    {
        KN_<R>::operator*=( a );
        return *this;
    }
    const KN& operator/=( const KN_<const_R>& a )
    {
        KN_<R>::operator/= ( a );
        return *this;
    }
    const KN& operator*=( const KN_<const_R>& a )
    {
        KN_<R>::operator*= ( a );
        return *this;
    }


    //  two opertor to cast to an array of constant
    //    operator KN_<const_R> & ()
    //          { return *  (KN_<const_R>*) this;}
    //    operator KN_<const_R> const & ()  const
    //          { return *(const KN_<const_R>*) this;}
    //    operator KN<const_R> & ()
    //          { return   (KN<const_R> &) *this;}
    //    operator KN<const_R> const & ()  const
    //          { return (const KN<const_R>& ) *this;}
};

//  Array with 2 indices
//  ---------------------

template <class R>
class KNM: public KNM_<R>
{
public:

    KNM( const int n, const int m )
            : KNM_<R>( new R[ n*m ], n, m )
    {
        assert( this->v != 0 );
    }

    /* Alain (28/06/02): I remove the explicit statment for allowing implicit
       conversion.
    explicit KNM(const KNM<R> & u)  // PB si stepi ou stepj nulle
          :KNM_<R>(new R[u.size()],u.N(),u.M())
         { KN_<R>::operator=(u);}
    explicit KNM(const KNM_<R> & u)
          :KNM_<R>(new R[u.size()],u.N(),u.M())
          { KNM_<R>::operator=(u);}
    */
    KNM( const KNM<R> & u )   // PB si stepi ou stepj nulle
            :
            KNM_<R>( new R[ u.size() ], u.N(), u.M() )
    {
        KN_<R>::operator=( u );
    }
    KNM( const KNM_<R> & u )
            : KNM_<R>( new R[ u.size() ], u.N(), u.M() )
    {
        KNM_<R>::operator=( u );
    }

    ~KNM()
    {
        delete [] this->v;
    }

    const KNM& operator=( const KNM_<const_R> & u )
    {
        KNM_<R>::operator=( u );
        return *this;
    }
    const KNM& operator=( const_R a )
    {
        KNM_<R>::operator=( a );
        return *this;
    }
    const KNM& operator+=( const_R a )
    {
        KNM_<R>::operator+=( a );
        return *this;
    }
    const KNM& operator-=( const_R a )
    {
        KNM_<R>::operator-=( a );
        return *this;
    }
    const KNM& operator/=( const_R a )
    {
        KNM_<R>::operator/=( a );
        return *this;
    }
    const KNM& operator*=( const_R a )
    {
        KNM_<R>::operator*=( a );
        return *this;
    }
    const KNM& operator+=( const KNM_<const_R> & u )
    {
        KNM_<R>::operator+=( u );
        return *this;
    }
    const KNM& operator-=( const KNM_<const_R> & u )
    {
        KNM_<R>::operator-=( u );
        return *this;
    }

    const KNM& operator/=( const KNM_<const_R> & u )
    {
        KNM_<R>::operator/=( u );
        return *this;
    }
    const KNM& operator*=( const KNM_<const_R> & u )
    {
        KNM_<R>::operator*=( u );
        return *this;
    }


    //  two opertors to cast to un array of constant
    //    operator KNM_<const_R> & ()
    //          { return *  (KNM_<const_R>*) this;}
    //    operator KNM_<const_R> const & ()  const
    //          { return *(const KNM_<const_R>*) this;}

    //    operator KNM<const_R> & ()
    //          { return *  (KNM<const_R>*) this;}
    //    operator KNM<const_R> const & ()  const
    //          { return *(const KNM<const_R>*) this;}

};

//  Array with 3 indices
//  ---------------------
template <class R>
class KNMK: public KNMK_<R>
{
public:

    KNMK( const int n, const int m, const int k )
            : KNMK_<R>( new R[ n*m*k ], n, m, k )
    {}
    explicit KNMK( const KNMK_<R> & u )
            : KNMK_<R>( new R[ u.size() ], u.N(), u.M(), u.K() )
    {
        KNMK_<R>::operator=( u );
    }
    explicit KNMK( const KNMK<R> & u )
            : KNMK_<R>( new R[ u.size() ], u.N(), u.M(), u.K() )
    {
        KNMK_<R>::operator=( u );
    }

    ~KNMK()
    {
        delete [] this->v;
    }

    KNMK& operator=( const KNMK_<const_R> & u )
    {
        KN_<R>::operator=( u );
        return *this;
    }
    KNMK& operator=( const_R a )
    {
        KN_<R>::operator=( a );
        return *this;
    }
    KNMK& operator+=( const_R a )
    {
        KN_<R>::operator+=( a );
        return *this;
    }
    KNMK& operator-=( const_R a )
    {
        KN_<R>::operator-=( a );
        return *this;
    }
    KNMK& operator/=( const_R a )
    {
        KN_<R>::operator/=( a );
        return *this;
    }
    KNMK& operator*=( const_R a )
    {
        KN_<R>::operator*=( a );
        return *this;
    }
    KNMK& operator+=( const KNMK_<const_R> & u )
    {
        KN_<R>::operator+=( u );
        return *this;
    }
    KNMK& operator-=( const KNMK_<const_R> & u )  // A. Veneziani: I guess here there was a bug
    {
        KN_<R>::operator-=( u );
        return *this;
    }

    KNMK& operator*=( const KNMK_<const_R> & u )   // A. Veneziani: I guess here there was a bug
    {
        KN_<R>::operator/=( u );
        return *this;
    }
    KNMK& operator/=( const KNMK_<const_R> & u )    // A. Veneziani: I guess here there was a bug
    {
        KN_<R>::operator/=( u );
        return *this;
    }

    //  two opertor to cast to un array of constant
    //    operator KNMK_<const_R> & ()
    //       { return *  (KNMK_<const_R>*) this;}
    //    operator KNMK_<const_R> const & ()  const
    //       { return *(const KNMK_<const_R>*) this;}

    //    operator KNMK<const_R> & ()
    //       { return *  (KNMK<const_R>*) this;}
    //    operator KNMK<const_R> const & ()  const
    //       { return *(const KNMK<const_R>*) this;}
};

//
// A. Veneziani, October, 30, 2002
//
//  Array with 4 indices
//  ---------------------
template <class R>
class KNMKL: public KNMKL_<R>
{
public:

    KNMKL( const int n, const int m, const int k, const int l )
            : KNMKL_<R>( new R[ n*m*k*l ], n, m, k, l )
    {}
    explicit KNMKL( const KNMKL_<R> & u )
            : KNMKL_<R>( new R[ u.size() ], u.N(), u.M(), u.K(), u.L() )
    {
        KNMKL_<R>::operator=( u );
    }
    explicit KNMKL( const KNMK<R> & u )
            : KNMKL_<R>( new R[ u.size() ], u.N(), u.M(), u.K(), u.L() )
    {
        KNMKL_<R>::operator=( u );
    }

    ~KNMKL()
    {
        delete [] this->v;
    }

    KNMKL& operator=( const KNMKL_<const_R> & u )
    {
        KN_<R>::operator=( u );
        return *this;
    }
    KNMKL& operator=( const_R a )
    {
        KN_<R>::operator=( a );
        return *this;
    }
    KNMKL& operator+=( const_R a )
    {
        KN_<R>::operator+=( a );
        return *this;
    }
    KNMKL& operator-=( const_R a )
    {
        KN_<R>::operator-=( a );
        return *this;
    }
    KNMKL& operator/=( const_R a )
    {
        KN_<R>::operator/=( a );
        return *this;
    }
    KNMKL& operator*=( const_R a )
    {
        KN_<R>::operator*=( a );
        return *this;
    }
    KNMKL& operator+=( const KNMKL_<const_R> & u )
    {
        KN_<R>::operator+=( u );
        return *this;
    }
    KNMKL& operator-=( const KNMKL_<const_R> & u )
    {
        KN_<R>::operator-=( u );
        return *this;
    }

    KNMKL& operator*=( const KNMKL_<const_R> & u )
    {
        KN_<R>::operator/=( u );
        return *this;
    }
    KNMKL& operator/=( const KNMKL_<const_R> & u )
    {
        KN_<R>::operator/=( u );
        return *this;
    }

    //  two opertor to cast to un array of constant
    //    operator KNMK_<const_R> & ()
    //       { return *  (KNMK_<const_R>*) this;}
    //    operator KNMK_<const_R> const & ()  const
    //       { return *(const KNMK_<const_R>*) this;}

    //    operator KNMK<const_R> & ()
    //       { return *  (KNMK<const_R>*) this;}
    //    operator KNMK<const_R> const & ()  const
    //       { return *(const KNMK<const_R>*) this;}
};


//  -------------  optimization ---------------------
template <class R>
class Add_KN_
{
public:
    const KN_<const_R> & a;
    const KN_<const_R> & b;
    Add_KN_( const KN_<const_R> & aa, const KN_<const_R> & bb )
            : a( aa ), b( bb )
    {
        K_assert( SameShape( a, b ) );
    }
};

template <class R>
class Sub_KN_
{
public:
    const KN_<const_R> & a;
    const KN_<const_R> & b;
    Sub_KN_( const KN_<const_R> & aa, const KN_<const_R> & bb )
            : a( aa ), b( bb )
    {
        K_assert( SameShape( a, b ) );
    }
};

template <class R>
class Mulc_KN_
{
public:
    const KN_<const_R> & a;
    const_R b;
    Mulc_KN_( const KN_<const_R> & aa, const_R bb ) : a( aa ), b( bb )
    {}
    Mulc_KN_( const Mulc_KN_<R> & aa, const_R bb ) : a( aa.a ), b( aa.b*bb )
    {}
}
;

template <class R>
class Add_Mulc_KN_
{
public:
    const KN_<const_R> a, b;
    const R ca, cb;
    Add_Mulc_KN_( const Mulc_KN_<R> & aa, const Mulc_KN_<R> & bb )
            : a( aa.a ), b( bb.a ), ca( aa.b ), cb( bb.b )
    {
        K_assert( SameShape( a, b ) );
    }
    Add_Mulc_KN_( const Mulc_KN_<R> & aa, const KN_<const_R> & bb, const R cbb )
            : a( aa.a ), b( bb ), ca( aa.b ), cb( cbb )
    {
        K_assert( SameShape( a, b ) );
    }
    Add_Mulc_KN_( const KN_<const_R> & aa, const R caa, const KN_<const_R> & bb, const R cbb )
            : a( aa ), b( bb ), ca( caa ), cb( cbb )
    {
        K_assert( SameShape( a, b ) );
    }
};


template <class R>
class Mul_KNM_KN_
{
public:
    const KNM_<const_R> A;
    const KN_<const_R> b;
    Mul_KNM_KN_( const KNM_<const_R> & aa, const KN_<const_R> & bb )
            : A( aa ), b( bb )
    {
        K_assert( SameShape( A.shapej, b ) );
    }
};


std::ostream & operator<<( std::ostream & f, const ShapeOfArray & s );

template <class R>
std::ostream & operator<<( std::ostream & f, const KN_<const_R> & v );
template <class R>
std::ostream & operator<<( std::ostream & f, const KNM_<const_R> & v );
template <class R>
std::ostream & operator<<( std::ostream & f, const KNMK_<const_R> & v );
template <class R>
inline std::ostream & operator<<( std::ostream & f, const KN<const_R> & v )
{
    return f << ( KN_<const_R> ) v;
}
template <class R>
inline std::ostream & operator<<( std::ostream & f, const KNM<const_R> & v )
{
    return f << ( KNM_<const_R> ) v;
}
template <class R>
inline std::ostream & operator<<( std::ostream & f, const KNMK<const_R> & v )
{
    return f << ( KNMK_<const_R> ) v;
}


template <class R>
inline Add_KN_<R> operator+( const KN_<const_R> &a, const KN_<const_R> &b )
{
    return Add_KN_<R>( a, b );
}
template <class R>
inline Sub_KN_<R> operator-( const KN_<const_R> &a, const KN_<const_R> &b )
{
    return Sub_KN_<R>( a, b );
}
template <class R>
inline Mulc_KN_<R> operator*( const KN_<const_R> &a, const R &b )
{
    return Mulc_KN_<R>( a, b );
}
template <class R>
inline Mulc_KN_<R> operator*( const R &b, const KN_<const_R> &a )
{
    return Mulc_KN_<R>( a, b );
}
template <class R>
inline Mulc_KN_<R> operator-( const KN_<const_R> &a )
{
    return Mulc_KN_<R>( a, -1 );
}



template <class R>
inline Add_Mulc_KN_<R> operator+( const Mulc_KN_<R>& a, const Mulc_KN_<R> &b )
{
    return Add_Mulc_KN_<R>( a, b );
}
template <class R>
inline Add_Mulc_KN_<R> operator-( const Mulc_KN_<R>& a, const Mulc_KN_<R> &b )
{
    return Add_Mulc_KN_<R>( a, b.a, -b.b );
}

template <class R>
inline Add_Mulc_KN_<R> operator+( const Mulc_KN_<R>& a, const KN_<const_R> &b )
{
    return Add_Mulc_KN_<R>( a, b, 1.0 );
}
template <class R>
inline Add_Mulc_KN_<R> operator-( const Mulc_KN_<R>& a, const KN_<const_R> &b )
{
    return Add_Mulc_KN_<R>( a, b, -1.0 );
}

template <class R>
inline Add_Mulc_KN_<R> operator+( const KN_<const_R> & b, const Mulc_KN_<R>& a )
{
    return Add_Mulc_KN_<R>( a, b, 1.0 );
}
template <class R>
inline Add_Mulc_KN_<R> operator-( const KN_<const_R> & b, const Mulc_KN_<R>& a )
{
    return Add_Mulc_KN_<R>( a, b, -1.0 );
}
template <class R>
inline Mul_KNM_KN_<R> operator*( const KNM_<const_R> A, const KN_<const_R> b )
{
    return Mul_KNM_KN_<R>( A, b );
}




template <class R>
inline int SameAdress( const KN_<R> &a, const KN_<R> &b )
{
    return & a[ 0 ] == &b[ 0 ];
}
}
#include <life/lifearray/RNM_tpl.hpp>
#ifdef K_assert
#undef K_assert
#endif
#endif

