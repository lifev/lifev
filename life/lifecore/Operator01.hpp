//=========================================================
// Operator PurePhi Diego Mastalli.
//=========================================================

template< typename FE >
class PurePhi
{ 
public: 

PurePhi< FE >( FE* fe ): _fe( fe ){} ;
 
FE* fe_ptr() { return _fe ; } ;
 
Real operator()( UInt k  , UInt iq )
{
	return _fe->phi( k , iq );
}

private:

FE* _fe ; 
};

//=========================================================
// Operator PureGrad Diego Mastalli.
//=========================================================

template< typename FE >
class PureGrad
{
public: 

PureGrad< FE >( FE* fe ): _fe( fe ){} ;
 
FE* fe_ptr() { return _fe ; } ;
 
Real operator()( UInt k , UInt icoor , UInt iq )
{
	return _fe->phiDer( k , icoor , iq );
}

private:

FE* _fe ; 
};

//=========================================================
// Operator PureDer2 Diego Mastalli.
//=========================================================

template< typename FE >
class PureDer2
{ 
public: 

PureDer2< FE >( FE* fe ): _fe( fe ){} ;
 
FE* fe_ptr() { return _fe ; } ;
 
Real operator()( UInt k , UInt icoor , UInt jcoor , UInt iq )
{
	return _fe->phiDer2( k , icoor , jcoor , iq );
}

private:

FE* _fe ; 
};

//=========================================================
// Operator Laplacian Diego Mastalli.
//=========================================================

template< typename PHD >
class Laplacian
{
public: 

Laplacian<  PHD > (PHD& phD): _phD(phD){};
 
Real operator()( UInt k , UInt iq )
{
	Real s = 0 ;
    
    for( UInt icoor = 0 ; icoor < _phD.fe_ptr()->nbCoor ; icoor++ )
    {
   		s += _phD( k , icoor , icoor , iq ) ;
    }

    return s ;
}

private:

PHD _phD;
};

//=========================================================
// Operator LS Diego Mastalli.
//=========================================================

template< typename L , typename P >
class LS
{
public: 

LS<L,P>(L& Lap,P& pBas,Real sigma,Real mu):_Lap(Lap),_pBas(pBas),_sigma(sigma),_mu(mu){};
 
Real operator()( UInt k  , UInt iq )
{
    return _mu * _Lap( k , iq ) + _sigma * _pBas( k , iq ) ;
}

private:

L _Lap ; 
P _pBas ;
Real _sigma ;
Real _mu ;

};

//=========================================================
// Operator LSS Diego Mastalli.
//=========================================================

template< typename PHD , typename VT >
class LSS
{
public: 

LSS<PHD,VT>(PHD& phD,VT& Beta): _phD(phD),_Beta(Beta){};

FE* fe_ptr() { return _fe ; } ;

Real operator()( UInt k , UInt iq )
{
	Real s = 0 ;
	
	for ( UInt icoor = 0 ; icoor < _phD.fe_ptr()->nbCoor ; icoor++ )
    {
    	s += _Beta[ icoor ] * _phD( k , icoor , iq ) ;
    }
    
    return s ;
}

private:
 
PHD _phD ;
VT _Beta ;
};

//=========================================================
// Operator Mass Diego Mastalli.
//=========================================================

template< typename PHB >
class Mass
{
public: 

Mass< PHB >( PHB& phB ): _phB(phB) {};
 
Real operator()( UInt i , UInt j , UInt iq )
{
	return _phB( i , iq ) * _phB( j , iq ) ;
}

private:

PHB _phB ;
};

//=========================================================
// Operator Stiff Diego Mastalli.
//=========================================================

template< typename PHD >
class Stiff
{ 
public:
  
Stiff< PHD >( PHD& phD ): _phD( phD ){} ;
 
Real operator()( UInt i, UInt j , UInt iq )
{
	Real s = 0 ;
	
	for ( UInt icoor = 0; icoor < _phD.fe_ptr()->nbCoor ; icoor++ )
	{
		s += _phD( i , icoor , iq ) * _phD( j , icoor , iq ) ;
	}
	
	return s ;
}

private:
 
PHD _phD ; 
};

//=========================================================
// Operator Grad Diego Mastalli.
//=========================================================

template< UInt coor, typename PHB , typename PHD >
class Grad
{
public: 

Grad< coor , PHB , PHD >( PHB& phB , PHD& phD ): _phB( phB ) , _phD( phD ) {} ;
 
Real operator()( UInt i , UInt j , UInt iq )
{
 	return _phB( i , iq ) * _phD( j , coor , iq ) ;
}

private:

PHB _phB ;
PHD _phD ;
};