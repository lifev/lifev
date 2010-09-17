/*
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file TimeAdvance_template.hpp
  \author F. Nobile,  M. Pozzoli,  C. Vergara

  File containing a class for an easy handling of different order time
  discretizations/extrapolations 

*/
#ifndef _TIMEADVANCE_TEMPLATE_H
#define _TIMEADVANCE_TEMPLATE_H
#include <string>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <life/lifecore/factory.hpp>
#include <life/lifearray/EpetraVector.hpp>

namespace LifeV
{
  
/*!
\class timeAdvance (templated)

This class define an astract method to build temporal discretitation schemes.
In particular we consider problems of the first order and the second  order in time;
that after space and time discretitation, and suitable linearization of non linear operator,
we obtain a linear system to solve at each time step (or each iteration ):

\f$ K U^{n+1} =F^{n+1}\f$
 
where K is an opportune matrix, U^{n+1} is the unknown vector and F^{n+1} is the
right hand side vector at the time tn+1 .
To determine F^{n+1} we deﬁne the state vector X^{n+1} that contained the infor-
mations about previous solutions.

1)   First order problems:

\f$ M \frac{u' }{\Delta t} + L( u) = f ,\f$

where L is a generic nonlinear operator.

We define U an approximation of u and the velocity vector V  an approximation of u'
at the time step $n+1$ as

\f$ V^{n+1}=\frac{\alpha_0}{\Delta t} U^{n+1}-  f_V^{n+1}  \f$,

where $\alpha_0$ is a suitable coefficient and  f_V  is a linear combination of the previous  solutions
 \Xbf^{n+1} with coefficients \alpha_i, that  will be specified in the following.\\

Then the time discrete version is

\f$ \frac{\alpha_0}{\Delta t}M U^{n+1}  +A (U^{n+1})= f^{n+1}+ M  f_V^{n+1}\f$,
that can be solved by any non-linear iterative solver (fixed point, Newton, ....).
The class \verb"timeAdvance" provides also a suitable extrapolation $U^*$ of  $\U^{n+1}$, given by a linear
combination of previous solutions and of order consistent with the time discretization formula.

2)  In this part we want to extend the previous approach to   second order
problems in time.

Second order problems:

\f$ u''  + L( u', u) = f \f$
wher L  is non linear operator.

We consider the following semidiscrete problem

\f$ M \frac{d^2U}{d t^2} + D ( U,  \frac{d U}{d t} ) + A ( U ) =\f  \f$,

where $M$ is the mass matrix, $\f$ the forcing term
 and $\U$ is the  vector of unknowns.
We define the following quantities

\f$ V :=\frac{d U}{d  t } &  W := \frac{d^2  U} {d  t^2}  \f$,
where V and W are the velocity and the  acceleration vectors, respectively.

At the time step $n+1$, we consider the following approssimations of V and W:

\f$ V^{n+1}=\frac{\alpha_0}{\Delta t} U^{n+1}- f_V^{n+1}  \f$

and 

\f$ W^{n+1}=\frac{\xi_0}{\Delta t^2} U^{n+1} - f_W^{n+1} \f$, 

where $f_V^{n+1}$ and $f_W^{n+1}$ are  linear combinations of the previous  solutions with
suitable coefficients alpha_i and  \xi_i, respectively. 
If  A and D depend on  U and V we can linearize the
problem using suitable extrapolations U^* V^* obtained   by linear combinations
 of previous solutions with coefficients  $\beta_i^U$ and $\beta_i^V$ respectively.

*/


template<typename VectorType = EpetraVector >

class TimeAdvance // T means template
{
public:

    /** @name Typedefs
     */
    //@{
    typedef VectorType                      vector_raw_type;
    typedef std::vector< vector_raw_type* > vector_type;
    typedef typename vector_type::iterator  vector_type_iterator;
    //@}

public:

    /*! Constructor
     *  @param n order of the BDF
     */
     TimeAdvance();

  //  TimeAdvance( const UInt order);
   
  // TimeAdvance( const UInt order,  const  UInt orderDev );
     
  // TimeAdvance( const double* coefficients,  const UInt orderDev );
    
    ~TimeAdvance();
  
    //! Initialize all the entries of the unknown vector to be derived with the
    //! vector u0 (duplicated)
    virtual void initialize_unk( VectorType u0 ) = 0;
  
    virtual void initialize_unk( VectorType u0, VectorType v0 ) = 0;
  
    virtual  void initialize_unk( VectorType u0, VectorType v0, VectorType const & w0 ) = 0;

  void initialize_rhs(const  VectorType & u0 ) ;


  //! Initialize all the entries of the unknown vector to be derived with a
  //! set of vectors uv0
  //! note: this is taken as a copy (not a reference), since uv0 is resized inside the method.
  virtual void initialize_unk(const  std::vector<VectorType> uv0 ) = 0;

  //!initialize parameters of time advance scheme
  void setup ( const  UInt orderDev ) { _M_orderDev = orderDev;} 

  virtual void setup ( const UInt order,  const  UInt orderDev ) = 0;

  virtual void setup ( const   std::vector<double>  coefficients, const  UInt orderDev ) = 0;

  /*! Update the vectors of the previous time steps by shifting on the right
     *  the old values.
     *  @param u_curr current (new) value of the state vector
     */
  virtual void shift_right(const VectorType & u_curr ) =0;
  
  //!spy unknown vector
   void spy() ;
  
   //!spy rhs vector
   void spy_rhs();

   //! initializa time step
  void setDeltaT( Real dt )  { _M_dt = dt; };

  //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
  //! formula
 virtual  VectorType time_der( Real dt = 1 ) /* const */= 0;
  
  //! Returns the right hand side \f$ \bar{p} \f$ of the time derivative
  //! formula
  virtual VectorType time_derOrder2( Real dt = 1 ) /*const*/ = 0;
  
  //! Compute the polynomial extrapolation approximation of order n-1 of
  //! u^{n+1} defined by the n stored state vectors
  virtual VectorType extrap() const  = 0;
  
  //! Compute the polynomial extrapolation approximation of order k-1 in iterative 
  //! methods  u^{k+1} defined by the unk and the n  stored state vectors  
  virtual VectorType extrap(const  VectorType unk ) = 0;

  //! Compute the polynomial extrapolation approximation of order n-1 of
  //! u^{n+1} defined by the n stored state vectors
   virtual VectorType extrapVelocity()  const = 0;
  
  //! Return the i-th coefficient of the time derivative alpha_i
   double coeff_der( UInt i )  const;
  
  //! Return the i-th coefficient of the second time derivative xi_i
  double coeff_derOrder2( UInt i ) const;
  
  //! Return the i-th coefficient of the time extrapolation beta_i
  virtual double coeff_ext( UInt i )  const = 0;
  
  //! Return the i-th coefficient of the time extrapolation betaV_i
  virtual double coeff_extVelocity( UInt i ) const =0;
  
  //! Return the last unknown vector
  const  VectorType unk()  const;

 //! Return the last unknown vector
  const  VectorType unk(const UInt i)  const;

  //! Return the velocity 
   virtual  VectorType vnk() const = 0;

 //! Return the accelerate
   virtual VectorType wnk() const = 0;

  // !Return velocity's right hand side 
  const VectorType rhsV() ;

  //! Return accelerate's right hand side
  const VectorType rhsW() ;

  //! Return the velocity associated to u
 // used for example in FSI to return the value of solid
  //  in the internal loop
  const VectorType  vnk(const  VectorType & u);
  
  //! Return the accelerate associated to u;
  // used for example in FSI to return the value of solid
  //  in the internal loop
  const VectorType  wnk(const  VectorType & u);

  //show the proprities of temporal scheme
  virtual void showMe()  const = 0;
  
  UInt order() const  {return _M_order;}

  //protected:
public :
    //! Order of the BDF derivative/extrapolation: the time-derivative
    //! coefficients vector has size n+1, the extrapolation vector has size n
    UInt _M_order;
   
    //! Order of temporal derivate: the time-derivative
    //! coefficients vector has size n+1, the extrapolation vector has size n
    UInt _M_orderDev;
  
    //! time step
    Real _M_dt;

   //! Size of the unknown vector
   UInt _M_size;

    //! Size for time_der loop (for bdf  equal M_order, for Newmark equal  _M_size/2)
   UInt _M_sizeTimeDer;

   //! Size for time_derOrder2 loop  (for bdf  equal M_order, for Newmark equal _M_size/2)
   UInt _M_sizeTimeDer2;
  
  //!Size of coefficients (for bdf equal M_order+_M_orderDev, for theta-method is 3, and Newmark is 4)
  UInt _M_sizeCoefficients;

  //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
  Vector _M_xi;
  
  //! Coefficients \f$ \alpha_i \f$ of the time bdf discretization
  Vector _M_alpha;
  
  //! Coefficients \f$ \beta_i \f$ of the extrapolation
  Vector _M_beta;
  
 //! Coefficients \f$ \beta_i \f$ of the extrapolation
  Vector _M_beta2;

  //! Last n state vectors
  vector_type _M_unknowns;
  
  //! Vector of rhs (rhsV and rhsW)
  vector_type _M_rhs;
  
};
 

 // ===================================================
//! MACROS
// ===================================================
 
typedef singleton< factory < TimeAdvance<>,  std::string> > TimeAdvanceFactory;
  
template<typename VectorType>
TimeAdvance<VectorType>::TimeAdvance( )
    :
  _M_unknowns(),
  _M_rhs()
{
  _M_unknowns.reserve( 1 );
  _M_rhs.reserve(2);
}
  
template<typename VectorType>
TimeAdvance<VectorType>::~TimeAdvance()
  {
    vector_type_iterator iter     = _M_unknowns.begin();
    vector_type_iterator iter_end = _M_unknowns.end();
    
    for ( ; iter != iter_end; iter++ )
        delete *iter;
  }

template<typename VectorType>
double
TimeAdvance<VectorType>::coeff_der( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i < _M_sizeCoefficients,
            "Error in specification of the time derivative coefficient for the time scheme" );
    return _M_alpha[ i ];
}

template<typename VectorType>
double
TimeAdvance<VectorType>::coeff_derOrder2( UInt i )  const
{
    // Pay attention: i is c-based indexed
  ASSERT( i < _M_sizeCoefficients,
            "Error in specification of the time derivative coefficient for the time scheme" );
    return _M_xi[ i ];
}

template<typename VectorType>
const VectorType 
TimeAdvance<VectorType>::unk() const
{
  VectorType u(*_M_unknowns[0]);
  return u;
}

template<typename VectorType>
const VectorType 
TimeAdvance<VectorType>::unk(const UInt i) const
{
// Pay attention: i is c-based indexed
    ASSERT( i < _M_size,
            "Error there isn't unk(i), i must be shorter than M_size" );

  VectorType u(*_M_unknowns[i]);
  return u;
}

//! Return the velocity associated to u
template<typename VectorType>
const VectorType 
TimeAdvance<VectorType>:: vnk(const  VectorType & u)
{
  VectorType v( u ); 
  v  *= _M_alpha[ 0 ] / _M_dt;
  v   -= (*this->_M_rhs[ 0 ]);
  return v;
}

//! Return the accelerate associated to u
template<typename VectorType>
const VectorType 
TimeAdvance<VectorType>:: wnk(const VectorType & u)
{
 VectorType w(u); 
 w  *= _M_xi[ 0 ] / (_M_dt*_M_dt);
 w  -= (*this->_M_rhs[1]);
  return w;
}

template<typename VectorType>
void
TimeAdvance<VectorType>:: initialize_rhs(const  VectorType & u0 )
{
    for (UInt i=0; i<2; i++ ) 
      { 
	_M_rhs.push_back(new VectorType(u0.getMap(), Unique));
	*_M_rhs[i] *=0;
      }
}


template<typename VectorType>
const VectorType
TimeAdvance<VectorType>::rhsV() { return *_M_rhs[0];}

template<typename VectorType>
const VectorType
TimeAdvance<VectorType>::rhsW() { return *_M_rhs[1];}



template<typename VectorType>
void
TimeAdvance<VectorType>::
spy()
{
  static UInt saveUnknowns=0;
  std::string unknowns="unknowns";
  
  for( UInt i=0 ; i< _M_size ; i++ )
    {
      std::ostringstream ii;
      ii<<saveUnknowns;     ii<<i;
      
      unknowns+ii.str();
     
      _M_unknowns[i]->spy(unknowns+ii.str());
    }
  saveUnknowns++;
}


template<typename VectorType>
void
TimeAdvance<VectorType>::
spy_rhs()
{
  static UInt saveRhs=0;
  std::string rhs="rhs";
  for( UInt i=0 ; i< 2 ; i++ )
    {
      std::ostringstream ii;
      ii<<saveRhs;     ii<<i;
      
      rhs+ii.str();
      _M_rhs[i]->spy(rhs+ii.str());
    }
  saveRhs++;
}

}
#endif
