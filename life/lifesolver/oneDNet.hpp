/* -*- mode: c++ -*-

   This file is part of the LifeV library

   Author(s): Tiziano Passerini <passerini@mate.polimi.it>
   Date: 2005-06-17

   Copyright (C) 2005 Politecnico di Milano

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

/*!
  \file oneDNet.hpp
  \author Tiziano Passerini <passerini@mate.polimi.it>
  \date 12/2006
  \version 1.5
*/

#ifndef _ONEDNET_H_
#define _ONEDNET_H_

#include <map>
#include <list>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <life/lifecore/debug.hpp>

#include <life/lifealg/clapack.hpp>

#include <life/lifearray/tab.hpp>

#include <life/lifefem/oneDBCFunctions.hpp>

#include <life/lifearray/EpetraVector.hpp>
#include <life/lifealg/EpetraMap.hpp>

namespace LifeV
{

/*!
  \brief A class to deal with 1D model networks.

  This class describes a network of tubes, each one
  represented by a 1 D model.
  The underlying data structure is a graph, so that
  one can regard the network as a set of vertices and edges
*/
template< class SOLVER1D, class PARAM1D >
class OneDNet{

public:

    struct OneDVesselsInterface;
    struct OneDVessel;

    /*! \name Typedefs
     */
    //@{
    //! Boost shared pointer to 1D parameter class
    typedef typename boost::shared_ptr< PARAM1D > OneDParamPtr;
    //! Map of pointers to 1D parameter class

    //! Boost shared pointer to 1D solver class
    typedef typename boost::shared_ptr< SOLVER1D > OneDSolverPtr;

    //! Network: based on boost::graph library
    /*!
      \param boost::listS Type for edge list (STL list)
      \param boost::vecS Type for vertex list (STL vector)
      \param boost::bidirectionalS Sets an oriented graph
      \param OneDVesselsInterface Attach vertex properties
      \param OneDVessel Attach edge properties
    */
    typedef boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        OneDVesselsInterface,
        OneDVessel
        >
    Network;

    //! Iterator type for the list of vertices
    typedef typename boost::graph_traits< Network >::vertex_iterator
    Vertex_Iter;
    //! Descriptor type for the list of vertices
    typedef typename boost::graph_traits< Network >::vertex_descriptor
    Vertex_Descr;

    //! Iterator type for the list of edges
    typedef typename boost::graph_traits< Network >::edge_iterator
    Edge_Iter;
    //! Descriptor type for the list of edges
    typedef typename boost::graph_traits< Network >::edge_descriptor
    Edge_Descr;

    //! Iterator type for the list or out edges (exiting from a vertex)
    typedef typename boost::graph_traits< Network >::out_edge_iterator
    Out_Edge_Iter;
    //! Iterator type for the list or in edges (exiting from a vertex)
    typedef typename boost::graph_traits< Network >::in_edge_iterator
    In_Edge_Iter;

    //! Type for degree size (number of edges connected to a vertex)
    typedef typename boost::graph_traits< Network >::degree_size_type
    Degree_Size;

    //! Map of pointers to 1D solver class
    typedef typename std::map< int, OneDSolverPtr > MapSolver;
    //! Value type for solver map
    typedef typename MapSolver::value_type MapSolverValueType;

    //! Pair type: 2-dimensional vector
    typedef std::pair< Real, Real > Vec2D;

    //@}

    /*! \name Public Members
     */
    //@{
    //! Vertex properties
    /*!
      Boost::graph classes allow to specify properties
      for edges and vertices through the definition of
      classes (bundled properties).

      OneDVesselsInterface contains the properties
      associated to the vertices of the graph
    */
    struct OneDVesselsInterface{
        int index; /*!< numerical label */
        bool internal;  /*!< is it an internal interface? */
        int type; /*!< the type of the interface */
//        boost::shared_ptr<LifeV::EpetraVector> values;
    };

    //! Edge properties
    /*!
      \sa OneDVesselsInterface

      OneDVessel contains the properties
      associated to the edges of the graph
    */
    struct OneDVessel{
        int index; /*!< numerical label */
        OneDParamPtr onedparam; /*!< pointer to 1D parameter class */
        OneDSolverPtr onedsolver; /*!< pointer to 1D solver class */
//        Vec2D;
    };
    //@}

    /*! \name Constructors
     */
    //@{
    //! Constructor taking the data file
    OneDNet( GetPot data_file );
    //@}


    /*! \name Getters
     */
    //@{
    //! Gets a const reference to the data structure
    Network const& network();


    //! Returns the "out" degree of one vertex of the net
    /*!
      \param ind vertex label (starting from 1)
      \return number of outgoing edges
    */
    Degree_Size outDegreeSize( const int& ind );

    //! Returns the degree of one vertex of the net
    /*!
      \param ind vertex label (starting from 1)
      \return number of outgoing + incoming edges
    */
    Degree_Size degreeSize( const int& ind );

    //! Returns the list of incoming edges
    std::vector<In_Edge_Iter> inEdges( const int& ind );
    //! Returns the list of outgoing edges
    std::vector<Out_Edge_Iter> outEdges( const int& ind );

    //! Gets an iterator to vertex with specified index
    Vertex_Iter Vertex( const int& ind );

    //! Gets an iterator to edge with specified index
    Edge_Iter Edge( const int& ind );

    //! Gets a reference to 1D Solver (associated to the edge) with specified index
    SOLVER1D & Solver( const int& ind );

    //! Gets a pointer to 1D Solver (associated to the edge) with specified index
    OneDSolverPtr SolverPtr( const int& ind );

    //! Gets a reference to 1D Param (associated to the edge) with specified index
    PARAM1D & Param( const int& ind );

    //! Gets a reference to 1D Param (associated to the edge) with specified index
    OneDParamPtr ParamPtr( const int& ind );

    //! Return the time step
    Real timestep() const;
    //! Return the initial time value
    Real inittime() const;
    //! Return the end time value
    Real endtime()  const;

    //@

    /*! \name Methods
     */
    //@{
    //! Initializes the solvers in the net
    /*!
      \param data_file data file
    */
    void initialize( GetPot data_file );

    //! Adds function for boundary condition evaluation
    /*!
      \param fun function for boundary condition evaluation
	  \param ind index of 1D solver
	  \param border which boundary
	  \param line first or second line in the associated 2x2 linear system
	  \param var W1, W2, A, Q
    */
    void setBC( OneDBCFunctionPointer fun, const int& ind,
                std::string const& border, std::string const& line,
                std::string const& var );
    //! Computes the unknown at the interface between different solvers
    /*!
      Imposes interface conditions according to the \b type of each vertex

      Default: continuity of total pressure and mass conservation
    */
    void computeInterfaceTubesValues( );

    void updateInterfaces( );

    //! Create a copy of the current solution
    /*!
      For use in the (implicit) coupling algorithm (eg. 3D - 1D)

      At the beginning of each k+1 fixed point sub-iterations,
      one needs to restart the solvers from the solution at time t_n.
      savesol() and loadsol() allow you to do that

      \sa OneDModelSover
    */
    void savesol(  );

    //! Load the solution previously saved with savesol()
    /*!
      \sa savesol()
    */
    void loadsol(  );

    //! Visits the graph and updates the right hand side of each solver for time advancing
    /*!
      \param time current time
    */
    void timeAdvance( const Real& time );

    //! Visits the graph and invokes iterate() on each solver
    /*!
      \param time current time
      \param count iteration number (requested by LifeV::OneDModelSolver)
    */
    void iterate( const Real& time , const int& count);

    //! Prints out results
    /*!
      \param time current time
    */
    void postProcess( const Real& time );

    //! Open some std::ostringstream variables as file buffers
    /*!
      When the network is huge you expect to have huge amounts of data
      to be written down during post processing.

      Instead of writing down to file, at each call of postprocessing()
      you add data to the buffers, which will be redirected to file
      by output2FileBuffers

      \sa output2FileBuffers
    */
    void openFileBuffers();

    //! Redirect buffer variables to files
    /*!
      \param time_val current time
      \sa openFileBuffers()
    */
    void output2FileBuffers( const Real& time_val );

    //! Empty buffer variables
    /*!
      \sa openFileBuffers()
    */
    void resetFileBuffers();

    //! Sets put pointer in the output files
    /*!
      \sa openFileBuffers()
    */
    void seekpFileBuffers();

    //! Store get pointer for the output files
    /*!
      \sa openFileBuffers()
    */
    void tellpFileBuffers();

    //! Erase buffer variables
    /*!
      \sa openFileBuffers()
    */
    void closeFileBuffers();
    //@

    //boost::shared_ptr<Epetra_Map> _M_epetra_interfaces;
    //boost::shared_ptr<EpetraMap> _M_epetra_vessels;
    //boost::shared_ptr<EpetraMap > _M_epetra_map_interface_values;
    //boost::shared_ptr<EpetraVector<double> > _M_epetra_vector_interface_values;
    
    Epetra_Map interface_map() { return *_M_epetra_interfaces; } 
    Epetra_Map tube_map() { return *_M_epetra_vessels->getRepeatedEpetra_Map(); } 


private:

    //! Method for imposing continuity interface conditions
    /*!
      The interface conditions for branching can be written in the form:
      \f[
      \mathbf{f}( \mathbf{x} ) = 0
      \f]
      where (being \f$n\f$ the number of branches)
      \f[
      \left\{
      \begin{array}{rclc}
      f_0( \mathbf{x} ) & = & \sum_{k=0}^{n} signum(k) * Q_k \\
      f_i( \mathbf{x} ) & = & P_{t,i} - P_{t,0}, & for\ i = 1 \cdots n-1 \\
      f_{n+i}( \mathbf{x} ) & = & W_{out,i} - W_{out,i}^{ext}, & for\ i = 0 \cdots n-1
      \end{array}
      \right.
      \f]

      Note that each branch has a signum: positive if it is an in-edge,
      negative for out-edges.

      The solution can be found with a Newton method of the form
      \f[
      x^{k+1} = x^{k} - \frac{\mathbf{f}(\mathbf{x})} {\mathbf{f}^{'}(\mathbf{x})}
      \f]

      The unknown vector x (of size 2*n) is defined as follows:
      \f[
      x_{2i} = A_i; \quad x_{2i+1} = Q_i,\quad for\ i = 0 \cdots n-1
      \f]
    */
    void interface_continuity_conditions( Vertex_Iter const& vertex );

    //! Method computing f and its gradient
    /*!
      \sa interface_continuity_conditions( )

      \param x unknown vector of size 2*n
      \param f vector function of size 2*n
      \param jac matrix function of size (2*n)^2 (the jacobian of \f$f\f$)
      \param intTubes a container for in-edges and out-edges of considered vertex; note that
      \f[
      n = intTubes.size()
      \f]
      \param sign associate each edge with its signum:
      - true if in-edge
      - false if out-edge
      This affects the way each tube is contibuting to the equations
    */
    void f_jac(  const Vector& x,
                 Vector& f,
                 Matrix& jac,
                 MapSolver& intTubes,
                 std::map< int, bool >& sign
                 );

    //! Scalar product for 2D vectors
    Real dot(const Vec2D& vec1, const Vec2D& vec2) const;

    //! Find the solution at the foot of the characteristic line
    /*!
      \param point_bound spatial coordinate for boundary node (first interpolation point)
      \param point_internal spatial coordinate for internal node (second interpolation point)
      \param deltaT time step
      \param eigenvalue the eigenvalue associated to the considered characteristic variable
      \param U_bound solution at boundary node (first interpolation value)
      \param U_intern solution at internal node (second interpolation value)
      \return Solution \f$U\f$ at the foot of the characteristic line
    */
    Vec2D interpolLinear(const Real& point_bound, const Real& point_internal,
                         const Real& deltaT, const Real& eigenvalue,
                         const Vec2D& U_bound, const Vec2D& U_intern) const;

    const UInt _M_num_vessels; //!< number of edges in the graph
    const UInt _M_num_interfaces; //!<  number of vertices in the graph

    //! boost::adjacency_list variable
    Network _M_network;

    //! time parameters
    Real _M_time_beg;  //!< starting time
    Real _M_time_end; //!< finishing time
    Real _M_time_step; //!< time step

    /*!
      Assumption: all 1D solvers have the same time parameters
      ( starting time, finishing time, time step).
      You can make sure of that by imposing to each edge the network's
      time parameters:
      - true -> the network sets time parameters for each solver
      - false -> the network does not set time parameters
    */
    bool _M_set_time_param;

    boost::shared_ptr<Epetra_Comm>   _M_comm;

    boost::shared_ptr<Epetra_Map> _M_epetra_interfaces;
    boost::shared_ptr<EpetraMap> _M_epetra_vessels;
    boost::shared_ptr<EpetraMap > _M_epetra_map_interface_values;
    boost::shared_ptr<EpetraVector<double> > _M_epetra_vector_interface_values;

};



/*********** IMPLEMENTATION **************/


template< class SOLVER1D, class PARAM1D >
OneDNet< SOLVER1D, PARAM1D >::OneDNet( GetPot data_file ):
    _M_num_vessels( data_file("1dnetwork/num_vessels", 1) ),
    _M_num_interfaces( data_file("1dnetwork/num_interfaces", 2) ),
    _M_network( _M_num_interfaces ),
    _M_time_beg( data_file("1dnetwork/timebeg", 0.) ),
    _M_time_end( data_file("1dnetwork/timeend", 1.) ),
    _M_time_step( data_file("1dnetwork/timestep", 1e-3) ),
    _M_set_time_param( data_file("1dnetwork/set_time_param", 1) )
{
#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    _M_comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int ntasks;
    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    _M_comm.reset( new Epetra_SerialComm() );
#endif

    if (!_M_comm->MyPID())
        std::cout << "My PID = " << _M_comm->MyPID() << " out of " << ntasks << " running." << std::endl;

    // build epetra map for graph edges
//    _M_epetra_vessels.reset( new Epetra_Map( (int)_M_num_vessels, 0, *_M_comm ) );
    // build epetra map for graph vertices
    _M_epetra_interfaces.reset( new Epetra_Map( (int)_M_num_interfaces, 1, *_M_comm) );

    // debug: Epetra
    Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] interfaces (indices):\n";
    for( int i = 0; i < _M_epetra_interfaces->NumMyElements(); ++i )
        Debug( 6030 ) << *(_M_epetra_interfaces->MyGlobalElements() + i) << "\n\t";
    Debug( 6030 ) << "\n";

    _M_comm->Barrier();
//    char ready;
//    do std::cin >> ready; while(!ready);

    // store vertex descriptors (you will need them)
    Vertex_Descr vd[ _M_num_interfaces ];
    //    Vertex_Descr vd[ _M_epetra_interfaces->NumMyElements ];

    // vertex iterators (for visiting the graph)
    std::pair< Vertex_Iter, Vertex_Iter > vertex_iter_pair;
    // edge descriptors (for visiting the graph)
    std::pair< Edge_Descr, bool > edge_descr_pair;

    // auxiliary variable for vertex indexes
    UInt i(0); // note that indexes start from 1

    // auxiliary variable (will contain interface label)
    std::string str;
    // point to the correct section in data file
    data_file.set_prefix( "1dnetwork/interfacemap/" );

    // assign properties to vertices
    for( vertex_iter_pair = vertices(_M_network);
         vertex_iter_pair.first != vertex_iter_pair.second;
         ++vertex_iter_pair.first, ++i ){
        //    for( int i = 0; i < _M_epetra_vessels->NumMyElements(); ++i )

        vd[i] = *vertex_iter_pair.first;
        // set index for the vertex
        _M_network[ vd[i] ].index = i+1;

        // Conversion using lexical_cast
        try {
            // store interface index in string variable
            str=boost::lexical_cast<std::string>(i+1);
        }
        catch (boost::bad_lexical_cast &) {
            // conversion failed
            std::cout << "\n[OneDNet::OneDNet]: lexical conversion failed!"
                      << std::endl;
        }
        // set type for the vertex
        _M_network[ vd[i] ].type = data_file(str.c_str(),1);
    }

    // labels for vertices at edges' left and right boundary
    UInt left, right;

    // access to data file
    /*
      GetPot files can be accessed at different sections;
      prefix contains the section name and is modified at each iteration
      (last three letters must be substituted by str)
    */
    std::string prefix("1dnetwork/tubenn/");
    std::string::iterator prefix_end = prefix.end();

    std::cout << "\n*** Class OneDNet ***\n" << std::flush;
    // debug: did I read the data file properly?
    Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] The net contains "
    							<< _M_num_vessels << " vessels and "
                  << _M_num_interfaces << " interfaces.\n";

    for( i=0; i < _M_num_vessels; ++i ){

        // Conversion using lexical_cast
        try {
            // store edge index in string variable
            str=boost::lexical_cast<std::string>(i+1);
        }
        catch (boost::bad_lexical_cast &) {
            // conversion failed
            std::cout << "\n[OneDNet::OneDNet proc " << _M_comm->MyPID() << "] lexical conversion failed!"
                      << std::endl;
        }

        // update prefix
        prefix.replace( prefix_end - 3, prefix_end - 1, str );
        // point to the correct section in data file
        data_file.set_prefix( prefix.c_str() );

        left = data_file("left_interface",1);
        right = data_file("right_interface",1);

        // debug: did I read the data file properly?
        Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] Tube number " << str << "\n";
        Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "]\t\tleft interface = " << left << "\n";
        Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "]\t\tright interface = " << right << "\n";

        // add (oriented) edge between connected vertices
        edge_descr_pair = add_edge( vd[left-1], vd[right-1], _M_network);

        // check if add_edge failed
        std::stringstream assert_sstr;
        assert_sstr << "\n[OneDNet::OneDNet process " << _M_comm->MyPID() << "] ERROR adding edge!"; 
        LIFEV_ASSERT( edge_descr_pair.second == true ).error( assert_sstr.str().c_str() );

        // set label for the edge
        _M_network[ edge_descr_pair.first ].index = i+1;

        // reset prefix to initial value
        prefix = "1dnetwork/tubenn/";
        prefix_end = prefix.end();

    }

    _M_comm->Barrier();
    // build epetra map for vessels
    std::set<int> my_vessels; 
    // build epetra map for interface values
    std::set<int> my_interface_values; 
    // pair of in-edge iterators
    std::pair< In_Edge_Iter, In_Edge_Iter > in_edge_iter_pair;
    // pair of out-edge iterators
    std::pair< Out_Edge_Iter, Out_Edge_Iter > out_edge_iter_pair;
    for( int i = 0; i < _M_epetra_interfaces->NumMyElements(); ++i ) {
        Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] consider vertex "
        							<< *(_M_epetra_interfaces->MyGlobalElements() + i);
        Debug( 6030 ) << "\n\t";

    	for( in_edge_iter_pair = in_edges(vd[*(_M_epetra_interfaces->MyGlobalElements() + i) - 1], _M_network);
         in_edge_iter_pair.first != in_edge_iter_pair.second;
         ++in_edge_iter_pair.first )
        {
    				Debug( 6030 ) << "\n\t[OneDNet::OneDNet process " << _M_comm->MyPID() << "] adding vessel "
    											<< _M_network[ *(in_edge_iter_pair).first ].index
    											<< " to my map";
//            Debug( 6030 ) << "\n\t";
    				my_vessels.insert( _M_network[ *(in_edge_iter_pair).first ].index );
    				Debug( 6030 ) << "\n\t[OneDNet::OneDNet process " << _M_comm->MyPID() << "] adding interface vector indexes "
    											<< 4 * (_M_network[ *(in_edge_iter_pair).first ].index - 1) + 1 << "\t"
    											<< 4 * (_M_network[ *(in_edge_iter_pair).first ].index - 1) + 2 << "\t"
    											<< 4 * (_M_network[ *(in_edge_iter_pair).first ].index - 1) + 3 << "\t"
    											<< 4 * (_M_network[ *(in_edge_iter_pair).first ].index - 1) + 4 << "\t"
    											<< " to my map";
//    				Debug( 6030 ) << "\n\t";
    				my_interface_values.insert( 4 * (_M_network[ *(in_edge_iter_pair).first ].index - 1) + 1 );
    				my_interface_values.insert( 4 * (_M_network[ *(in_edge_iter_pair).first ].index - 1) + 2 );
    				my_interface_values.insert( 4 * (_M_network[ *(in_edge_iter_pair).first ].index - 1) + 3 );
    				my_interface_values.insert( 4 * (_M_network[ *(in_edge_iter_pair).first ].index - 1) + 4 );
        }
    for( out_edge_iter_pair = out_edges(vd[*(_M_epetra_interfaces->MyGlobalElements() + i) - 1], _M_network);
         out_edge_iter_pair.first != out_edge_iter_pair.second;
         ++out_edge_iter_pair.first )
        {
    				Debug( 6030 ) << "\n\t[OneDNet::OneDNet process " << _M_comm->MyPID() << "] adding vessel "
    											<< _M_network[ *(out_edge_iter_pair).first ].index
    											<< " to my map";
//            Debug( 6030 ) << "\n\t";
						my_vessels.insert( _M_network[ *(out_edge_iter_pair).first ].index );
    				Debug( 6030 ) << "\n\t[OneDNet::OneDNet process " << _M_comm->MyPID() << "] adding interface vector indexes "
    											<< 4 * (_M_network[ *(out_edge_iter_pair).first ].index - 1) + 1 << "\t"
    											<< 4 * (_M_network[ *(out_edge_iter_pair).first ].index - 1) + 2 << "\t"
    											<< 4 * (_M_network[ *(out_edge_iter_pair).first ].index - 1) + 3 << "\t"
    											<< 4 * (_M_network[ *(out_edge_iter_pair).first ].index - 1) + 4 << "\t"
    											<< " to my map";
//            Debug( 6030 ) << "\n\t";
    				my_interface_values.insert( 4 * (_M_network[ *(out_edge_iter_pair).first ].index - 1) + 1 );
    				my_interface_values.insert( 4 * (_M_network[ *(out_edge_iter_pair).first ].index - 1) + 2 );
    				my_interface_values.insert( 4 * (_M_network[ *(out_edge_iter_pair).first ].index - 1) + 3 );
    				my_interface_values.insert( 4 * (_M_network[ *(out_edge_iter_pair).first ].index - 1) + 4 );
        }
    }
    int my_vessels_array[my_vessels.size()], num(0);
		Debug( 6030 ) << "\n[OneDNet::OneDNet process " << _M_comm->MyPID()
									<< "] this list should be non repeated (vessels index set, size "
									<< my_vessels.size() << ")";
    for( std::set<int>::iterator it = my_vessels.begin();
    	it != my_vessels.end(); ++it, ++num )
    {
    	Debug( 6030 ) << *it << " \t";
    	my_vessels_array[num] = *it;
    }
    _M_epetra_vessels.reset( new EpetraMap( -1,
                                             my_vessels.size(),
                                             my_vessels_array,
                                             1, *_M_comm ) );

    Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] vessels indices (repeated map):\n";
    for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i )
        Debug( 6030 ) << *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i) << "\n\t";
    Debug( 6030 ) << "\n";

    Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] vessels indices (unique map):\n";
    for( int i = 0; i < _M_epetra_vessels->getUniqueEpetra_Map()->NumMyElements(); ++i )
        Debug( 6030 ) << *(_M_epetra_vessels->getUniqueEpetra_Map()->MyGlobalElements() + i) << "\n\t";
    Debug( 6030 ) << "\n";

    int my_interface_values_array[my_interface_values.size()];
    num = 0;
		Debug( 6030 ) << "\n[OneDNet::OneDNet process " << _M_comm->MyPID()
									<< "] this list should be non repeated (interface values index set, size "
									<< my_interface_values.size() << ")";
    for( std::set<int>::iterator it = my_interface_values.begin();
    	it != my_interface_values.end(); ++it, ++num )
    {
    	Debug( 6030 ) << *it << " \t";
    	my_interface_values_array[num] = *it;
    }
    _M_epetra_map_interface_values.reset( new EpetraMap( -1,
                                                         my_interface_values.size(),
                                                         my_interface_values_array,
                                                         1, *_M_comm ) );
    _M_epetra_vector_interface_values.reset( new EpetraVector<double>( *_M_epetra_map_interface_values ) );
    // be careful, _M_epetra_interface_values has indexes starting from 1
    /*
     interface values for tube i are
     left A : _M_epetra_interface_values[4(i-1) + 1]
     left Q : _M_epetra_interface_values[4(i-1) + 2]
     right A : _M_epetra_interface_values[4(i-1) + 3]
     right Q : _M_epetra_interface_values[4(i-1) + 4]
    */
    
    Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] interface values (indices):\n";
    for( int i = 0; i < _M_epetra_map_interface_values->getRepeatedEpetra_Map()->NumMyElements(); ++i )
        Debug( 6030 ) << *(_M_epetra_map_interface_values->getRepeatedEpetra_Map()->MyGlobalElements() + i) << "\n\t";
    Debug( 6030 ) << "\n";
    
    // debug: Epetra
    for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
//        Debug( 6030 ) << *(_M_epetra_vessels->MyGlobalElements() + i) << "\n\t";

    		// Conversion using lexical_cast
      	try {
          // store edge index in string variable
          str=boost::lexical_cast<std::string>(*(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i));
      	}
      	catch (boost::bad_lexical_cast &) {
          // conversion failed
          std::cout << "\n[OneDNet::OneDNet proc " << _M_comm->MyPID() << "] lexical conversion failed!"
                    << std::endl;
      	}

      	// update prefix
      	prefix.replace( prefix_end - 3, prefix_end - 1, str );
      	// point to the correct section in data file
      	data_file.set_prefix( prefix.c_str() );

            // allocate an object of class 1D param
        typename OneDNet< SOLVER1D, PARAM1D >::OneDParamPtr
            param ( new PARAM1D( data_file ) );

        // attach 1D param pointer to network edge
        _M_network[ *this->Edge(*(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i)) ].onedparam = param;

        // allocate an object of class 1D solver
        typename OneDNet< SOLVER1D, PARAM1D >::OneDSolverPtr
            solver( new SOLVER1D( data_file, *param ) );

        // attach 1D solver pointer to network edge
        _M_network[ *this->Edge(*(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i)) ].onedsolver = solver;

        if( _M_set_time_param ) /* if true the network sets time data */
            {
                Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID()
                							<< "] 0- Setting time parameters for onedModelSolver "
                              << *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i) << "...\n";
                // set network time step in each solver
                solver->settimestep( _M_time_step );
                // set network end time in each solver
                solver->setendtime( _M_time_end );
                // set network initial time in each solver
                solver->setinittime( _M_time_beg );
                Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "]\t...done!\n";
            }

        /*
          Check if all solvers have the same time parameters.
          If not, just abort (the net is not (yet?) able to manage
          different time parameters)
        */
        if( solver->timestep() != _M_time_step )
            {
                std::cout << "\n[OneDNet::OneDNet process " << _M_comm->MyPID() << "]:\nWARNING! Different time step: "
                          << "oneDNet dt != onedModelSolver (number "
                          << *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i) << ") dt."
                          << " Not yet implemented!" << std::endl;
                abort();
            }

        if( solver->endtime() != _M_time_end )
            {
                std::cout << "\n[OneDNet::OneDNet process " << _M_comm->MyPID() << "]: WARNING! Different time end: "
                          << "oneDNet T != onedModelSolver (number "
                          << *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i) << ")."
                          << " Not yet implemented!" << std::endl;
                abort();
            }

        if( solver->inittime() != _M_time_beg )
            {
                std::cout << "\n[OneDNet::OneDNet process " << _M_comm->MyPID() << "]: WARNING! Different beginning time: "
                          << "oneDNet t0 != onedModelSolver (number "
                          << *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i) << ")."
                          << " Not yet implemented!" << std::endl;
                abort();
            }

        // Print out to screen info on 1D model
        param->showMeData(std::cout);
        solver->showMeData(std::cout);


        // reset prefix to initial value
        prefix = "1dnetwork/tubenn/";
        prefix_end = prefix.end();

    }
    // reset getpot pointer to root section
    prefix = "";
    data_file.set_prefix( prefix.c_str() );


    // auxiliary variables (number of edges connected to a vertex)
    Degree_Size degree_size, out_degree_size;
    // in-edges and out-edges lists for given interface
    std::vector<In_Edge_Iter> in_edges;
    std::vector<Out_Edge_Iter> out_edges;

    for( i=0; i<_M_num_interfaces; ++i){
        // get degree informations from boost::adjacency_list functions
        degree_size = degree(vd[i], _M_network);
        out_degree_size = out_degree(vd[i], _M_network);

        // visit the list of vertices and mark terminal ones
        switch( degree_size ){
            //      case 0:
            //	break;
        case 1: // only one connected edge -> not an internal vertex
            _M_network[ vd[i] ].internal = false;
            break;
        default: // internal vertex
        		Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] setting "
        									<< _M_network[ vd[i] ].index << " as an internal interface \n";
            _M_network[ vd[i] ].internal = true;
            in_edges = inEdges( _M_network[ vd[i] ].index );
            for( UInt ii = 0; ii < in_edges.size(); ++ii )
            	if( _M_epetra_vessels->getRepeatedEpetra_Map()->MyGID(_M_network[ *(in_edges[ii]) ].index) )
            	{
                Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] setting "
                							<< _M_network[ *(in_edges[ii]) ].index << " as an internal right BC vessel \n";
                _M_network[ *(in_edges[ii]) ].onedsolver->setBCRight_internalnode();
            	}
            out_edges = outEdges( _M_network[ vd[i] ].index );
            for( UInt oi = 0; oi < out_edges.size(); ++oi )
            	if( _M_epetra_vessels->getRepeatedEpetra_Map()->MyGID(_M_network[ *(out_edges[oi]) ].index) )
            	{
                Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] setting "
                							<< _M_network[ *(out_edges[oi]) ].index << " as an internal left BC vessel \n";
                _M_network[ *out_edges[oi] ].onedsolver->setBCLeft_internalnode();
            	}
        }
    }
    Debug( 6030 ) << "[OneDNet::OneDNet process " << _M_comm->MyPID() << "] leaving constructor\n";
}


template< class SOLVER1D, class PARAM1D >
typename OneDNet<SOLVER1D,PARAM1D>::Vertex_Iter
OneDNet<SOLVER1D,PARAM1D>::Vertex( const int& ind )
{
    // iterators to visit vertex list
    std::pair<Vertex_Iter, Vertex_Iter> vertex_iter_pair;

    for( vertex_iter_pair = vertices(_M_network);
         vertex_iter_pair.first != vertex_iter_pair.second;
         ++vertex_iter_pair.first){
        // check if the vertex has the requested index
        if ( _M_network[ *vertex_iter_pair.first ].index == ind )
            return *vertex_iter_pair.first; // return a Vertex_Iter object
    }
    // the function failed if it's printing out this!
    std::cout << "\n[OneDNet::Vertex] failed! Vertex " << ind << " not found..."
              << std::endl;
    abort();
}


template< class SOLVER1D, class PARAM1D >
typename OneDNet<SOLVER1D,PARAM1D>::Edge_Iter
OneDNet<SOLVER1D,PARAM1D>::Edge( const int& ind )
{
    // iterators to visit edge list
    std::pair<Edge_Iter, Edge_Iter> edge_iter_pair;

    for( edge_iter_pair = edges(_M_network);
         edge_iter_pair.first != edge_iter_pair.second;
         ++edge_iter_pair.first ){
        // check if the edge has the requested index
        if ( _M_network[ *edge_iter_pair.first ].index == ind )
            return edge_iter_pair.first; // return an Edge_Iter object
    }
    // the function failed if it's printing out this!
    std::cout << "\n[OneDNet::Edge] failed! Edge " << ind << " not found..."
              << std::endl;
    abort();
}


template< class SOLVER1D, class PARAM1D >
SOLVER1D &
OneDNet<SOLVER1D,PARAM1D>::Solver( const int& ind )
{
    return *_M_network[ *(this->Edge(ind)) ].onedsolver; // return a SOLVER1D object
}


template< class SOLVER1D, class PARAM1D >
typename OneDNet<SOLVER1D, PARAM1D>::OneDSolverPtr
OneDNet<SOLVER1D,PARAM1D>::SolverPtr( const int& ind )
{
    return _M_network[ *(this->Edge(ind)) ].onedsolver; // return a OneDSolverPtr object
}


template< class SOLVER1D, class PARAM1D >
PARAM1D &
OneDNet<SOLVER1D,PARAM1D>::Param( const int& ind )
{
    return *_M_network[ *(this->Edge(ind)) ].onedparam; // return a PARAM1D object
}


template< class SOLVER1D, class PARAM1D >
typename OneDNet<SOLVER1D, PARAM1D>::OneDParamPtr
OneDNet<SOLVER1D,PARAM1D>::ParamPtr( const int& ind )
{
    return _M_network[ *(this->Edge(ind)) ].onedparam; // return a OneDParamPtr object
}


template< class SOLVER1D, class PARAM1D >
std::vector<typename OneDNet<SOLVER1D,PARAM1D>::In_Edge_Iter>
OneDNet<SOLVER1D,PARAM1D>::inEdges( const int& ind )
{
    // return parameter
    std::vector<In_Edge_Iter> in_edge_iter_vec;
    // pair of in-edge iterators
    std::pair< In_Edge_Iter, In_Edge_Iter > in_edge_iter_pair;
    // visit vertex's in-edges list
    for( in_edge_iter_pair = in_edges( *this->Vertex(ind), _M_network);
         in_edge_iter_pair.first != in_edge_iter_pair.second;
         ++in_edge_iter_pair.first )
        {
            // add the index of each in-edge to the return variable
            in_edge_iter_vec.push_back( in_edge_iter_pair.first );
        }

    // return an empty vector if vertex not found (or vertex without in-edges)
    return in_edge_iter_vec;
}


template< class SOLVER1D, class PARAM1D >
std::vector<typename OneDNet<SOLVER1D,PARAM1D>::Out_Edge_Iter>
OneDNet<SOLVER1D,PARAM1D>::outEdges( const int& ind )
{
    // return parameter
    std::vector<Out_Edge_Iter> out_edge_iter_vec;
    // vertices() method returns a pair of vertex iterators (begin(), end())
    std::pair< Vertex_Iter, Vertex_Iter > vertex_iter_pair = vertices(_M_network);
    // pair of out-edge iterators
    std::pair< Out_Edge_Iter, Out_Edge_Iter > out_edge_iter_pair;
    // visit vertex list
    while( ( vertex_iter_pair.first != vertex_iter_pair.second ) )
        {
            // find vertex with index == ind
            if( _M_network[ *vertex_iter_pair.first ].index == ind )
                // visit vertex's out-edges list
                for( out_edge_iter_pair = out_edges(*vertex_iter_pair.first, _M_network);
                     out_edge_iter_pair.first != out_edge_iter_pair.second;
                     ++out_edge_iter_pair.first )
                    {
                        // add the index of each out-edge to the return variable
                        out_edge_iter_vec.push_back( out_edge_iter_pair.first );
                    }

            ++vertex_iter_pair.first;
        }
    // return an empty vector if vertex not found (or vertex without in-edges)
    return out_edge_iter_vec;
}


template< class SOLVER1D, class PARAM1D >
typename OneDNet<SOLVER1D,PARAM1D>::Degree_Size
OneDNet<SOLVER1D,PARAM1D>::outDegreeSize( const int& ind )
{
    // vertices() method returns a pair of vertex iterators (begin(), end())
    std::pair< Vertex_Iter, Vertex_Iter > vertex_iter_pair = vertices(_M_network);
    // visit vertex list
    while( ( vertex_iter_pair.first != vertex_iter_pair.second ) )
        {
            // find vertex with index == ind
            if( _M_network[ *vertex_iter_pair.first ].index == ind )
                // get out_degree information from boost::adjacency_list function
                return out_degree(*vertex_iter_pair.first, _M_network);

            ++vertex_iter_pair.first;
        }
    // the function failed if it's printing out this!
    std::cout << "\n[OneDNet::outDegreeSize] failed! Vertex "
              << ind << " not found..."
              << std::endl;
    abort();
}


template< class SOLVER1D, class PARAM1D >
typename OneDNet<SOLVER1D,PARAM1D>::Degree_Size
OneDNet<SOLVER1D,PARAM1D>::degreeSize( const int& ind )
{
    // vertices() method returns a pair of vertex iterators (begin(), end())
    std::pair< Vertex_Iter, Vertex_Iter > vertex_iter_pair = vertices(_M_network);
    // visit vertex list
    while( ( vertex_iter_pair.first != vertex_iter_pair.second ) )
        {
            // find vertex with index == ind
            if( _M_network[ *vertex_iter_pair.first ].index == ind )
                // get degree information from boost::adjacency_list function
                return degree(*vertex_iter_pair.first, _M_network);

            ++vertex_iter_pair.first;
        }
    // the function failed if it's printing out this!
    std::cout << "\n[OneDNet::degreeSize] failed! Vertex "
              << ind << " not found..."
              << std::endl;
    abort();
}


template< class SOLVER1D, class PARAM1D >
typename OneDNet<SOLVER1D,PARAM1D>::Network const&
OneDNet<SOLVER1D,PARAM1D>::network()
{
    return _M_network;
}


template< class SOLVER1D, class PARAM1D >
Real OneDNet<SOLVER1D,PARAM1D>::timestep() const
{
    return _M_time_step;
}


template< class SOLVER1D, class PARAM1D >
Real OneDNet<SOLVER1D,PARAM1D>::inittime() const
{
    return _M_time_beg;
}


template< class SOLVER1D, class PARAM1D >
Real OneDNet<SOLVER1D,PARAM1D>::endtime() const
{
    return _M_time_end;
}


// Call setBC() for each solver
template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::setBC( OneDBCFunctionPointer fun,
                                       const int& ind,
                                       std::string const& border,
                                       std::string const& line,
                                       std::string const& var )
{
    this->Solver( ind ).bcH().setBC( fun, border, line, var );
}


// Call initialize() for each solver
template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::initialize( GetPot data_file )
{
    // read data file (at each 1D model section, see also OneDNet constructor)
    std::string str, init_var;
    std::string prefix("1dnetwork/tubenn/");
    std::string::iterator prefix_end = prefix.end();

    // iterators to visit edge list
    std::pair<Edge_Iter, Edge_Iter> edge_iter_pair;

    int globalID;
    //    for( edge_iter_pair = edges(_M_network);
    //         edge_iter_pair.first != edge_iter_pair.second;
    //         ++edge_iter_pair.first){
    // debug: Epetra
    for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
        globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
        Debug( 6030 ) << "[OneDNet::initialize process " << _M_comm->MyPID() << "] initializing vessels:\n";
        Debug( 6030 ) << globalID << "\n";

        // Conversion using lexical_cast
        try {
            // store edge index in string variable
            str=boost::lexical_cast<std::string>(_M_network[ *this->Edge(globalID) ].index);
        }
        catch (boost::bad_lexical_cast &) {
            // conversion failed
            std::cout << "\n[OneDNet::initialize process " << _M_comm->MyPID() << "] lexical conversion failed!"
                      << std::endl;
        }

        // update prefix
        prefix.replace( prefix_end - 3, prefix_end - 1, str );
        // point to the correct section in data file
        data_file.set_prefix( prefix.c_str() );

        Debug( 6030 ) << "[OneDNet::initialize process " << _M_comm->MyPID() << "] 0- Initializing tube "
                      << _M_network[ *this->Edge(globalID) ].index << "\n";

        // call OneDModelSolver::initialize
        _M_network[ *this->Edge(globalID) ].onedsolver->initialize( data_file );
        //      solver->initialize(w1_0, w2_0, area0, 0, 0);
        //      solver->initialize(u1_0, u2_0);

        prefix = "1dnetwork/tubenn/";
    }

    // reset getpot pointer to root section
    data_file.set_prefix( "" );
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::savesol(  )
{
  for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
      int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
      Debug( 6030 ) << "[OneDNet::savesol process " << _M_comm->MyPID() << "] saving sol for tube:\n";
      Debug( 6030 ) << globalID << "\n";

        // call OneDModelSolver method
        _M_network[*this->Edge(globalID)].onedsolver->savesol();
    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::loadsol(  )
{
  for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
      int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
      Debug( 6030 ) << "[OneDNet::loadsol process " << _M_comm->MyPID() << "] loading sol for tube:\n";
      Debug( 6030 ) << globalID << "\n";

        // call OneDModelSolver method
        _M_network[*this->Edge(globalID)].onedsolver->loadsol();
    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::updateInterfaces( )
{
	_M_epetra_vector_interface_values->GlobalAssemble();

//	EpetraVector<double> passaggio(*_M_epetra_map_interface_values);
//	passaggio.Import( *_M_epetra_vector_interface_values, Insert );
//	EpetraVector<double> passaggio( *_M_epetra_vector_interface_values, _M_comm->MyPID() );
		EpetraVector<double> passaggio(*_M_epetra_vector_interface_values,
		                               *_M_epetra_map_interface_values->getRepeatedEpetra_Map());
		for( int i = 0; i < _M_epetra_map_interface_values->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
			int globalID = *(_M_epetra_map_interface_values->getRepeatedEpetra_Map()->MyGlobalElements() + i);
			Debug( 6030 ) << "\n[OneDNet::updateInterfaces process " << _M_comm->MyPID()
			<< "] passaggio[ " << globalID << " ] = " << passaggio[globalID]; } 
	// pair of out-edge iterators
//	std::pair<Out_Edge_Iter, Out_Edge_Iter> out_edge_iter_pair;
	// pair of in-edge iterators
//	std::pair<In_Edge_Iter, In_Edge_Iter> in_edge_iter_pair;
	// you expect this interface to have both in-edges and out-edges
	
//	for( int i = 0; i < _M_epetra_interfaces->NumMyElements(); ++i ) {
//		int globalID = *(_M_epetra_interfaces->MyGlobalElements() + i);
    _M_comm->Barrier();

		for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
			int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
			
//		for( in_edge_iter_pair = in_edges(*this->Vertex(globalID), _M_network);
//					in_edge_iter_pair.first != in_edge_iter_pair.second;
//					++in_edge_iter_pair.first )
//		{
			int tubeIndex = _M_network[ *this->Edge(globalID) ].index;

			Debug( 6030 ) << "[OneDNet::updateInterfaces process " << _M_comm->MyPID()
										<< "] 0- Updating tube " << tubeIndex //<< " at interface " << globalID;
										<< " with (indexes, value) " 
										<< 4*(tubeIndex-1) + 1 << ", "
//										<< _M_epetra_vector_interface_values->operator[](4*(tubeIndex-1) + 1) << "; "
										<< passaggio[4*(tubeIndex-1) + 1] << "; "
										<< 4*(tubeIndex-1) + 2 << ", "
//										<< _M_epetra_vector_interface_values->operator[](4*(tubeIndex-1) + 2) << "; "
										<< passaggio[4*(tubeIndex-1) + 2] << "; "
										<< 4*(tubeIndex-1) + 3 << ", "
//										<< _M_epetra_vector_interface_values->operator[](4*(tubeIndex-1) + 3) << "; "
										<< passaggio[4*(tubeIndex-1) + 3] << "; "
										<< 4*(tubeIndex-1) + 4 << ", "
//										<< _M_epetra_vector_interface_values->operator[](4*(tubeIndex-1) + 4) << "\n";
										<< passaggio[4*(tubeIndex-1) + 4] << "\n";

			_M_network[*this->Edge(globalID) ].onedsolver->setBCValuesLeft(
//			                  _M_epetra_vector_interface_values->operator[](4*(tubeIndex-1) + 1),
//			                  _M_epetra_vector_interface_values->operator[](4*(tubeIndex-1) + 2));
			                  passaggio[4*(tubeIndex-1) + 1],
							          passaggio[4*(tubeIndex-1) + 2]);
			_M_network[*this->Edge(globalID) ].onedsolver->setBCValuesRight(
//			                  _M_epetra_vector_interface_values->operator[](4*(tubeIndex-1) + 3),
//			                  _M_epetra_vector_interface_values->operator[](4*(tubeIndex-1) + 4));
			                                       passaggio[4*(tubeIndex-1) + 3],
			                                       passaggio[4*(tubeIndex-1) + 4] );

//		}
//		for( out_edge_iter_pair = out_edges(*this->Vertex(globalID), _M_network);
//					out_edge_iter_pair.first != out_edge_iter_pair.second;
//					++out_edge_iter_pair.first )
//		{
//			int tubeIndex = _M_network[ *(out_edge_iter_pair).first ].index;
//
//			Debug( 6030 ) << "[OneDNet::updateInterfaces process " << _M_comm->MyPID()
//										<< "] 0- Updating tube " << tubeIndex << " at interface " << globalID;
//			Debug( 6030 ) << " with values (indexes) " << 4*(tubeIndex-1) + 1
//										<< ", "  << passaggio[4*(tubeIndex-1) + 1]
//										<< "; " << 4*(tubeIndex-1) + 2
//										<< ", "  << passaggio[4*(tubeIndex-1) + 2] << "\n";
//
//			_M_network[*(out_edge_iter_pair).first ].onedsolver->setBCValuesLeft(
//          passaggio[4*(tubeIndex-1) + 1],
//          passaggio[4*(tubeIndex-1) + 2]);
//
//		}
	}
    *(_M_epetra_vector_interface_values) *= 0; //->getEpetraVector().ReplaceGlobalValues( 2, comp, val );
}
    
template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::timeAdvance( const Real& time_val )
{
		// impose interface conditions at internal vertices
    computeInterfaceTubesValues();
    
    updateInterfaces();

    for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
        int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
        Debug( 6030 ) << "[OneDNet::timeAdvance process " << _M_comm->MyPID() << "]\n";

        // tell me what I am doing
        Debug( 6030 ) << "[OneDNet::timeAdvance process " << _M_comm->MyPID() << "] 0- Time advancing tube "
        							<< _M_network[ *this->Edge(globalID) ].index << "\n";
        // call OneDModelSolver method
        _M_network[*this->Edge(globalID)].onedsolver->timeAdvance(time_val);

    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::iterate( const Real& time_val , const int& count)
{
  	for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
      	int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);

        // tell me what I am doing
        Debug( 6030 ) << "[OneDNet::iterate process " << _M_comm->MyPID() << "] 0- Iterating tube "
        							<< _M_network[ *this->Edge(globalID) ].index << "\n";
        // call OneDModelSolver method
        _M_network[ *this->Edge(globalID) ].onedsolver->iterate(time_val, count);

    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::postProcess( const Real& time_val )
{
		Debug( 6030 ) << "[OneDNet::postProcess process " << _M_comm->MyPID() << "] inside postprocessing ";
		for( int i = 0; i < _M_epetra_vessels->getUniqueEpetra_Map()->NumMyElements(); ++i ) {
    		int globalID = *(_M_epetra_vessels->getUniqueEpetra_Map()->MyGlobalElements() + i);

        // tell me what I am doing
        Debug( 6030 ) << "[OneDNet::postProcess process " << _M_comm->MyPID() << "] 0- Postprocessing tube "
        							<< _M_network[ *this->Edge(globalID) ].index << "\n";
        // call OneDModelSolver method
        _M_network[ *this->Edge(globalID) ].onedsolver->postProcess(time_val);
    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::openFileBuffers()
{
		for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
  			int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
  			Debug( 6030 ) << "[OneDNet::openFileBuffers process " << _M_comm->MyPID() << "]\n";

        // call OneDModelSolver method
        _M_network[ *this->Edge(globalID) ].onedsolver->openFileBuffers();
    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::output2FileBuffers( const Real& time_val )
{
		for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
				int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
				Debug( 6030 ) << "[OneDNet::output2FileBuffers process " << _M_comm->MyPID() << "]\n";

        // call OneDModelSolver method
        _M_network[ *this->Edge(globalID) ].onedsolver->output2FileBuffers( time_val );
    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::resetFileBuffers()
{
		for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
				int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
				Debug( 6030 ) << "[OneDNet::resetFileBuffers process " << _M_comm->MyPID() << "]\n";

        // call OneDModelSolver method
        _M_network[ *this->Edge(globalID) ].onedsolver->resetFileBuffers();
    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::seekpFileBuffers()
{
		for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
				int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
				Debug( 6030 ) << "[OneDNet::seekpFileBuffers process " << _M_comm->MyPID() << "]\n";

        // call OneDModelSolver method
        _M_network[ *this->Edge(globalID) ].onedsolver->seekpFileBuffers();
    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::tellpFileBuffers()
{
		for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
				int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
				Debug( 6030 ) << "[OneDNet::tellpFileBuffers process " << _M_comm->MyPID() << "]\n";

        // call OneDModelSolver method
        _M_network[ *this->Edge(globalID) ].onedsolver->tellpFileBuffers();
    }
}


template< class SOLVER1D, class PARAM1D >
void OneDNet<SOLVER1D,PARAM1D>::closeFileBuffers()
{
		for( int i = 0; i < _M_epetra_vessels->getRepeatedEpetra_Map()->NumMyElements(); ++i ) {
				int globalID = *(_M_epetra_vessels->getRepeatedEpetra_Map()->MyGlobalElements() + i);
				Debug( 6030 ) << "[OneDNet::closeFileBuffers process " << _M_comm->MyPID() << "]\n";

        // call OneDModelSolver method
        _M_network[ *this->Edge(globalID) ].onedsolver->closeFileBuffers();
    }
}


template< class SOLVER1D, class PARAM1D >
void
OneDNet<SOLVER1D,PARAM1D>::computeInterfaceTubesValues( )
{
//    // vertex iterators (for visiting the graph)
//    std::pair< Vertex_Iter, Vertex_Iter > vertex_iter_pair;
//    // visit vertex list
//    for( vertex_iter_pair = vertices(_M_network);
//         vertex_iter_pair.first != vertex_iter_pair.second;
//         ++vertex_iter_pair.first ){

		for( int i = 0; i < _M_epetra_interfaces->NumMyElements(); ++i ) {
 				int globalID = *(_M_epetra_interfaces->MyGlobalElements() + i);

 				Debug( 6030 ) << "[OneDNet::computeInterfaceTubesValues process "
 											<< _M_comm->MyPID() << "] 0- Computing Interface "
                      << _M_network[ *this->Vertex(globalID) ].index << "\n";
        /*
          In principle, it is possible to implement different behaviour for
          interfaces (e. g. energy dissipation due to branching angles,
          vasculare valves etc). (not done yet)
        */
        if( _M_network[ *this->Vertex(globalID) ].internal ){ // internal interface
            switch( _M_network[ *this->Vertex(globalID) ].type ){
            case 0:
                /*
                  inflow tube: no "interface" conditions
                  (actual boundary conditions instead)
                */
                break;
            case 1:
                interface_continuity_conditions( this->Vertex(globalID) );
                break;
            case 99:
                /*
                  outflow tube: no "interface" conditions
                  (actual boundary conditions instead)
                */
                break;

                // other cases can be added here!

            default:
                std::cout << "\n[OneDNet::computeInterfaceTubesValues process " << _M_comm->MyPID() << "] Unknown type "
                          << _M_network[*this->Vertex(globalID)].type
                          << " for vertex " << _M_network[*this->Vertex(globalID)].index
                          << std::endl;
                break;
            } // switch
        } // if
    } // for
}


template< class SOLVER1D, class PARAM1D >
void
OneDNet<SOLVER1D,PARAM1D>::interface_continuity_conditions( Vertex_Iter const& vertex )
{
    // pair of out-edge iterators
    std::pair<Out_Edge_Iter, Out_Edge_Iter> out_edge_iter_pair;
    // pair of in-edge iterators
    std::pair<In_Edge_Iter, In_Edge_Iter> in_edge_iter_pair;
    // list of edges connected to the considered vertex
    MapSolver interfaceTubes;
    // map identifying in-edges(+, true) and out-edges(-, false)
    std::map< int, bool > signum;
    // map identifying edges indexes
    std::map< int, int > index;
    // edges index
    UInt i( 0 );
    // number of iterations
    UInt niter(0);
    // vector and matrix dimension
    UInt f_size;
    // unknown of non linear equation f(x) = 0
    Vector x;
    // non linear function f
    Vector f;
    // jacobian of the non linear function
    Matrix jac;
    // transpose of the jacobian of the non linear function
    Matrix jac_trans;
    // tmp matrix for lapack lu inversion
    boost::numeric::ublas::vector<Int> ipiv;
    // lapack variable
    int INFO[1];
    int NBRHS[1]; // nb columns of the rhs := 1.
    int NBU[1];
    
    // take into account interface type

    // each edge (1D solver) needs 2 boundary conditions
    f_size = 2 * degree(*vertex, _M_network);

    // create the map for interface EpetraVector
//    int _vector_map[f_size];
    
    // unknown vector
    /*
      x contains the (unknown) boundary values for edges connected
      to the considered interface.
    */
    x.resize( f_size );
    x.clear();
    Debug( 6030 ) << "[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
    							<< "] before newton iterations x =\n" ;

    // you expect this interface to have both in-edges and out-edges
    for( in_edge_iter_pair = in_edges(*vertex, _M_network);
         in_edge_iter_pair.first != in_edge_iter_pair.second;
         ++in_edge_iter_pair.first )
        {
    				// this piece of code is executed only for my vertices, for which I have all inedges and outedges 
    				
    				// add edges to the list
            interfaceTubes.insert
                ( MapSolverValueType(i, _M_network[*(in_edge_iter_pair).first ].onedsolver ) );
            // in-edges have "+" signum (boolean value true)
            signum.insert( std::map< int, bool >::value_type(i, true) );
            index.insert( std::map< int, int >::value_type(i, _M_network[*(in_edge_iter_pair).first ].index) );
            Debug( 6030 ) << "[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
            							<< "] \ttube " << _M_network[*(in_edge_iter_pair).first ].index
                          << " has signum " << signum[i] << " (map index) " << i << "\n";
            i++;
        }
    for( out_edge_iter_pair = out_edges(*vertex, _M_network);
         out_edge_iter_pair.first != out_edge_iter_pair.second;
         ++out_edge_iter_pair.first )
        {
            // add edges to the list
            interfaceTubes.insert
                ( MapSolverValueType(i, _M_network[ *(out_edge_iter_pair).first ].onedsolver ) );
            // in-edges have "-" signum (boolean value false)
            signum.insert( std::map< int, bool >::value_type(i, false) );
            index.insert( std::map< int, int >::value_type(i, _M_network[*(out_edge_iter_pair).first ].index) );
            Debug( 6030 ) << "[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
            							<< "] \ttube " << _M_network[*(out_edge_iter_pair).first ].index
                          << " has signum " << signum[i] << " (map index) " << i << "\n";
            i++;
        }

//    M_matrMass.reset(new EpetraVector(M_localMap));

    // i goes from 0 to the number of considered edges (n_edges)
    for( i = 0; i < interfaceTubes.size(); ++i )
        {
            // first value: A at the boundary for i-th edge
            x[2*i] = signum[i] ? interfaceTubes[i]->BCValuesRight().first :
                interfaceTubes[i]->BCValuesLeft().first;
            // second value: Q at the boundary for i-th edge
            x[2*i+1] = signum[i] ? interfaceTubes[i]->BCValuesRight().second :
                interfaceTubes[i]->BCValuesLeft().second;
            Debug( 6030 ) << "[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
            							<< "] \tcomponent " << 2*i << " " << x[2*i] << "\n";
            Debug( 6030 ) << "[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
            							<< "] \tcomponent " << 2*i+1 << " " << x[2*i+1] << "\n";
        }
    /*
      prepare the data structures needed for solving non linear equations
      (continuity + compatibility)
    */
    f.resize( f_size );
    f.clear();
    jac.resize( f_size, f_size );
    jac.clear();
    jac_trans.resize( f_size, f_size );
    jac_trans.clear();
    ipiv.resize( f_size );
    ipiv.clear();
    INFO[0] = 0;
    NBRHS[0] = 1; // nb columns of the rhs := 1.
    NBU[0] = f_size;
    i=0;
    // newton raphson iteration
    do
        {
            // fill f and its jacobian matrix
            f_jac( x, f, jac, interfaceTubes, signum );
            // transpose to pass to fortran storage (lapack!)
            jac_trans = trans(jac);
            // Compute f <-  ( df(x)^{-1} f(x) ) (lu dcmp)
            dgesv_(NBU, NBRHS, &jac_trans(0,0), NBU , &ipiv(0), &f(0), NBU, INFO);
            ASSERT_PRE(!INFO[0],"Lapack LU resolution of y = df(x)^{-1} f(x) is not achieved.");
            // x = x - df(x)^{-1} f(x)
            x += - f;
        }
    // convergence check
    while( (std::fabs( f(0) ) > 1e-12) && (++niter < 100) );
    /*
      f(0) contains the mass flow balance: by minimizing its value we want
      to ensure mass conservation.
      Moreover, we expect a low number of iterations, since the initial
      guess is given from the boundary conditions at previous time step,
      which are likely to be a good estimation of the current solution
      (time step is small and we expect solutions to be continuous
      in time).
    */
    // no convergence
    if( niter == 100 )
        {
            std::cout << "\nOneDNet::computeInterfaceTubesValues process " << _M_comm->MyPID() << ": newton iterations"
                      << " did not converge after " << niter << " iterations."
                      << " Expected tolerance for mass flow balance = " << 1e-12
                      << ", error = " << std::fabs( f(0) )
                      << ", interface = " << _M_network[*vertex].index
                      << std::endl;
            abort();
        }

    Debug( 6030 ) << "[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
            							<< "] after newton iterations x =\n" ;
    // i goes from 0 to the number of considered edges (n_edges)
    for( i = 0; i < interfaceTubes.size(); ++i )
        {
            Debug( 6030 ) << "[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
            							<< "] \tcomponent " << 2*i << " " << x[2*i] << "\n";
            Debug( 6030 ) << "[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
            							<< "] \tcomponent " << 2*i+1 << " " << x[2*i+1] << "\n";
        }

    // be careful, _M_epetra_interface_values has indexes starting from 1
    /*
     interface values for tube i are
     left A : _M_epetra_interface_values[4(i-1) + 1]
     left Q : _M_epetra_interface_values[4(i-1) + 2]
     right A : _M_epetra_interface_values[4(i-1) + 3]
     right Q : _M_epetra_interface_values[4(i-1) + 4]
    */
    
// 		EpetraVector<double> passaggio(*_M_epetra_map_interface_values);

		for( int i = 0; i < _M_epetra_vector_interface_values->Map().NumMyElements(); ++i ) {
 				Debug( 6030 ) << "\n[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
				<< "] show epetra vector interface values "
				<< *(_M_epetra_vector_interface_values->Map().MyGlobalElements() + i); }

    // set boundary conditions to 1D solvers
    for( i = 0; i < interfaceTubes.size(); ++i )
        {
			int comp[2];
			double val[2];
            // in-edges need conditions to the right boundary
            if( signum[i] )
                {
//                    interfaceTubes[i]->setBCValuesRight( x[2*i], x[2*i+1] );
            				comp[0] = 4*(index[i]-1) + 3;
            				comp[1] = 4*(index[i]-1) + 4;
            				val[0] = x[2*i];
            				val[1] = x[2*i+1];
                    Debug( 6030 ) << "\n[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
                    							<< "] modifying epetra vector interface for tube "
                    							<< index[i] << " (component, value) "
                    							<< comp[0] << ", " << val[0]
//                    Debug( 6030 ) << _M_epetra_vector_interface_values->operator()(4*(index[i]-1) + 3)
                    							<< "; " << comp[1] << ", " << val[1]; // << ", "
//                    							<< _M_epetra_vector_interface_values->operator()(4*(index[i]-1) + 4) << "\n";
                }
            // out-edges need conditions to the left boundary
            else
                {
//                    interfaceTubes[i]->setBCValuesLeft( x[2*i], x[2*i+1] );
            				comp[0] = 4*(index[i]-1) + 1;
            				comp[1] = 4*(index[i]-1) + 2;
            				val[0] = x[2*i];
            				val[1] = x[2*i+1];
                    Debug( 6030 ) << "\n[OneDNet::interface_continuity_conditions process " << _M_comm->MyPID()
                    							<< "] modifying epetra vector interface for tube "
                    							<< index[i] << " (component, value) "
                    							<< comp[0] << ", " << val[0]
//                    Debug( 6030 ) << _M_epetra_vector_interface_values->operator()(4*(index[i]-1) + 1)
                    							<< "; " << comp[1] << ", " << val[1]; // << ", "
//                    							<< _M_epetra_vector_interface_values->operator()(4*(index[i]-1) + 2) << "\n";
//                    _M_epetra_vector_interface_values->operator()(4*(index[i]-1) + 1) = x[2*i];
//                    _M_epetra_vector_interface_values->operator()(4*(index[i]-1) + 2) = x[2*i+1];
                }
          _M_epetra_vector_interface_values->getEpetraVector().ReplaceGlobalValues( 2, comp, val );
        }
//    _M_epetra_vector_interface_values->Import( passaggio, Insert );
}


template< class SOLVER1D, class PARAM1D >
void
OneDNet<SOLVER1D,PARAM1D>::f_jac( const Vector& x, Vector& f, Matrix& jac,
                                  typename OneDNet<SOLVER1D,PARAM1D>::MapSolver& intTubes,
                                  std::map< int, bool >& sign )
{
    // eigen values of the jacobian diffFlux (= dF/dU)
    Real  eigval1, eigval2, eigval;
    // left eigen vectors for the eigen values eigval1 and eigval2
    Vec2D left_eigvec1, left_eigvec2, left_eigvec;
    // right hand side for the 2x2 linear system to be solved for each 1D solver
    Real rhsBC = 0.;
    // quasi linear source term
    Vec2D qlSource;
    // value of U at the boundary, at the neighboring internal node
    //            and at the foot of the characteristics line
    Vec2D U_boundary, U_internalBd, U_charact_pt;
    // axial coordinate of boundary and neighboring internal node
    Real  boundaryPoint, internalBdPoint;
    // number of degrees of freedom
    UInt dof;

    /***** f(0): non linear equation
           Mass flow conservation    *****/
    f( 0 ) = 0.;
    for( UInt i = 0; i < intTubes.size(); ++i )
        {
            /*
              in-edges give a positive contribution to mass flow balance;
              out-edges give a negative contribution
            */
            f( 0 ) += sign[i] ? x[2*i+1] : ( - x[2*i+1] );
            // the derivative if trivial (+1 for in-edges, -1 for out-edges)
            jac( 0, 2*i+1 ) =  sign[i] ? 1 : ( - 1 );
        }

    /***** f(1) ... f(n_edges-1): set of non linear equations
           Total pressure conservation                        *****/
    for( UInt i = 1; i < intTubes.size(); ++i )
        {
            /*
              the total pressure computed from each 1D model in the interface list
              (starting from the second) is compared to which
              computed from the first one
            */
            f(i) = ( intTubes[i]->oneDParam().totalPressure(
                                                            x[2*i],  x[2*i+1], sign[i] ?
                                                            intTubes[i]->RightNodeId() :
                                                            intTubes[i]->LeftNodeId() )
                     - intTubes[0]->oneDParam().totalPressure(
                                                              x[0], x[1], sign[0] ?
                                                              intTubes[0]->RightNodeId() :
                                                              intTubes[0]->LeftNodeId() ) );

            // Jacobian (exploit Lifev::vectorFunction1D methods)
            // df_i / dA_0:
            jac( i, 0 ) = - intTubes[0]->oneDParam().totalPressureDiff(
                                                                       x[0], x[1], 1,  sign[0] ?
                                                                       intTubes[0]->RightNodeId() :
                                                                       intTubes[0]->LeftNodeId() );
            // df_i / dQ_0:
            jac( i, 1 ) = - intTubes[0]->oneDParam().totalPressureDiff(
                                                                       x[0], x[1], 2,  sign[0] ?
                                                                       intTubes[0]->RightNodeId() :
                                                                       intTubes[0]->LeftNodeId() );
            // df_i / dA_i:
            jac( i, 2*i ) = intTubes[i]->oneDParam().totalPressureDiff(
                                                                       x[2*i], x[2*i+1], 1,
                                                                       sign[i] ?
                                                                       intTubes[i]->RightNodeId() :
                                                                       intTubes[i]->LeftNodeId() );
            // df_i / dQ_i:
            jac( i, 2*i+1 ) = intTubes[i]->oneDParam().totalPressureDiff(
                                                                         x[2*i], x[2*i+1], 2,
                                                                         sign[i] ?
                                                                         intTubes[i]->RightNodeId() :
                                                                         intTubes[i]->LeftNodeId() );
        }

    /***** f(n_edges) ... f(2*n_edges -1): set of non linear equations
           Compatibility conditions                                    *****/
    // consider edges in the interface list
    for( UInt i = 0; i < intTubes.size(); ++i )
        {
            // in-edges --> right boundary
            if( sign[i] )
                {
                    // store coordinates of boundary and internal neighbouring node
                    boundaryPoint = intTubes[i]->RightEdge().pt2().x();
                    internalBdPoint = intTubes[i]->RightEdge().pt1().x();
                    // store solution values at those nodes
                    U_boundary   = Vec2D ( intTubes[i]->BCValuesRight() );
                    U_internalBd = Vec2D ( intTubes[i]->BCValuesInternalRight() );
                    // store the label of the boundary node
                    dof = intTubes[i]->RightNodeId();
                }
            // out-edges --> right boundary
            else
                {
                    // store coordinates of boundary and internal neighbouring node
                    boundaryPoint = intTubes[i]->LeftEdge().pt1().x();
                    internalBdPoint = intTubes[i]->LeftEdge().pt2().x();
                    // store solution values at those nodes
                    U_boundary   = Vec2D ( intTubes[i]->BCValuesLeft() );
                    U_internalBd = Vec2D ( intTubes[i]->BCValuesInternalLeft() );
                    // store the label of the boundary node
                    dof = intTubes[i]->LeftNodeId();
                }
            // compute eigenvalues and eigenvectors at boundary node
            intTubes[i]->FluxFun().jacobian_EigenValues_Vectors(
                                                                U_boundary.first,
                                                                U_boundary.second,
                                                                eigval1, eigval2,
                                                                left_eigvec1.first,
                                                                left_eigvec1.second,
                                                                left_eigvec2.first,
                                                                left_eigvec2.second,
                                                                dof);
            ASSERT( eigval1 > 0. && eigval2 < 0. ,
                    "The eigenvalues do not have the expected signs.");
            /*
              Find which eigenvalue/eigenvector is associated to the
              characteristic variable exiting from the domain at
              the interface boundary
            */
            if( sign[i] )
                {
                    // in-edges: first eigenvalue/eigenvector
                    eigval = eigval1;
                    left_eigvec = left_eigvec1;
                }
            else
                {
                    // out-edges: second eigenvalue/eigenvector
                    eigval = eigval2;
                    left_eigvec = left_eigvec2;
                }
            // find solution at the foot of characteristic line (by linear interpolation)
            U_charact_pt = interpolLinear(boundaryPoint, internalBdPoint,
                                          _M_time_step, eigval,
                                          U_boundary, U_internalBd);
            // find the pseudo-characteristic associated to U_charact_pt
            rhsBC = dot( left_eigvec , U_charact_pt );
            // compute the (linearized) source term of Euler equations
            qlSource.first = intTubes[i]->SourceFun().QuasiLinearSource(
                                                                        U_charact_pt.first,
                                                                        U_charact_pt.second,
                                                                        1, dof);
            qlSource.second = intTubes[i]->SourceFun().QuasiLinearSource(
                                                                         U_charact_pt.first,
                                                                         U_charact_pt.second,
                                                                         2, dof);
            // extrapolate the pseudo-characteristic along (exiting) characteristic line
            rhsBC -= _M_time_step * dot( left_eigvec , qlSource );

            // f( n_edges+i ): compatibility condition
            f( intTubes.size()+i ) = left_eigvec.first * x[2*i]
                + left_eigvec.second * x[2*i+1] - rhsBC;

            // Jacobian
            // the only non-zero values are the following
            jac( intTubes.size()+i, 2*i ) =  left_eigvec.first; //< df_(n_edges-1+i)/dA_(2*i)
            jac( intTubes.size()+i, 2*i+1 ) =  left_eigvec.second;//< df_(n_edges-1+i)/dQ_(2*1+i)
        }
}


template< class SOLVER1D, class PARAM1D >
Real OneDNet<SOLVER1D,PARAM1D>::dot(const Vec2D& vec1, const Vec2D& vec2) const
{
    // scalar product of 2D vectors
    return vec1.first * vec2.first + vec1.second * vec2.second;
}


template< class SOLVER1D, class PARAM1D >
typename OneDNet<SOLVER1D,PARAM1D>::Vec2D
OneDNet<SOLVER1D,PARAM1D>::interpolLinear(const Real& point_bound,
                                          const Real& point_internal,
                                          const Real& deltaT,
                                          const Real& eigenvalue,
                                          const Vec2D& U_bound,
                                          const Vec2D& U_intern) const
{
    // size of space interval
    Real deltaX = std::abs(point_bound - point_internal);
    // cfl number
    Real cfl =  eigenvalue * deltaT / deltaX;
    // weight in the linear approximation
    Real weight;

    if ( point_bound < point_internal ) { // the edge is on the left of the domain
        ASSERT( -1. < cfl && cfl < 0. ,
                "This characteristics is wrong!\nEither it is not outcoming \
                (eigenvalue>0 at the left of the domain),\n or CFL is too high.");

        weight = - cfl;
    }
    else {  // the edge is on the right of the domain
        ASSERT( 0. < cfl && cfl < 1. ,
                "This characteristics is wrong!\nEither it is not outcoming \
                (eigenvalue<0 at the right of the domain),\n or CFL is too high.");

        weight = cfl;
    }
    // convex linear combination
    Vec2D u_interp( ( 1 - weight ) * U_bound.first  + weight * U_intern.first ,
                    ( 1 - weight ) * U_bound.second + weight * U_intern.second );
    return u_interp;
}

}

#endif
