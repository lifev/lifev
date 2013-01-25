//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief File constaining a solver for the level set equation

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 05-11-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    A more detailed description of the file (if necessary)
 */

#ifndef LEVELSETSOLVER_H
#define LEVELSETSOLVER_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <lifev/core/algorithm/SolverAztecOO.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/solver/ADRAssemblerIP.hpp>

#include <lifev/level_set/solver/LevelSetData.hpp>

#include <vector>
#include <limits>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

// LevelSetSolver - This class contains a solver for the level set problem usually created when dealing with free surface problems.
/*!

  <b> Scope </b>

  For solving free surface problems (free surface flows for example), there exists many methods to
  track or capture the moving interface. The level set method, popularized by Osher and Sethian, is
  one of them.

  One has to deal with a scalar function \f$ \phi \f$ such that the interface is described by the
  zero level set of \f$ \phi \f$, i.e. the interface is \f$ \{\mathbf{x} | \phi(\mathbf{x})=0 \} \f$.

  Two key elements are then needed to make the interface evolve. The first one is the PDE for \f$ \phi \f$:

  \f[ \frac{\partial \phi}{\partial t} + \beta \cdot \nabla \phi = 0 \f]

  The second one is the reinitialization, essential to keep the distance property of $\f \phi \f$ and
  so the accuracy on the zero level set location.

  The scope of this class is to provide both the time and space discretization of the level set PDE
  (and the stabilization needed by this pure transport equation), as well as some tools inherants to the
  level set method, including reinitialization procedures.

  <b> PDE discretization </b>

  For the time discretization, a BDF scheme is used. Space discretization
  is handled in the ADR framework (see ADRAssembler) with finite elements.
  As the standard approximation leads to instabilities, this solver gives
  the possibility to use stabilized schemes.

  Right now, only IP stabilization is enabled (see ADRAssemblerIP).

  The LevelSetData structure provides the convenient way to set all the
  parameters for both discretizations.

  <b> Reinitialization </b>

  In the moment, only direct reinitialization using geometric computations
  of all the distances is available. This makes the choice of the element
  used for the space discretization to be restricted to the P1.

  <b> Usage </b>

  The best usage consists in, first of all, build a LevelSetData stucture

  \code
  boost::shared_ptr<DataLevelSet> data_level_set(new DataLevelSet);
  data_level_set->setup(...);
  \endcode

  Then, the level set solver can be defined and setup:

  \code
  LevelSetSolver<mesh_type> level_set(...);
  level_set.setup(data_level_set);
  \endcode

  Do not forget to initialize the solver with the initial data!

  \code
  level_set.initialize(...);
  \endcode

  Finally, before starting the time loop, you have to setup the solver

  \code
  level_set.setupLinearSolver(...);
  \endcode

  In the time loop, the most important things are:

  \code
  data_level_set->dataTime()->updateTime(); // Update the time in the data
  level_set.updateSystem(...);              // Update the linear system
  level_set.iterate();                      // Solve the system
  level_set.reinitializationDirect();       // Reinitialize if needed
  \endcode


  <b> Futur improvements </b>

  The current reinitialization is quite primitive: it is
  not very accurate, volume is not very well conserved
  and the computational cost is far too high (it consists
  more or less to a broadcast of the mesh...). Better
  reinitialization procedures should be added.

    @author Samuel Quinodoz
    @version 1.0
*/

template< typename mesh_type, typename solver_type = LifeV::SolverAztecOO>
class LevelSetSolver
{
public:

    //! @name Public Types
    //@{

    typedef MapEpetra map_type;

    typedef FESpace<mesh_type, map_type> fespace_type;
    typedef boost::shared_ptr<fespace_type>              fespace_ptrType;

    typedef typename solver_type::vector_type vector_type;

    typedef typename solver_type::matrix_type matrix_type;
    typedef boost::shared_ptr<matrix_type> matrix_ptrType;

    typedef DataLevelSet data_type;
    typedef boost::shared_ptr<data_type> data_ptrType;

    typedef TimeAdvanceBDF<vector_type> bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrType;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty constructor
    LevelSetSolver( fespace_ptrType fespace, fespace_ptrType betaFESpace );

    //! Destructor
    ~LevelSetSolver() {};

    //@}

    //! @name Methods
    //@{

    //! Set the internal data with this data
    void setup(const data_ptrType& data);

    //! Set the initial solution
    /*!
      @Warning: the setup method has to be called before
     */
    void initialize(const vector_type& init);

    //! Rebuild the system
    /*!
      @Warning: the initialize method has to be called before
     */
    void updateSystem(const vector_type& beta, BCHandler& bch, const Real& time);

    //! Setup the linear solver
    void setupLinearSolver(const GetPot& dataFile, const string& section = "level-set");

    //! Solve the problem with the built system
    void iterate();

    //! Reinitialization via direct geometrical computations
    void reinitializationDirect();

    //@}

    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    //! Return the map of the unknown FE space
    inline map_type map() const { return M_fespace->map(); }

    //! Return the solution
    inline vector_type solution() const { return M_solution; }

    //! Return the data
    inline data_ptrType data() const {return M_data; }

    //@}

private:

    // Private typedefs

    // Assemblers
    typedef ADRAssembler<mesh_type,matrix_type,vector_type> adrAssembler_type;
    typedef ADRAssemblerIP<mesh_type,matrix_type,vector_type> ipAssembler_type;

    // Geometry
    typedef std::vector<Real>       point_type;
    typedef std::vector<point_type> face_type;


    //! @name Private Methods
    //@{

    void updateFacesNormalsRadius();
    Real computeUnsignedDistance(const std::vector< Real >& point);
    void cleanFacesData();
    inline Real distanceBetweenPoints(const point_type& P1, const point_type& P2) const
    {
        return std::sqrt((P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1]) + (P1[2]-P2[2])*(P1[2]-P2[2]));
    };
    inline point_type crossProd(const point_type& P1, const point_type& P2) const
    {
        point_type n(3,0);

        n[0]=P1[1]*P2[2] - P1[2]*P2[1];
        n[1]=P1[2]*P2[0] - P1[0]*P2[2];
        n[2]=P1[0]*P2[1] - P1[1]*P2[0];

        return n;
    };
    inline Real scalProd(const point_type& P1, const point_type& P2) const
    {
        return P1[0]*P2[0] + P1[1]*P2[1] + P1[2]*P2[2];
    };
    inline Real norm(const point_type& V) const
    {
        return std::sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
    };

    //@}


    // Private members
    fespace_ptrType M_fespace;
    fespace_ptrType M_betaFESpace;

    vector_type M_solution;

    adrAssembler_type M_adrAssembler;
    ipAssembler_type M_ipAssembler;

    matrix_ptrType M_systemMatrix;
    matrix_ptrType M_massMatrix;
    matrix_ptrType M_rhsMatrix;

    vector_type M_rhs;

    data_ptrType M_data;

    bdf_ptrType M_bdf;

    solver_type M_linearSolver;

    std::vector<face_type> M_faces;
    std::vector<point_type> M_normals;
    std::vector<Real> M_radius;

};

// ===================================================
// Constructors & Destructor
// ===================================================

template< typename mesh_type, typename solver_type>
LevelSetSolver<mesh_type,solver_type>::
LevelSetSolver( fespace_ptrType fespace, fespace_ptrType betaFESpace )
        :
        M_fespace(fespace),
        M_betaFESpace(betaFESpace),

        M_solution(fespace->map()),

        M_adrAssembler(),
        M_ipAssembler(),

        M_systemMatrix(),
        M_massMatrix(),
        M_rhsMatrix(),

        M_rhs(fespace->map()),

        M_data(),

        M_bdf(),

        M_linearSolver()
{
    M_adrAssembler.setup(fespace,betaFESpace);
    M_ipAssembler.setup(fespace,betaFESpace);
    M_linearSolver.setCommunicator(fespace->map().commPtr());
}

// ===================================================
// Methods
// ===================================================

template< typename mesh_type, typename solver_type>
void
LevelSetSolver<mesh_type,solver_type>::
setup(const data_ptrType& data)
{
    M_data=data;
    M_bdf.reset(new bdf_type());
    M_bdf->setup(data->dataTimeAdvance()->orderBDF());
}


template< typename mesh_type, typename solver_type>
void
LevelSetSolver<mesh_type,solver_type>::
initialize(const vector_type& init)
{
    ASSERT(M_bdf!=0, "Bdf structure not initialized. Use the setup method before initialize.");
    M_solution = init;
    M_bdf->setInitialCondition(init);
}

template< typename mesh_type, typename solver_type>
void
LevelSetSolver<mesh_type,solver_type>::
updateSystem(const vector_type& beta, BCHandler& bcHandler, const Real& time)
{
    ASSERT(M_bdf!=0, "Bdf structure not initialized. Use the setup method before updateSystem.");
    ASSERT(M_data!=0, "No data available. Use the setup method before updateSystem.");

    // Check the existence of the matrices
    if (M_massMatrix == 0)
    {
        M_massMatrix.reset(new matrix_type(M_fespace->map()));
        M_adrAssembler.addMass(M_massMatrix,1.0);
        M_massMatrix->globalAssemble();
    }

    M_systemMatrix.reset(new matrix_type(M_fespace->map()));
    M_rhsMatrix.reset(new matrix_type(M_fespace->map()));

    // Erase the old matrices
    *M_systemMatrix *= 0.0;
    *M_rhsMatrix *=0.0;


    // Mass matrix
    *M_systemMatrix += (*M_massMatrix)*( M_bdf->coefficientFirstDerivative(0)/M_data->dataTime()->timeStep() );

    // advection matrix
    M_adrAssembler.addAdvection(M_systemMatrix,beta);

    // ip stabilization
    if (M_data->stabilization() == DataLevelSet::IP)
    {
        if (M_data->IPTreatment() == DataLevelSet::IMPLICIT)
        {
            M_ipAssembler.addIPStabilization(M_systemMatrix,beta,M_data->IPCoef());
        }
        else if (M_data->IPTreatment() == DataLevelSet::SEMI_IMPLICIT)
        {
            M_ipAssembler.addIPStabilizationStencil(M_systemMatrix,M_rhsMatrix,beta,M_data->IPCoef());
        }
        else if (M_data->IPTreatment() == DataLevelSet::EXPLICIT)
        {
            M_ipAssembler.addIPStabilization(M_rhsMatrix,beta,M_data->IPCoef());
        }
    }

    // Global Assembly

    M_systemMatrix->globalAssemble();
    M_rhsMatrix->globalAssemble();


    // Erase old rhs
    M_rhs *=0.0;

    // Rhs time
     M_bdf->updateRHSContribution(M_data->dataTime()->timeStep());
     M_rhs += *M_massMatrix *M_bdf->rhsContributionFirstDerivative();

    // Rhs stab
    M_rhs += *M_rhsMatrix*M_solution;


    // Manage the BCs
    bcHandler.bcUpdate(*M_fespace->mesh(),M_fespace->feBd(),M_fespace->dof());
    bcManage(*M_systemMatrix,M_rhs,*M_fespace->mesh(),M_fespace->dof(),bcHandler,M_fespace->feBd(),1.0,time);
}


template< typename mesh_type, typename solver_type>
void
LevelSetSolver<mesh_type,solver_type>::
setupLinearSolver(const GetPot& dataFile, const string& section)
{
    M_linearSolver.setDataFromGetPot( dataFile, (section+"/solver").data() );
    M_linearSolver.setupPreconditioner( dataFile, (section+"/prec").data() );
}

template< typename mesh_type, typename solver_type>
void
LevelSetSolver<mesh_type,solver_type>::
iterate()
{
    ASSERT(M_systemMatrix!= 0, "No system matrix for the linear system! Build it before.");
    ASSERT(M_bdf!=0, "No bdf structure. Use the setup method before the iterate.");

    M_linearSolver.setMatrix(*M_systemMatrix);

    M_linearSolver.solveSystem(M_rhs, M_solution, M_systemMatrix);

    M_bdf->shiftRight(M_solution);
}



template<typename mesh_type, typename solver_type>
void
LevelSetSolver<mesh_type,solver_type>::
reinitializationDirect()
{
    updateFacesNormalsRadius();


    // Initialization of the MPI things
    int nb_proc(0);
    int my_rank(0);

    const Epetra_MpiComm* my_comm = dynamic_cast<Epetra_MpiComm const*>(&(M_fespace->map().comm()));

    unsigned int dim(3);
    int nPt (M_fespace->mesh()->storedPoints());

    nb_proc = my_comm->NumProc();
    my_rank = my_comm->MyPID();
    MPI_Status status;

    std::vector< std::vector< Real > > my_points(nPt, std::vector<Real>(dim));
    for (int iter_pt(0); iter_pt < nPt; ++iter_pt)
    {
        my_points[iter_pt][0]=M_fespace->mesh()->point(iter_pt).x();
        my_points[iter_pt][1]=M_fespace->mesh()->point(iter_pt).y();
        my_points[iter_pt][2]=M_fespace->mesh()->point(iter_pt).z();
    };



    // Storage for the distance
    std::vector < std::vector < Real > > temp_distances(nb_proc, vector< Real>(nPt,0.0));

    // Communications

    my_comm->Barrier();

    // Start the communications

    for (int i(0); i<nb_proc; ++i)
    {
        // Sender
        if (my_rank==i)
        {
            //std::cout << " Reinitialize part " << my_rank << std::endl;

            for (int j(0); j<nb_proc; ++j)
            {
                if (j!=i) MPI_Send(&nPt,1,MPI_INT,j,1,my_comm->Comm());
            };

            //std::cout << " Send points " << std::endl;
            for (int d(0); d<nPt; ++d)
            {
                for (int j(0); j<nb_proc; ++j)
                {
                    if (j!=i)
                    {
                        //MPI_Request req;
                        //MPI_Isend(&my_points[d][0],dim,MPI_DOUBLE,j,1,my_comm->Comm(),&req);
                        //MPI_Wait(&req,&status);
                        MPI_Ssend(&my_points[d][0],dim,MPI_DOUBLE,j,1,my_comm->Comm());
                    };
                };
            };

            //std::cout << " Compute " << std::endl;
            for (int iter_pt(0); iter_pt<nPt; ++iter_pt)
            {
                vector<Real> dof_coord;
                dof_coord.push_back(my_points[iter_pt][0]);
                dof_coord.push_back(my_points[iter_pt][1]);
                dof_coord.push_back(my_points[iter_pt][2]);
                temp_distances[i][iter_pt]=computeUnsignedDistance(dof_coord);
            };

            //std::cout << " Get the distances " << std::endl;
            for (int distToRecieve(0); distToRecieve<nPt; ++distToRecieve)
            {
                for (int j(0); j<nb_proc; ++j)
                {
                    if (j!=i)
                    {
                        //MPI_Request req;
                        //MPI_Irecv(&(temp_distances[j][distToRecieve]),1,MPI_DOUBLE,j,1,my_comm->Comm(),&req);
                        //MPI_Wait(&req,&status);
                        MPI_Recv(&(temp_distances[j][distToRecieve]),1,MPI_DOUBLE,j,1,my_comm->Comm(),&status);
                    }
                };
            };


        }
        else
        {
            // Reciever

            int n_pt_to_recieve(0);
            MPI_Recv(&n_pt_to_recieve,1,MPI_INT,i,1,my_comm->Comm(),&status);
            //std::cout << " Points (R) : " << n_pt_to_recieve << std::endl;

            std::vector< std::vector< Real > > current_points(n_pt_to_recieve, std::vector< Real> (dim,0));
            std::vector< Real> current_distances(n_pt_to_recieve,0);

            for (int d(0); d< n_pt_to_recieve; ++d)
            {
                //MPI_Request req;
                //MPI_Irecv(&current_points[d][0],dim,MPI_DOUBLE,i,1,my_comm->Comm(),&req);
                //MPI_Wait(&req,&status);
                MPI_Recv(&current_points[d][0],dim,MPI_DOUBLE,i,1,my_comm->Comm(),&status);
            };

            for (int iter_pt(0); iter_pt<n_pt_to_recieve; ++iter_pt)
            {
                std::vector< Real > current_pt(dim,0);
                current_pt[0]=current_points[iter_pt][0];
                current_pt[1]=current_points[iter_pt][1];
                current_pt[2]=current_points[iter_pt][2];

                current_distances[iter_pt]=computeUnsignedDistance(current_pt);
            };

            for (int iter_pt(0); iter_pt<n_pt_to_recieve; ++iter_pt)
            {
                //MPI_Request req;
                //MPI_Isend(&current_distances[iter_pt],1,MPI_DOUBLE,i,1,my_comm->Comm(),&req);
                //MPI_Wait(&req,&status);
                MPI_Ssend(&current_distances[iter_pt],1,MPI_DOUBLE,i,1,my_comm->Comm());
            };

        }; // end if

        my_comm->Barrier();

    }; // end for i


    // Post processing
    // Also give the sign!

    vector_type repSol(M_solution,Repeated);

    for (int iter_pt(0); iter_pt<nPt; ++iter_pt)
    {
        // Compute the minimum over the points
        Real abs_dist(temp_distances[0][iter_pt]);

        for (int j(1); j<nb_proc; ++j)
        {
            if (temp_distances[j][iter_pt] < abs_dist)
            {
                abs_dist = temp_distances[j][iter_pt];
            };
        };

        ID my_id(M_fespace->mesh()->point(iter_pt).id());
        int sign(1);
        if (repSol(my_id) < 0) {sign = -1;};
        repSol(my_id) = abs_dist*sign;
    };

    my_comm->Barrier();

    M_solution = vector_type(repSol,Unique,Zero);

}

// ===================================================
// Private Methods
// ===================================================

template<typename mesh_type, typename solver_type>
void
LevelSetSolver<mesh_type,solver_type>::
updateFacesNormalsRadius()
{
    // Use a repeated vector, as we need the boundary
    // values on both processor domains
    vector_type rep_solution(M_solution, Repeated);

    cleanFacesData();

    UInt n_el(M_fespace->mesh()->numElements());
    UInt n_dof(M_fespace->dof().numLocalDof());
    for (ID iter_el(0); iter_el < n_el; ++iter_el)
    {

        // Storage for the points on the face
        std::vector< point_type> intersection_points;


        // Find the edges that are crossed by the interface
        for (UInt iter_dof_1(0); iter_dof_1 < n_dof; ++iter_dof_1)
        {
            ID id_dof_1( M_fespace->dof().localToGlobalMap(iter_el, iter_dof_1));
            Real val1 (rep_solution(id_dof_1));

            for (UInt iter_dof_2(iter_dof_1+1); iter_dof_2 < n_dof ; ++iter_dof_2)
            {
                ID id_dof_2( M_fespace->dof().localToGlobalMap(iter_el, iter_dof_2));
                Real val2 (rep_solution(id_dof_2));

                // Check if there is a change in the sign
                // DO NOT TAKE INTO ACCOUNT THE DEGENERATED EDGES
                if ((val1*val2 <=0) && ((val1!=0) || (val2!=0)))
                {
                    Real lambda(-val1/(val2-val1));

                    point_type inter(3,0);

                    inter[0]= (1-lambda)* M_fespace->mesh()->element(iter_el).point(iter_dof_1).x()
                              + lambda* M_fespace->mesh()->element(iter_el).point(iter_dof_2).x() ;
                    inter[1]= (1-lambda)* M_fespace->mesh()->element(iter_el).point(iter_dof_1).y()
                              + lambda* M_fespace->mesh()->element(iter_el).point(iter_dof_2).y() ;
                    inter[2]= (1-lambda)* M_fespace->mesh()->element(iter_el).point(iter_dof_1).z()
                              + lambda* M_fespace->mesh()->element(iter_el).point(iter_dof_2).z() ;

                    intersection_points.push_back(inter);

                };
            };
        };

        // Now we know the points on the faces
        // We want to store them and compute
        // the radius for each of them

        if (intersection_points.size() ==3)
        {
            // Store the face

            face_type my_face;
            my_face.push_back(intersection_points[0]);
            my_face.push_back(intersection_points[1]);
            my_face.push_back(intersection_points[2]);

            M_faces.push_back(my_face);

            // Compute the radius

            Real radius(0);
            Real current_dist(0);
            current_dist=distanceBetweenPoints(intersection_points[0],intersection_points[1]);
            if (current_dist > radius) radius = current_dist;
            current_dist=distanceBetweenPoints(intersection_points[0],intersection_points[2]);
            if (current_dist > radius) radius = current_dist;
            current_dist=distanceBetweenPoints(intersection_points[1],intersection_points[2]);
            if (current_dist > radius) radius = current_dist;

            M_radius.push_back(radius);

            // Compute the normal

            point_type v1(3,0);
            v1[0]=intersection_points[1][0]-intersection_points[0][0];
            v1[1]=intersection_points[1][1]-intersection_points[0][1];
            v1[2]=intersection_points[1][2]-intersection_points[0][2];

            point_type v2(3,0);
            v2[0]=intersection_points[2][0]-intersection_points[0][0];
            v2[1]=intersection_points[2][1]-intersection_points[0][1];
            v2[2]=intersection_points[2][2]-intersection_points[0][2];

            point_type normal_face(crossProd(v1,v2));
            Real my_norm(norm(normal_face));
            normal_face[0]=normal_face[0]/my_norm;
            normal_face[1]=normal_face[1]/my_norm;
            normal_face[2]=normal_face[2]/my_norm;

            M_normals.push_back(normal_face);


        }
        else if (intersection_points.size() ==4)
        {
            // Here we have to beware to choose the right one!
            // Store the face

            face_type my_face;
            my_face.push_back(intersection_points[0]);
            my_face.push_back(intersection_points[1]);
            my_face.push_back(intersection_points[2]);

            M_faces.push_back(my_face);

            // Compute the radius

            Real radius(0);
            Real current_dist(0);
            current_dist=distanceBetweenPoints(intersection_points[0],intersection_points[1]);
            if (current_dist > radius) radius = current_dist;
            current_dist=distanceBetweenPoints(intersection_points[0],intersection_points[2]);
            if (current_dist > radius) radius = current_dist;
            current_dist=distanceBetweenPoints(intersection_points[1],intersection_points[2]);
            if (current_dist > radius) radius = current_dist;

            M_radius.push_back(radius);

            point_type v1(3,0);
            v1[0]=intersection_points[1][0]-intersection_points[0][0];
            v1[1]=intersection_points[1][1]-intersection_points[0][1];
            v1[2]=intersection_points[1][2]-intersection_points[0][2];

            point_type v2(3,0);
            v2[0]=intersection_points[2][0]-intersection_points[0][0];
            v2[1]=intersection_points[2][1]-intersection_points[0][1];
            v2[2]=intersection_points[2][2]-intersection_points[0][2];

            point_type normal_face(crossProd(v1,v2));
            Real my_norm(norm(normal_face));
            normal_face[0]=normal_face[0]/my_norm;
            normal_face[1]=normal_face[1]/my_norm;
            normal_face[2]=normal_face[2]/my_norm;

            M_normals.push_back(normal_face);


            // Find the second face
            face_type my_second_face;

            my_second_face.push_back(intersection_points[3]);

            for (int i(0); i<3; ++i)
            {
                int other_1((i+1)%3);
                int other_2((i+2)%3);

                point_type v(3,0);
                v[0]=intersection_points[other_2][0]-intersection_points[other_1][0];
                v[1]=intersection_points[other_2][1]-intersection_points[other_1][1];
                v[2]=intersection_points[other_2][2]-intersection_points[other_1][2];

                point_type normal_v(crossProd(v,normal_face));

                point_type o1_to_i(3,0);
                o1_to_i[0]=intersection_points[i][0]-intersection_points[other_1][0];
                o1_to_i[1]=intersection_points[i][1]-intersection_points[other_1][1];
                o1_to_i[2]=intersection_points[i][2]-intersection_points[other_1][2];

                point_type o1_to_inter3(3,0);
                o1_to_inter3[0]=intersection_points[3][0]-intersection_points[other_1][0];
                o1_to_inter3[1]=intersection_points[3][1]-intersection_points[other_1][1];
                o1_to_inter3[2]=intersection_points[3][2]-intersection_points[other_1][2];

                if ( scalProd(normal_v,o1_to_i) * scalProd(normal_v,o1_to_inter3) > 0)
                    my_second_face.push_back(intersection_points[i]);

            };

            ASSERT(my_second_face.size()==3,"Internal error while computing the faces");
            M_faces.push_back(my_second_face);

            // Compute its radius

            current_dist=0;
            current_dist=distanceBetweenPoints(my_second_face[0],my_second_face[1]);
            if (current_dist > radius) radius = current_dist;
            current_dist=distanceBetweenPoints(my_second_face[0],my_second_face[2]);
            if (current_dist > radius) radius = current_dist;
            current_dist=distanceBetweenPoints(my_second_face[1],my_second_face[2]);
            if (current_dist > radius) radius = current_dist;

            M_radius.push_back(radius);

            // Normal is the same
            // (already supposed in the computations)
            M_normals.push_back(normal_face);



        }
        else if (intersection_points.size()>0)
        {
            std::cout << intersection_points.size() << std::endl;
            ERROR_MSG("Internal error in the localization of the interface");
        };
    }; // end for iter_el

    ASSERT(M_faces.size() == M_radius.size(),"Internal non coerance");
    ASSERT(M_faces.size() == M_normals.size(),"Internal non coerance");

}

template<typename mesh_type, typename solver_type>
Real
LevelSetSolver<mesh_type,solver_type>::
computeUnsignedDistance(const std::vector< Real >& point)
{
    // We use here a 2 steps algorithm.
    // The first step is used to erase quickly all the
    //  the faces that are not good candidates
    // The second step then compute the exact distance
    //  for the remaining faces.

    // Initialisation
    Real absdist(std::numeric_limits<double>::max());

    // Informations
    int nbFaces (this->M_faces.size());

    // Storage of the distances with the center of gravity
    std::vector<Real> dist_lower_bound(nbFaces);


    // ==> STEP 1
    // A first run to get a rough estimations and lower bounds
    // We just test the first point as this is fast and already
    // eliminates lot of faces

    for (int iter_face(0); iter_face < nbFaces; ++iter_face)
    {
        Real current_dist(distanceBetweenPoints(point,M_faces[iter_face][0]));
        dist_lower_bound[iter_face] =  current_dist - M_radius[iter_face];

        if (current_dist < absdist) absdist = current_dist;
    };


    // ==> STEP 2
    // Compute the distance using only the candidates
    // selected by the lower bound
    //
    for (int iter_face(0); iter_face < nbFaces; ++iter_face)
    {
        // Exclude cells that are too far
        if (dist_lower_bound[iter_face] <= absdist)
        {

            Real current_abs_dist(std::numeric_limits<double>::max());


            // First of all, project the point on the plan
            //  of the face

            point_type point_to_face(3,0);
            point_to_face[0]=point[0] - M_faces[iter_face][0][0];
            point_to_face[1]=point[1] - M_faces[iter_face][0][1];
            point_to_face[2]=point[2] - M_faces[iter_face][0][2];

            Real lambda = scalProd(point_to_face,M_normals[iter_face]);
            point_type projected_point(3,0);
            projected_point[0] = point[0]-lambda*M_normals[iter_face][0];
            projected_point[1] = point[1]-lambda*M_normals[iter_face][1];
            projected_point[2] = point[2]-lambda*M_normals[iter_face][2];

            // Compute the distance
            // For debugging purpose, this should be equal to std::abs(lambda)
            current_abs_dist = distanceBetweenPoints(point,projected_point);

            // This is the distance only if the projected point
            //  belongs to the present facel. We test this.

            for (int iter_edge(0); iter_edge<3 ; ++iter_edge)
            {
                int vertex1(iter_edge), vertex2((iter_edge+1)%3), vertex_out((iter_edge+2)%3);

                // Check if the projected point is on the right side of the edge
                point_type edge_vector(3,0);
                edge_vector[0]=M_faces[iter_face][vertex2][0]-M_faces[iter_face][vertex1][0];
                edge_vector[1]=M_faces[iter_face][vertex2][1]-M_faces[iter_face][vertex1][1];
                edge_vector[2]=M_faces[iter_face][vertex2][2]-M_faces[iter_face][vertex1][2];

                // Compute the "normal to the edge", the vector that is orthogonal to
                // the edge and to the normal of the face
                point_type edge_normal(crossProd(M_normals[iter_face],edge_vector));

                Real edge_normal_norm(norm(edge_normal));
                edge_normal[0]=edge_normal[0]/edge_normal_norm;
                edge_normal[1]=edge_normal[1]/edge_normal_norm;
                edge_normal[2]=edge_normal[2]/edge_normal_norm;


                // Compute the vectors that link point and vertex_out to the edge
                point_type projpoint_to_edge(3,0);
                projpoint_to_edge[0] = projected_point[0] - M_faces[iter_face][vertex1][0];
                projpoint_to_edge[1] = projected_point[1] - M_faces[iter_face][vertex1][1];
                projpoint_to_edge[2] = projected_point[2] - M_faces[iter_face][vertex1][2];

                point_type vertex_out_to_edge(3,0);
                vertex_out_to_edge[0] = M_faces[iter_face][vertex_out][0] - M_faces[iter_face][vertex1][0];
                vertex_out_to_edge[1] = M_faces[iter_face][vertex_out][1] - M_faces[iter_face][vertex1][1];
                vertex_out_to_edge[2] = M_faces[iter_face][vertex_out][2] - M_faces[iter_face][vertex1][2];

                // Check whether the point and the vertex_out are on the same side
                // of the edge
                if (scalProd(edge_normal,projpoint_to_edge) * scalProd(edge_normal,vertex_out_to_edge) < 0)
                {
                    // The projected point is on the wrong side of the edge, so we project it
                    //  again, but this time on the edge (inside the face plan)
                    // The plan of the projection is determined by the vector of the edge
                    //  and the normal of the face, what means that the normal of the edge
                    //  is already the projection direction.

                    Real lambda_2(scalProd(edge_normal,projpoint_to_edge));
                    point_type projproj_point(3,0);
                    projproj_point[0] = projected_point[0] - lambda_2*edge_normal[0];
                    projproj_point[1] = projected_point[1] - lambda_2*edge_normal[1];
                    projproj_point[2] = projected_point[2] - lambda_2*edge_normal[2];

                    point_type projprojpoint_to_edge(3,0);
                    projprojpoint_to_edge[0] = projproj_point[0] - M_faces[iter_face][vertex1][0];
                    projprojpoint_to_edge[1] = projproj_point[1] - M_faces[iter_face][vertex1][1];
                    projprojpoint_to_edge[2] = projproj_point[2] - M_faces[iter_face][vertex1][2];

                    Real report(scalProd(projprojpoint_to_edge,edge_vector)/scalProd(edge_vector,edge_vector));

                    if (report<0)
                    {
                        current_abs_dist = distanceBetweenPoints(point,M_faces[iter_face][vertex1]);
                    }
                    else if (report<=1)
                    {
                        current_abs_dist = distanceBetweenPoints(point,projproj_point);
                    }
                    else
                    {
                        current_abs_dist = distanceBetweenPoints(point,M_faces[iter_face][vertex2]);
                    };

                    // To avoid useless computations, we
                    //  break the "for iter_edge" loop
                    break;
                };

            };

            // Update if needed

            if (current_abs_dist <= absdist)
            {
                absdist = current_abs_dist;
            };
        };
    };

    return absdist;
}

template<typename mesh_type, typename solver_type>
void
LevelSetSolver<mesh_type,solver_type>::
cleanFacesData()
{
    M_faces.clear();
    M_normals.clear();
    M_radius.clear();
}

} // Namespace LifeV

#endif /* LEVELSETSOLVER_H */
