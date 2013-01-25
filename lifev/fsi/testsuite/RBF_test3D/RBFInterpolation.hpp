#ifndef RBF_INTERPOLATION_HPP
#define RBF_INTERPOLATION_HPP

#include <lifev/core/LifeV.hpp>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include "EpetraExt_CrsMatrixIn.h"
#include <lifev/core/array/MatrixEpetra.hpp>
#include "AztecOO.h"
#include <lifev/core/array/GhostHandler.hpp>

namespace LifeV 
{
  template <typename Mesh>
  class RBFInterpolation
  {
    
  public:
    
    typedef Mesh                                              mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                      meshPtr_Type;
    typedef VectorEpetra                                      vector_Type;
    typedef boost::shared_ptr<vector_Type >                   vectorPtr_Type;
    typedef Epetra_CrsMatrix                                  matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                    matrixPtr_Type;
    typedef std::vector<int>                                  flagContainer_Type;
    typedef std::set<ID>                                      idContainer_Type;
    typedef Epetra_Vector                                     vectorEpetra_Type;
    typedef boost::shared_ptr<Epetra_Vector>                  vectorEpetraPtr_Type;
    typedef boost::shared_ptr<Epetra_Map>                     mapEpetraPtr_type;
    typedef GhostHandler<mesh_Type>                           neighbors_Type;
    typedef boost::shared_ptr<neighbors_Type>                 neighborsPtr_Type;

    //! Constructor
    RBFInterpolation( meshPtr_Type fullMeshKnown, 
		      meshPtr_Type localMeshKnown, 
		      meshPtr_Type fullMeshUnknown, 
		      meshPtr_Type localMeshUnknown, 
		      flagContainer_Type flags);
    
    //! Destructor
    ~RBFInterpolation(){}

    //! Setup the RBF data
    void setupRBFData(vectorPtr_Type KnownField, vectorPtr_Type UnknownField /*, double radius*/);

    //! Set the value of the RBF radius 
    void setRBFradius(double const & radius ){ M_RBFRadius = radius; }

    //! Build the RBF Operators, namely the interpolant and the projection ones.
    void buildOperators();
    
    //! Build the RBF interpolant
    void InterpolationOperator();

    //! Identifies nodes with an assigned markerID 
    void identifyNodes(meshPtr_Type LocalMesh, std::set<ID> & GID_nodes, vectorPtr_Type CheckVector);

     //! Check if the point with markerID pointMarker has to be selected
    bool isInside(ID pointMarker, flagContainer_Type Flags);

    //! Evaluation of the RBF
    double rbf(double x1, double y1, double z1, double x2, double y2, double z2, double radius);

    //! Evaluate the RBF radius
    double computeRBFradius(meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID);

    //! Build the projection operator
    void ProjectionOperator();

     //! Spy of an EpetraCrsMatrix
    void spy( std::string const &fileName, matrixPtr_Type A);

    //! Prepare Rhs
    void buildRhs();

    //! Manages the solution of the interpolation problem
    void interpolate();

    //! It solves the interpolation problem, namely it solves a system and performs a matrix-vector multiplication
    void solve(matrixPtr_Type A, vectorEpetra_Type & x, vectorEpetraPtr_Type b);

    //! Getter for the solution
    void solution(vectorPtr_Type & Solution);

    //! Getter for the solution obtained by a pure RBF approach
    void solutionrbf(vectorPtr_Type & Solution_rbf);

  private:
    
    meshPtr_Type          M_fullMeshKnown;
    meshPtr_Type          M_localMeshKnown;
    meshPtr_Type          M_fullMeshUnknown;
    meshPtr_Type          M_localMeshUnknown;
    matrixPtr_Type        M_interpolationOperator;
    matrixPtr_Type        M_projectionOperator;
    flagContainer_Type    M_flags;
    double                M_RBFRadius;
    vectorPtr_Type        M_knownField;
    vectorPtr_Type        M_unknownField;
    idContainer_Type      M_GIdsKnownMesh;
    idContainer_Type      M_GIdsUnknownMesh;
    vectorEpetraPtr_Type  M_RhsF;
    vectorEpetraPtr_Type  M_RhsOne;
    mapEpetraPtr_type     M_interpolationOperatorMap;
    mapEpetraPtr_type     M_projectionOperatorMap;
    neighborsPtr_Type     M_neighbors;
    vectorPtr_Type        M_unknownField_rbf;
  };

  template <typename Mesh>
  RBFInterpolation<Mesh>::RBFInterpolation( meshPtr_Type fullMeshKnown, 
					    meshPtr_Type localMeshKnown, 
					    meshPtr_Type fullMeshUnknown, 
					    meshPtr_Type localMeshUnknown, 
					    flagContainer_Type flags):
    M_fullMeshKnown( fullMeshKnown ),
    M_localMeshKnown( localMeshKnown ),
    M_fullMeshUnknown( fullMeshUnknown ),
    M_localMeshUnknown( localMeshUnknown ),
    M_flags( flags )
  {
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::setupRBFData(vectorPtr_Type KnownField, vectorPtr_Type UnknownField )
  {
    M_knownField   = KnownField;
    M_unknownField = UnknownField;
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::buildOperators()
  {
    this->InterpolationOperator();
    this->ProjectionOperator();
    this->buildRhs();

    spy("PHII", M_interpolationOperator);
    spy("PHII_f", M_projectionOperator);

  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::InterpolationOperator()
  {
    this->identifyNodes(M_localMeshKnown, M_GIdsKnownMesh, M_knownField);
    M_neighbors.reset( new neighbors_Type( M_fullMeshKnown, M_localMeshKnown, M_knownField->mapPtr(), M_knownField->mapPtr()->commPtr() ) );
    M_neighbors->mysetUp();
    
    int LocalNodesNumber = M_GIdsKnownMesh.size();
    int TotalNodesNumber = 0;

    MPI_Allreduce(&LocalNodesNumber, &TotalNodesNumber, 1,  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    vector<double>   RBF_radius(LocalNodesNumber);
    vector<set<ID> > MatrixGraph(LocalNodesNumber);
    int * ElementsPerRow = new int[LocalNodesNumber];
    int * GlobalID = new int[LocalNodesNumber];
    int k = 0;
    int Max_entries = 0;
    
    for(set<ID>::iterator it = M_GIdsKnownMesh.begin(); it != M_GIdsKnownMesh.end(); ++it)
      {
	GlobalID[k] = *it;
	MatrixGraph[k] = M_neighbors->nodeNodeNeighborsList()[GlobalID[k]];
	MatrixGraph[k].insert(GlobalID[k]);
	RBF_radius[k] = computeRBFradius( M_fullMeshKnown, M_fullMeshKnown, MatrixGraph[k], GlobalID[k]);
	ElementsPerRow[k] = MatrixGraph[k].size();
	if(ElementsPerRow[k]>Max_entries)
	  Max_entries = ElementsPerRow[k];
	++k;
      }
    
    M_interpolationOperatorMap.reset(new Epetra_Map(TotalNodesNumber, LocalNodesNumber, GlobalID, 0, *(M_knownField->mapPtr()->commPtr())));
    M_interpolationOperator.reset(new matrix_Type(Copy, *M_interpolationOperatorMap, ElementsPerRow));
    
    int * Indices = new int[Max_entries];
    double * Values = new double[Max_entries];
    
    for( int i = 0 ; i < LocalNodesNumber; ++i )
      {
	k = 0;
	for( set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
	  {
	      Indices[k] = *it;
	      Values[k]  = rbf( M_fullMeshKnown->point(GlobalID[i]).x(),
				M_fullMeshKnown->point(GlobalID[i]).y(),
				M_fullMeshKnown->point(GlobalID[i]).z(),
				M_fullMeshKnown->point(*it).x(),
				M_fullMeshKnown->point(*it).y(),
				M_fullMeshKnown->point(*it).z(),
				RBF_radius[i]);
	      ++k;
	    }
	M_interpolationOperator->InsertGlobalValues(GlobalID[i], k, Values, Indices);
      }      
    M_interpolationOperator->FillComplete();
  
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::ProjectionOperator()
  {
    
    this->identifyNodes(M_localMeshUnknown, M_GIdsUnknownMesh, M_unknownField);
    
    int LocalNodesNumber = M_GIdsUnknownMesh.size();
    int TotalNodesNumber = 0;

    vector<double>   RBF_radius(LocalNodesNumber);
    vector<set<ID> > MatrixGraph(LocalNodesNumber);
    int * ElementsPerRow = new int[LocalNodesNumber];
    int * GlobalID = new int[LocalNodesNumber];
    int k = 0;
    int Max_entries = 0;
    double d;
    double d_min;
    int nearestPoint;
    
    for(set<ID>::iterator it = M_GIdsUnknownMesh.begin(); it != M_GIdsUnknownMesh.end(); ++it)
      {
	GlobalID[k] = *it;
	d_min = 100;
	for (int j = 0; j <  M_fullMeshKnown->numVertices(); ++j)
	  {
	    if( M_flags[0] == -1 || this->isInside(M_fullMeshKnown->point(j).markerID(), M_flags) )
	      {
		d = std::sqrt( pow(M_fullMeshKnown->point(j).x()-M_fullMeshUnknown->point(GlobalID[k]).x(),2)
			       + pow(M_fullMeshKnown->point(j).y()-M_fullMeshUnknown->point(GlobalID[k]).y(),2)
			       + pow(M_fullMeshKnown->point(j).z()-M_fullMeshUnknown->point(GlobalID[k]).z(),2) );
		if (d < d_min)
		  {
		    d_min = d;
		    nearestPoint = M_fullMeshKnown->point(j).id();
		  }
	      }
	  }
	MatrixGraph[k] = M_neighbors->nodeNodeNeighborsList()[nearestPoint];
	MatrixGraph[k].insert(nearestPoint);       
        RBF_radius[k] = computeRBFradius( M_fullMeshKnown, M_fullMeshUnknown, MatrixGraph[k], GlobalID[k]);	     
	ElementsPerRow[k] = MatrixGraph[k].size();
	if(ElementsPerRow[k] > Max_entries)
	  Max_entries = ElementsPerRow[k];
	++k;
      }
	
    MPI_Allreduce(&LocalNodesNumber, &TotalNodesNumber, 1,  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    M_projectionOperatorMap.reset(new Epetra_Map(TotalNodesNumber, LocalNodesNumber, GlobalID, 0, *(M_unknownField->mapPtr()->commPtr())));
    M_projectionOperator.reset(new matrix_Type(Copy, *M_projectionOperatorMap, ElementsPerRow));
    
    int * Indices = new int[Max_entries];
    double * Values = new double[Max_entries];

    for( int i = 0 ; i < LocalNodesNumber; ++i )
      {
        k = 0;
        for( set<ID>::iterator it = MatrixGraph[i].begin(); it != MatrixGraph[i].end(); ++it)
	  {
            Indices[k] = *it;
            Values[k]  = rbf( M_fullMeshUnknown->point(GlobalID[i]).x(),
			      M_fullMeshUnknown->point(GlobalID[i]).y(),
			      M_fullMeshUnknown->point(GlobalID[i]).z(),
			      M_fullMeshKnown->point(*it).x(),
			      M_fullMeshKnown->point(*it).y(),
			      M_fullMeshKnown->point(*it).z(),
			      RBF_radius[i]);
            ++k;
	  }
        M_projectionOperator->InsertGlobalValues(GlobalID[i], k, Values, Indices);
      }
    M_projectionOperator->FillComplete(*M_interpolationOperatorMap, *M_projectionOperatorMap); 
    spy("PHI_f", M_projectionOperator);
    
    delete Indices;
    delete Values;
    delete ElementsPerRow;
    delete GlobalID;

  }

  template <typename Mesh>
  double RBFInterpolation<Mesh>::computeRBFradius(meshPtr_Type MeshNeighbors, meshPtr_Type MeshGID, idContainer_Type Neighbors, ID GlobalID)
  {
    double r = 0;
    double r_max = 0;
    for(idContainer_Type::iterator it = Neighbors.begin(); it != Neighbors.end(); ++it)
	  {
	    r = std::sqrt( pow( MeshGID->point( GlobalID ).x() - MeshNeighbors->point( *it ).x(), 2 )
			   + pow( MeshGID->point( GlobalID ).y() - MeshNeighbors->point( *it ).y(), 2 )
			   + pow( MeshGID->point( GlobalID ).z() - MeshNeighbors->point( *it ).z(), 2 ) );
	    r_max = ( r > r_max ) ? r : r_max;
	  }
    return r_max;
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::buildRhs()
  {
    int k = 0;
    M_RhsF.reset(new vectorEpetra_Type(M_interpolationOperator->RowMap()));
    M_RhsOne.reset(new vectorEpetra_Type(M_interpolationOperator->RowMap()));
    for(set<ID>::iterator it = M_GIdsKnownMesh.begin(); it != M_GIdsKnownMesh.end(); ++it) 
      if(M_knownField->blockMap().LID(*it)!=-1)
	{
	  (*M_RhsF)[k] = (*M_knownField)[*it];
	  (*M_RhsOne)[k] = 1;
	  ++k;
	}
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::interpolate()
  {
    vectorEpetra_Type gamma_f(*M_interpolationOperatorMap);
    vectorEpetra_Type gamma_one(*M_interpolationOperatorMap);

    this->solve(M_interpolationOperator, gamma_f, M_RhsF);
    this->solve(M_interpolationOperator, gamma_one, M_RhsOne);

    vectorEpetra_Type rbf_f(*M_projectionOperatorMap);
    vectorEpetra_Type rbf_one(*M_projectionOperatorMap);
    vectorEpetra_Type solution(*M_projectionOperatorMap);

    M_projectionOperator->Multiply(false, gamma_f, rbf_f);	  
    M_projectionOperator->Multiply(false, gamma_one, rbf_one);

    for ( UInt i = 0; i < rbf_f.MyLength(); ++i )
      solution[i] = rbf_f[i]/rbf_one[i];

    M_unknownField_rbf.reset(new vector_Type(M_unknownField->map()));

    int k = 0;
    for(set<ID>::iterator it = M_GIdsUnknownMesh.begin(); it != M_GIdsUnknownMesh.end(); ++it)
      {
	(*M_unknownField_rbf)[*it] = rbf_f[k];
	(*M_unknownField)[*it] = solution[k];
	++k;
      }
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::solve(matrixPtr_Type A, vectorEpetra_Type & x, vectorEpetraPtr_Type b)
  {
    Epetra_LinearProblem Problem(&(*A), &x, &(*b));
    AztecOO solver (Problem);
    solver.SetAztecOption(AZ_solver, AZ_gmres);
    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    solver.SetAztecOption(AZ_output, AZ_last);
    int Niters = 1500;
    double tol = 1e-10;
    solver.Iterate(Niters, tol);
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::identifyNodes(meshPtr_Type LocalMesh, std::set<ID> & GID_nodes, vectorPtr_Type CheckVector)
  {
    
    if(M_flags[0] == -1)
      {
	for ( UInt i = 0; i < LocalMesh->numVertices(); ++i )
	  if(CheckVector->blockMap().LID(LocalMesh->point(i).id()) != -1)
	    GID_nodes.insert(LocalMesh->point(i).id()); 
      }
    else
      {
	for ( UInt i = 0; i < LocalMesh->numVertices(); ++i )
	  if( this->isInside(LocalMesh->point(i).markerID(), M_flags) )
	    if(CheckVector->blockMap().LID(LocalMesh->point(i).id()) != -1)
	      GID_nodes.insert(LocalMesh->point(i).id()); 
      }
    
  }
  
  template <typename Mesh>
  bool RBFInterpolation<Mesh>::isInside(ID pointMarker, flagContainer_Type flags)
  {
    int check = 0;
    for(UInt i = 0; i < flags.size(); ++i)
      if(pointMarker == flags[i])
	++check;
    return (check > 0) ? true : false;
  }

  template <typename Mesh>
  double RBFInterpolation<Mesh>::rbf(double x1, double y1, double z1, double x2, double y2, double z2, double radius)
  {
    double distance = sqrt( pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2));
    return pow(1-distance/radius,4)*(4*distance/radius+1); 
  }
  
  template <typename Mesh>
  void RBFInterpolation<Mesh>::solution(vectorPtr_Type & Solution)
  {
    Solution = M_unknownField;
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::solutionrbf(vectorPtr_Type & Solution_rbf)
  {
    Solution_rbf = M_unknownField_rbf;
  }

  template <typename Mesh>
  void RBFInterpolation<Mesh>::spy( std::string const &fileName, matrixPtr_Type A)
  {
    std::string name = fileName, uti = " , ";
    Int  me = A->Comm().MyPID();
    std::ostringstream myStream;
    myStream << me;
    name = fileName + ".m";
    EpetraExt::RowMatrixToMatlabFile( name.c_str(), *A );
  }

} // namespace LifeV

#endif // RBF_INTERPOLATION_HPP
