#ifndef DATA_INTERPOLATION_HPP
#define DATA_INTERPOLATION_HPP

#include <lifev/core/LifeV.hpp>
#include <Epetra_Vector.h>
#include "Epetra_Map.h"
#include <Epetra_CrsMatrix.h>
#include <lifev/core/array/GhostHandler.hpp>

namespace LifeV 
{
  template <typename Mesh>
  class DataInterpolation
  {

  public:
  
    typedef Mesh                             mesh_Type;
    typedef boost::shared_ptr<mesh_Type>     meshPtr_Type;
    typedef Epetra_CrsMatrix                 matrix_Type;
    typedef boost::shared_ptr<matrix_Type>   matrixPtr_Type;
    typedef VectorEpetra                     vector_Type;
    typedef boost::shared_ptr<vector_Type >  vectorPtr_Type;

    explicit DataInterpolation();
    
    virtual ~DataInterpolation(){}

    //void BuildInterpolationOperator(matrixPtr_Type & InterpolationOperator, meshPtr_Type Mesh1, meshPtr_Type LocalMesh1, vectorPtr_Type KnownField, double radius);

    //void BuildProjectionOperator(matrixPtr_Type ProjectionOperator, meshPtr_Type Mesh2, meshPtr_Type LocalMesh2, double radius);

    //double RBF(double x1, double y1, double x2, double y2, double radius);

  private:
    
    double M_RBFRadius;
    
  };
  
  template <typename Mesh>
  DataInterpolation<Mesh>::DataInterpolation(){}
  
  /*
  template <typename Mesh>
  void DataInterpolation<Mesh>::BuildInterpolationOperator(matrixPtr_Type & InterpolationOperator, meshPtr_Type Mesh1, meshPtr_Type LocalMesh1, std::vector<int> flags, double radius)
  {
    M_RBFRadius = radius;
    GhostHandler<mesh_Type> ghostObj( Mesh1, LocalMesh1, KnownField->mapPtr(), KnownField->mapPtr()->commPtr() );
    ghostObj.setUp();
    
    vector<set<ID> >     Grafo(KnownField->epetraVector().MyLength());
    int * NumeroVicini = new int[KnownField->epetraVector().MyLength()];
    int * GlobalID     = new int[KnownField->epetraVector().MyLength()];
    int max_entries = 0;

    for ( UInt i = 0; i < KnownField->epetraVector().MyLength(); ++i )
      if(KnownField->blockMap().LID(KnownField->blockMap().GID(i)) != -1)
      {
	GlobalID[i] = KnownField->blockMap().GID(i);
	Grafo[i] = ghostObj.createNodeNodeNeighborsMapWithinRadius(radius, KnownField->blockMap().GID(i));  
	Grafo[i].insert(KnownField->blockMap().GID(i));
	NumeroVicini[i] = Grafo[i].size();
	if(NumeroVicini[i]>max_entries)
	  max_entries = NumeroVicini[i];
      }

    Epetra_Map Matrix_Map(Mesh1->numVertices(), KnownField->epetraVector().MyLength(), GlobalID, 0, *(KnownField->mapPtr()->commPtr()));
    InterpolationOperator.reset(new matrix_Type(Copy, Matrix_Map, NumeroVicini));
    
    int * Indices = new int[max_entries];
    double * Values = new double[max_entries];
    int k = 0;
	  
    for( int i = 0 ; i < KnownField->epetraVector().MyLength(); ++i )
    {
        k = 0;
        for( set<ID>::iterator it = Grafo[i].begin(); it != Grafo[i].end(); ++it)
        {
            Indices[k] = *it;
            Values[k]  = RBF( Mesh1->point(GlobalID[i]).x(),
                              Mesh1->point(GlobalID[i]).y(),
                              Mesh1->point(*it).x(),
                              Mesh1->point(*it).y(),
                              M_RBFRadius);
            ++k;
        }
        InterpolationOperator->InsertGlobalValues(GlobalID[i], k, Values, Indices);
    }
    
    InterpolationOperator->FillComplete();
  
    delete Indices;
    delete Values;
    delete NumeroVicini;
    delete GlobalID;

  }

  template <typename Mesh>
  void DataInterpolation<Mesh>::BuildProjectionOperator(matrixPtr_Type ProjectionOperator, meshPtr_Type Mesh2, meshPtr_Type LocalMesh2, double radius)
  {
    
  }

  template <typename Mesh>
  double DataInterpolation<Mesh>::RBF(double x1, double y1, double x2, double y2, double radius)
  {
    double distance = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
    return pow(1-distance/radius,4)*(4*distance/radius+1);
  }
  */

} // namespace LifeV

#endif // DATA_INTERPOLATION_HPP
