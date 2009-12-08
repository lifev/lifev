/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

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

#include <lifemc/lifesolver/Monolithic.hpp>

#include <lifemc/lifealg/eigenSolver.hpp>

#include <EpetraExt_MatrixMatrix.h>

#include <boost/shared_ptr.hpp>

#include <EpetraExt_Reindex_MultiVector.h>

namespace LifeV
{
// Constructors
Monolithic::Monolithic():
    super(),
    M_monolithicMap(),
    M_couplingMatrix(),
    M_interface(0),
    M_interfaceMap(),
    M_monolithicInterfaceMap(),
    M_beta(),
    M_monolithicMatrix(),
    M_precMatrPtr(),
    M_precPtr(),
    M_DDBlockPrec(),
    M_rhsFull(),
    M_BCh_flux(new fluid_bchandler_raw_type),
    M_BCh_Robin(new fluid_bchandler_raw_type),
    M_fluxes(0),
    M_BChWSS(),
    M_bdMass(),
    M_bcfWss(),
    M_robinCoupling(),
    M_alphaf(0.),
    M_alphas(0.),
    M_offset(0),
    M_solidAndFluidDim(0),
    M_solidOper(),
    M_fluidOper(),
    M_meshOper(),
    M_fluidBlock(),
    M_solidBlock(),
    M_solidBlockPrec(),
    M_meshBlock(),
    M_linearSolver(),
    M_wss(),
    //end of protected attributes
    M_PAAP(),
    M_numerationInterface(),
#ifdef OBSOLETE
    M_rhsShapeDerivatives(),
#endif
    M_rhsNew(),
    M_fullMonolithic(),
    M_entry(),
    M_diagonalScale(false),
    M_reusePrec(true),
    M_resetPrec(true),
    M_maxIterSolver(-1)
{}

// Destructor
Monolithic::~Monolithic()
{
}

void
Monolithic::setupFEspace()
{
	super::setupFEspace();

	// Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset( new FESpace<mesh_type, EpetraMap>( M_dataSolid->mesh(),
                                                         M_dataSolid->order(),
                                                         3,
                                                         *M_epetraComm));
}

void
Monolithic::setupDOF( void )
{
	M_dofStructureToHarmonicExtension->setup(   M_uFESpace->refFE(), M_uFESpace->dof(),
											    M_dFESpace->refFE(), M_dFESpace->dof() );
	M_dofStructureToHarmonicExtension->update( *M_uFESpace->mesh(),  M_structureInterfaceFlag,
											   *M_dFESpace->mesh(),  M_harmonicInterfaceFlag,
											    M_interfaceTolerance );

    createInterfaceMaps(M_dofStructureToHarmonicExtension);
}

void Monolithic::setupFluidSolid()
{
    super::setupFluidSolid();


    // Added here the code to build the monolithicMap!!
    // Note: up to now it works only with matching grids (and poly order) on the interface
    assert(M_fluidInterfaceMap->getMap(Unique)->NumGlobalElements() == M_solidInterfaceMap->getMap(Unique)->NumGlobalElements());

    M_interfaceMap = *M_solidInterfaceMap;

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    M_monolithicMap.reset(new EpetraMap(M_uFESpace->map()));
    M_fluxes=M_BCh_flux->size( );
    *M_monolithicMap+= M_pFESpace->map();
    *M_monolithicMap+= M_fluxes;
    *M_monolithicMap+= M_dFESpace->map();

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    int* numInterfaceDof(new int[numtasks]);
    int pid=M_epetraWorldComm->MyPID();
    int numMyElements = M_interfaceMap.getMap(Unique)->NumMyElements();
    numInterfaceDof[pid]=numMyElements;
    EpetraMap subMap(*M_interfaceMap.getMap(Unique), (UInt)0, M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions);

    M_numerationInterface.reset(new vector_type(subMap,Unique));
    //should be an int vector instead of double
    //                    M_numerationInterfaceInt.reset(new Epetra_IntVector(*M_interfaceMap.getMap(Unique)));

    for(int j=0; j<numtasks; ++j)
        M_epetraWorldComm->Broadcast( &numInterfaceDof[j], 1, j);

    for(int j=numtasks-1; j>0 ; --j)
    {
        numInterfaceDof[j] = numInterfaceDof[j-1];
    }
    numInterfaceDof[0]=0;
    for(int j=1; j<numtasks ; ++j)
        numInterfaceDof[j] += numInterfaceDof[j-1];

    UInt k=1;
    UInt l=0;

    M_interface = M_interfaceMap.getMap(Unique)->NumGlobalElements()/nDimensions;
    UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;

    for(l=0, ITrow=locDofMap.begin(); ITrow!=locDofMap.end() ; ++ITrow)
    {
        if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/)>=0)
        {
            (*M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
            //                                    (*M_numerationInterfaceInt)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(M_numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
            if((int)(*M_numerationInterface)(ITrow->second )!=floor(l+1+ numInterfaceDof[pid]/nDimensions+0.2 /*+ dim*localInterface*/) )
                std::cout<<"ERROR! the numeration of the coupling map is not correct"<<std::endl;

            ++l;
        }

    }

    std::vector<int> couplingVector;
    couplingVector.reserve((int)(M_interfaceMap.getMap(Unique)->NumMyElements()));

    for(int dim=0; dim<nDimensions; ++dim)
    {
        for( ITrow=locDofMap.begin(); ITrow!=locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second)>=0)
            {
                couplingVector.push_back((*M_numerationInterface)(ITrow->second /*+ dim * solidDim*/)+ dim * M_interface );
                //couplingVector.push_back((*M_numerationInterfaceInt)[ITrow->second /*+ dim * solidDim*/]+ dim * M_interface );
            }
        }
    }// so the map for the coupling part of the matrix is just Unique

    M_monolithicInterfaceMap = EpetraMap(-1, couplingVector.size(), &couplingVector[0], M_monolithicMap->getMap(Repeated)->IndexBase()/*1*/, *M_epetraWorldComm);
    *M_monolithicMap  += M_monolithicInterfaceMap;

    //std::cout<<"map global elements : "<<M_monolithicMap->getMap(Unique)->NumGlobalElements()<<std::endl;
    //            std::cout<<"map My elements : "<<M_monolithicMap.getMap(Unique)->NumMyElements()<<std::endl;

    //            std::cout<<"newmap global elements : "<<newMap.getMap(Unique)->NumGlobalElements()<<std::endl;
    //          std::cout<<"newmap My elements : "<<newMap.getMap(Unique)->NumMyElements()<<std::endl;

    //          std::cout<<"repeated newmap global elements : "<<newMap.getMap(Repeated)->NumGlobalElements()<<std::endl;
    //          std::cout<<"repeated newmap My elements : "<<newMap.getMap(Repeated)->NumMyElements()<<std::endl;


    //the map for the interface coupling matrices should be done with respect to the coarser mesh.
    M_beta.reset  (new vector_type(/*M_monolithicMap*/M_uFESpace->map()));

    M_offset = M_uFESpace->dof().numTotalDof()*nDimensions + M_fluxes +  M_pFESpace->dof().numTotalDof();
    M_solidAndFluidDim= M_offset + M_dFESpace->dof().numTotalDof()*nDimensions;
    M_BCh_d->setOffset(M_offset);
    M_BCh_flux->setOffset(M_offset-M_fluxes);
    std::vector<BCBase>::iterator fluxIt = M_BCh_flux->begin( );
    for ( UInt i = 0; i < M_fluxes; ++i, ++fluxIt )
        fluxIt->setOffset( i );

    M_BCh_Robin->setOffset(M_offset);

    if(!M_fullMonolithic)
    {
        M_meshMotion.reset(new FSIOperator::meshmotion_raw_type(*M_mmFESpace,
                                                                *M_epetraComm));

        M_fluid.reset(new FSIOperator::fluid_raw_type(dataFluid(),
                                                      *M_uFESpace,
                                                      *M_pFESpace,
                                                      *M_epetraComm,
                                                      *M_monolithicMap));

        //             if (isLinearFluid())// to be implemented
        //                 M_fluidLin.reset(new FSIOperator::fluidlin_raw_type(dataFluid(),
        //                                                                    *M_uFESpace,
        //                                                                    *M_pFESpace,
        //                                                                    *M_epetraComm));


        M_rhs.reset(new vector_type(*this->M_monolithicMap));
        M_rhsFull.reset(new vector_type(*this->M_monolithicMap));
        M_un.reset (new vector_type(*this->M_monolithicMap));
        M_beta.reset  (new vector_type(/*M_monolithicMap*/M_uFESpace->map()));

        M_solid.reset(new FSIOperator::solid_raw_type(dataSolid(),
                                                      *M_dFESpace,
                                                      *M_epetraComm,
                                                      *M_monolithicMap,
                                                      M_offset
                                                      ));

        //             if (isLinearSolid())// to be implemented with the offset
        //                 M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
        //                                                                    *M_dFESpace,
        //
        //                                                      *M_epetraComm));
    }
}//end setup


void
Monolithic::setDataFromGetPot( GetPot const& data_file )
{
    super::setDataFromGetPot(data_file);

    this->M_dataFluid->setSemiImplicit( data_file("problem/semiImplicit", false) );
    this->M_dataFluid->setUseShapeDerivatives( data_file("fluid/useShapeDerivatives", false) );
    M_DDBlockPrec = data_file( "interface/DDBlockPrec",  0 );
    M_fullMonolithic  = !(M_method.compare("fullMonolithic"));
    M_diagonalScale           = data_file( "solid/prec/diagonalScaling",  false );
    M_entry           = data_file( "solid/prec/entry",  0. );
    M_alphaf           = data_file( "interface/alphaf",  0.5 );
    M_alphas           = data_file( "interface/alphas",  0.5 );
}

void
Monolithic::buildSystem()
{
    M_couplingMatrix.reset(new matrix_type(*M_monolithicMap, 1));//since it is constant, we keep this throughout the simulation
    couplingMatrix(M_couplingMatrix);
    M_couplingMatrix->GlobalAssemble();
    M_solidBlock.reset(new matrix_type(*M_monolithicMap, 1));//since it is constant, we keep this throughout the simulation
    M_solid->buildSystem(M_solidBlock);
    M_solidBlock->GlobalAssemble();
    *M_solidBlock *= (M_dataSolid->getTimeStep()*M_solid->rescaleFactor());
    M_solid->rescaleMatrices();
    M_couplingMatrix->GlobalAssemble();

    if(M_DDBlockPrec==7 || M_DDBlockPrec==8)
    {
        EpetraMap pressureMap(M_pFESpace->map());
        if(M_fluxes)
            pressureMap+=M_fluxes;
        M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, 1));
        *M_solidBlockPrec += *M_solidBlock;
        addDiagonalEntries(1., M_solidBlockPrec, M_uFESpace->map() );
        addDiagonalEntries(1., M_solidBlockPrec, pressureMap, M_uFESpace->dof().numTotalDof()*nDimensions );
        addDiagonalEntries(1., M_solidBlockPrec, M_monolithicInterfaceMap, M_solidAndFluidDim);
        if(M_DDBlockPrec==8)
            couplingMatrix(M_solidBlockPrec, 8);
    }

}

void
Monolithic::couplingMatrix(matrix_ptrtype & bigMatrix, int coupling) // not working with non-matching grids
{// coupling from 1 to 15, working as chmod

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    int flag(coupling);

    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow, flag=coupling)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                if(flag-8>=0)//right low
                {
                    bigMatrix->set_mat_inc( M_offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(),(int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (1.) );//right low
                    flag -= 8;
                }
                if(flag-4>=0)// right up
                {
                    bigMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (-1.) );//right up
                    flag -= 4;
                }
                if(flag-2>=0)//low left
                {
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), 1.);//low left
                    flag -= 2;
                }
                if(flag-1>=0)//low right
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), (-1.*M_solid->rescaleFactor()/*/M_dataFluid->timestep()*/));//low right

                bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim , (int)(*this->M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, 0.0);
            }
        }
    }
}


void
Monolithic::robinCoupling(matrix_ptrtype & IdentityMatrix, Real& alphaf, Real& alphas) // not working with non-matching grids
{
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //Real* entry(new Real(1.));

    addDiagonalEntries(1., IdentityMatrix, *M_monolithicMap);
    //addDiagonalEntries(alphaf, IdentityMatrix, *M_solidInterfaceMap, M_offset, true);
    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), alphas);//low right
                IdentityMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, alphaf );//right up
                IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), alphas);//low left
            }
        }
    }
}


void Monolithic::zeroBlock( matrix_ptrtype matrixPtr,vector_type& colNumeration , const std::map<ID, ID>& map, UInt rowOffset, UInt colOffset)
{//to improve
    std::map<ID, ID> const& Map = M_dofStructureToHarmonicExtension->locDofMap();
    Real* entry(new Real(0.));
    std::map<ID, ID>::const_iterator ITrow;
    int err;
    if(colOffset>=M_offset && colOffset<M_solidAndFluidDim && rowOffset<M_solidAndFluidDim)
    {
        for( ITrow=map.begin(); ITrow != map.end(); ++ITrow)// scalable loops,
        {

            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second) >= 0 )//to avoid repeated stuff
            {
                int index=(colOffset + ITrow->second);// second is column
                err=matrixPtr->getMatrixPtr()->ReplaceGlobalValues(rowOffset + ITrow->first,1, entry ,&index);
            }
        }
        return;
    }
    if(rowOffset>=M_offset && rowOffset<M_solidAndFluidDim && colOffset<M_solidAndFluidDim)
    {
        for( ITrow=map.begin(); ITrow != map.end(); ++ITrow)// scalable loops,
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second) >= 0 )//to avoid repeated stuff
            {
                int index2= (colOffset + ITrow->first);
                err=matrixPtr->getMatrixPtr()->ReplaceGlobalValues(rowOffset + ITrow->second,1, entry ,&index2);
            }
        }
        return;
    }

    if(colOffset>=M_solidAndFluidDim && rowOffset < M_offset)
    {
        for( ITrow=map.begin(); ITrow != map.end(); ++ITrow)// scalable loops,
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second) >= 0 )//to avoid repeated stuff
            {
                int* index(new int(colOffset + colNumeration[ITrow->second]));
                err= matrixPtr->getMatrixPtr()->ReplaceGlobalValues(rowOffset + ITrow->first,1, entry ,index);

            }
        }
        return;
    }

    if((colOffset>=M_solidAndFluidDim && rowOffset >= M_offset))
    {
        for( ITrow=map.begin(); ITrow != map.end(); ++ITrow)// scalable loops,
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                int* index2(new int(colOffset + colNumeration[ITrow->second]));
                err=matrixPtr->getMatrixPtr()->ReplaceGlobalValues( rowOffset + ITrow->second,1, entry ,index2);
            }
            if (err != 0)
            {
                M_solid->getDisplayer().leaderPrint("error ", err);
                M_solid->getDisplayer().leaderPrint("Go to fill your matrix!!");
            }
        }
        return;
    }

}


void
Monolithic::couplingRhs(vector_ptrtype rhs, vector_ptrtype un) // not working with non-matching grids
{
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //    UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;

    vector_type lambda(M_interfaceMap, Unique);
    this->monolithicToInterface(lambda, *un);

    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {

                if(rhs.get() != 0)
                    (*rhs)[  (int)(*this->M_numerationInterface)[ITrow->second ] + dim*M_interface +M_solidAndFluidDim ] = -lambda( ITrow->second + dim*M_dFESpace->dof().numTotalDof() )*(M_solid->rescaleFactor());
            }
        }
    }
}


void
Monolithic::addDiagonalEntries(Real entry, matrix_ptrtype bigMatrix, const EpetraMap& Map, UInt offset, bool replace)
{
    if(!replace)
        for(UInt i=0 ; i<Map.getMap(Unique)->NumMyElements(); ++i)//num from 1
        {
            bigMatrix->set_mat_inc(  offset + Map.getMap(Unique)->GID(i)-1 ,   offset + Map.getMap(Unique)->GID(i)-1, entry);
        }
    else
    {
        for(UInt i=0 ; i<Map.getMap(Unique)->NumMyElements(); ++i)
        {
            //diagonal[Map.getMap(Repeated)->GID(i)+offset]=entry;
            int* index(new int(offset + Map.getMap(Unique)->GID(i)));
            bigMatrix->getMatrixPtr()->ReplaceGlobalValues(offset + Map.getMap(Repeated)->GID(i), 1, &entry, index );
        }
    }
}
void
Monolithic::addDiagonalEntries(Real entry,  matrix_ptrtype bigMatrix, std::map<ID, ID> const& Map, UInt offset, bool replace)
{
    std::map<ID, ID>::const_iterator ITrow;
    if(!replace)
        for(ITrow=Map.begin() ; ITrow != Map.end(); ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                bigMatrix->set_mat_inc(  offset + (*M_numerationInterface)[ITrow->second]-1 ,   offset + (*M_numerationInterface)[ITrow->second]-1, entry);
        }
    else
    {
        for(ITrow=Map.begin() ; ITrow != Map.end(); ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                int* index(new int(offset + (*M_numerationInterface)[ITrow->second]));
                int err = bigMatrix->getMatrixPtr()->ReplaceGlobalValues(offset + (*M_numerationInterface)[ITrow->second], 1, &entry, index );
                if (err != 0)
                {
                    M_solid->getDisplayer().leaderPrint("error ", err);
                    M_solid->getDisplayer().leaderPrint("Go to fill your matrix!!");
                }
            }
        }
    }
}

void
Monolithic::updateSystem(const vector_type& displacement)
{

    vector_type solution(*this->M_monolithicMap);
    monolithicToX(displacement, solution, M_uFESpace->map(), UInt(0));
    this->M_bdf->shift_right(solution);

    M_meshMotion->updateSystem();

    *this->M_un                = displacement;
    this->fluid().updateUn(*this->M_un);
    *M_rhs*=0;
    *M_rhsFull*=0;
    this->M_fluid->resetStab();
}

void
Monolithic::monolithicToX(const vector_type& disp, vector_type& dispFluid, EpetraMap& map, UInt offset)
{
    if(disp.getMaptype()== Repeated)
    {
        vector_type dispUnique(disp, Unique);
        monolithicToX(dispUnique, dispFluid, map, offset);
        dispFluid = dispUnique;
        return;
    }
    dispFluid.subset(disp, map, offset, offset);
}

void
Monolithic::monolithicToInterface(vector_type& lambdaSolid, const vector_type& disp)
{
    if (disp.getMaptype() == Repeated)
    {
        vector_type const  dispUnique(disp, Unique);
        monolithicToInterface(lambdaSolid, dispUnique);
        return;
    }
    if (lambdaSolid.getMaptype() == Repeated)
    {
        vector_type  lambdaSolidUn(lambdaSolid.getMap(), Unique);
        monolithicToInterface( lambdaSolidUn, disp);
        lambdaSolid = lambdaSolidUn;
        return;
    }
    /* UInt MyOffset(M_uFESpace->map().getMap(Unique)->NumMyElements()+M_pFESpace->map().getMap(Unique)->NumMyElements());
       vector_type subDisp(this->M_dFESpace->map(), Unique);
       subDisp.mySubset(disp, MyOffset);
       lambdaSolid=subDisp;*/

    EpetraMap subMap(*disp.getMap().getMap(Unique), M_offset,disp.getMap().getMap(Unique)->NumGlobalElements() );
    vector_type subDisp(subMap, Unique);
    subDisp.subset(disp, M_offset);
    lambdaSolid=subDisp;
}


void Monolithic::setDispSolid(const vector_type &sol)
{
    vector_type disp(*M_monolithicMap);
    monolithicToX(sol, disp, M_dFESpace->map(), M_offset);
    this->M_solid->setDisp(disp);
}

void
Monolithic::evalResidual( vector_type&       res,
                          const vector_type& disp,
                          const UInt          iter )
{
    //            eval(disp, iter);


    if((iter==0)|| !this->M_dataFluid->isSemiImplicit())
    {
        setDispSolid(disp);

        vector_type lambdaFluid(this->M_interfaceMap, Unique);

        monolithicToInterface(lambdaFluid, disp);

        lambdaFluid *= (M_dataFluid->getTimeStep()*(M_solid->rescaleFactor()));//because of the matrix scaling
        this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface
        M_meshMotion->iterate(*M_BCh_mesh);
        M_meshMotion->updateDispDiff();

        M_beta.reset(new vector_type(M_uFESpace->map()));
        vector_type meshDispDiff( M_meshMotion->disp(), Repeated );

        this->moveMesh(meshDispDiff);//initialize the mesh position with the total displacement

        meshDispDiff=M_meshMotion->dispDiff();//repeating the mesh dispDiff
        this->interpolateVelocity(meshDispDiff, *this->M_beta);

        double alpha = 1./M_dataFluid->getTimeStep();//mesh velocity w

        *this->M_beta *= -alpha;
        vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
        fluid->subset(*M_un, 0);
        *this->M_beta += *fluid/*M_un*/;//relative velocity beta=un-w


        M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap));
        if(M_DDBlockPrec==1)
        {
            M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
            M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, M_fluidBlock );
            M_fluidBlock->GlobalAssemble();
            *M_monolithicMatrix += *M_fluidBlock;
        }
        else
            if( M_DDBlockPrec==7 || M_DDBlockPrec==8 )
            {
                matrix_ptrtype newMatrix(new matrix_type(*M_monolithicMap));
                M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, newMatrix );
                newMatrix->GlobalAssemble();
                *M_monolithicMatrix += *newMatrix;

                M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
                *M_fluidBlock += *newMatrix;
                this->M_fluid->updateStab( *M_fluidBlock);//applies the stabilization terms
                addDiagonalEntries(1., M_fluidBlock, M_dFESpace->map(), M_offset);
                if(M_DDBlockPrec==7)
                    couplingMatrix(M_fluidBlock, 7);
                else
                    if(M_DDBlockPrec==8)
                        couplingMatrix(M_fluidBlock, 6);

                if ( !M_BCh_flux->bdUpdateDone() )
                    M_BCh_flux->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
                bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_flux, M_uFESpace->feBd(), 1., dataSolid().getTime() );
                M_fluidBlock->GlobalAssemble();
                M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            }
            else
            {
                M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, M_monolithicMatrix );//here it assembles the fluid matrices
            }



        M_fluid->updateStab( *M_monolithicMatrix);//applies the stabilization terms

        updateMatrix(*M_monolithicMatrix);

        //M_monolithicMatrix->GlobalAssemble();
        if(iter==0)
        {
            M_nbEval = 0; // new time step
            M_resetPrec=true;
            *this->M_rhs               += M_fluid->matrMass()*M_bdf->time_der( M_dataFluid->getTimeStep() );
            couplingRhs(this->M_rhs, M_un);
            this->M_solid->updateStuff();
            updateSolidSystem(this->M_rhs);
        }

        M_rhsFull.reset(new vector_type(*M_rhs));

        M_BCh_flux->setOffset(M_offset-M_fluxes);
        if ( !M_BCh_flux->bdUpdateDone() )
            M_BCh_flux->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
        bcManage( *M_monolithicMatrix, *this->M_rhsFull, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_flux, M_uFESpace->feBd(), 1., dataSolid().getTime() );

        M_BCh_Robin->setOffset(M_offset);
        if ( !M_BCh_Robin->bdUpdateDone() )
            M_BCh_Robin->bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
        bcManage( *M_monolithicMatrix, *this->M_rhsFull, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime() );

        evalResidual( *M_BCh_u, *M_BCh_d, disp, M_rhsFull, res, M_diagonalScale);

        if(M_DDBlockPrec>=2 && M_DDBlockPrec!=5 && M_DDBlockPrec!=7 && M_DDBlockPrec!=8)
        {
            if(!M_robinCoupling.get())
            {
                M_robinCoupling.reset( new matrix_type(*M_monolithicMap));
                robinCoupling(M_robinCoupling, M_alphaf, M_alphas);
                M_robinCoupling->GlobalAssemble();
            }
            this->applyPreconditioner(M_robinCoupling, *M_rhsFull);
        }

        if(M_DDBlockPrec!=5 && M_DDBlockPrec!=7 &&M_DDBlockPrec!=8)
            M_precMatrPtr.reset(new matrix_type(*M_monolithicMap));
        M_nbEval++ ;
    }
    evalResidual( disp,  M_rhsFull, res, false);
}

#if OBSOLETE
void Monolithic::shapeDerivatives(vector_ptrtype rhs, vector_ptrtype meshDeltaDisp, const vector_type& sol)
{
    double alpha = 1./M_dataFluid->getTimeStep();
    vector_type derVeloFluidMsh(*meshDeltaDisp, Repeated);
    derVeloFluidMsh *= alpha;
    vector_ptrtype rhsNew(new vector_type(*M_monolithicMap));
    vector_type un(M_uFESpace->map()/*+M_pFESpace->map()*/);
    vector_type uk(M_uFESpace->map()+M_pFESpace->map());
    vector_type meshVel(M_meshMotion->dispDiff(), Repeated);
    un.subset(*this->M_un, 0);
    uk.subset(sol, 0);
    //    if(!M_rhsNew.get())
    //    M_rhsNew.reset(new vector_type(M_uFESpace->map()+M_pFESpace->map()));//useless
    vector_type dvfm(M_uFESpace->map()/*+M_pFESpace->map()*/, Repeated);
    vector_type vfm(M_uFESpace->map()/*+M_pFESpace->map()*/, Repeated);
    vector_type vmdisp(M_uFESpace->map()/*+M_pFESpace->map()*/, Repeated);
    this->transferMeshMotionOnFluid(derVeloFluidMsh, dvfm);
    this->transferMeshMotionOnFluid(*meshDeltaDisp, vmdisp);
    this->transferMeshMotionOnFluid(meshVel, vfm);

    //    rhs->spy("rhsSd");
    //    meshDeltaDisp->spy("deltaDispSd");
    //    sol.spy("solSd");
    //    un.spy("unSd");
    //    uk.spy("ukSd");
    //    vmdisp.spy("deltax");
    //    vfm.spy("velFluidMesh");
    //    dvfm.spy("derVelFluidMesh");

    M_fluid->updateLinearSystem(*M_monolithicMatrix,
                                alpha,
                                un,//un
                                uk,//uk
                                vmdisp,//unknown x.
                                vfm, //(xk-xn)/dt //Repeated
                                dvfm, // (x)/dt //Repeated
                                *rhsNew// useless
                                );
    //    M_rhsShapeDerivatives.reset(new vector_type(*this->M_rhs)) ;
    *rhs = M_fluid->rhsLinNoBC();// Import
    //    rhs->spy("rhsShape");
    //    std::cout<<"rhs lin no bc : "<<M_fluid->rhsLinNoBC().NormInf()<<std::endl;
}
#endif

void Monolithic::shapeDerivatives(matrix_ptrtype sdMatrix, const vector_type& sol, bool domainVelImplicit, bool convectiveTermDer)
{
    double alpha = 1./M_dataFluid->getTimeStep();
    vector_ptrtype rhsNew(new vector_type(*M_monolithicMap));
    vector_type un(M_uFESpace->map()/*+M_pFESpace->map()*/);
    vector_type uk(M_uFESpace->map()+M_pFESpace->map());

    //vector_type meshVel(M_meshMotion->dispDiff(), Repeated);
    vector_ptrtype meshVel;//(M_mmFESpace->map(), Repeated);
    meshVel.reset(new vector_type(M_mmFESpace->map()));

    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);
    if(domainVelImplicit)
    {
        vector_type meshDispOld(M_mmFESpace->map());
        meshVel->subset(sol, offset); //if the conv. term is to be condidered implicitly
        meshDispOld.subset(*M_un, offset);
        *meshVel -= meshDispOld;
    }
    else
    {
        meshVel->subset(*M_un, offset); //if the conv. term is to be condidered partly explicitly
        *meshVel -= M_meshMotion->dispOld();
    }

    if(convectiveTermDer)
        un.subset(sol, 0);
    else
        un.subset(*this->M_un, 0);


    *meshVel *= alpha;
    vector_ptrtype meshVelRep(new vector_type(M_mmFESpace->map(), Repeated));
    *meshVelRep = *meshVel;

    uk.subset(sol, 0);
    vector_type dvfm(M_uFESpace->map(), Repeated);
    vector_type vfm(M_uFESpace->map(), Repeated);
    //vector_type vmdisp(M_uFESpace->map(), Repeated);
    this->transferMeshMotionOnFluid(*meshVelRep, vfm);

    M_fluid->updateShapeDerivatives(*sdMatrix,
                                    alpha,
                                    un,//un if !domainVelImplicit, otherwise uk
                                    uk,//uk
                                    vfm, //(xk-xn)/dt (FI), or (xn-xn-1)/dt (CE)//Repeated
                                    M_solidAndFluidDim+M_interface*nDimensions,
                                    *M_uFESpace,
                                    domainVelImplicit,
                                    convectiveTermDer
                                    );
}

void Monolithic::
evalResidual( const vector_type& sol,const vector_ptrtype& rhs, vector_type& res, bool diagonalScaling)
{
    if(diagonalScaling)
        diagonalScale(*rhs, M_monolithicMatrix);
    res = *M_monolithicMatrix*sol;
    res -= *rhs;
    // Ax-b
}

void Monolithic::
evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, vector_ptrtype& rhs, vector_type& res, bool diagonalScaling, matrix_ptrtype preconditioner)
{

    /* New version to multiply matrices: Paolo should test this and replace the code

       matrix_ptrtype tmpMatPtr(new matrix_type(M_solid->getMap()));
       M_monolithicMatrix->swapCrsMatrix(*tmpMatPtr);

       tmpMatPtr -> GlobalAssemble();

       preconditioner->Multiply( false,
       *tmpMatPtr,
       false,
       *M_monolithicMatrix);

       */


    M_monolithicMatrix->GlobalAssemble();

    matrix_ptrtype tmpMatPtr(new matrix_type(*M_monolithicMatrix));
    tmpMatPtr->GlobalAssemble();
    M_monolithicMatrix.reset(new matrix_type(M_solid->getMap()));

    int err = EpetraExt::MatrixMatrix::
        Multiply( *preconditioner->getMatrixPtr(),
                  false,
                  *tmpMatPtr->getMatrixPtr(),
                  false,
                  *M_monolithicMatrix->getMatrixPtr()
                  );

    *rhs = (*preconditioner)*(*rhs);

    evalResidual( bchFluid, bchSolid, sol, rhs, res, diagonalScaling);
}

void Monolithic::
evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, vector_ptrtype& rhs, vector_type& res, bool diagonalScaling)
{
    vector_type rhsFullSolid(*rhs, Unique); // ignoring non-local entries, Otherwise they are summed up lately

    if(M_solid->offset())
        bchSolid.setOffset(M_solid->offset());
    if ( !bchSolid.bdUpdateDone() )
        bchSolid.bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );

    bcManage( *M_monolithicMatrix, rhsFullSolid, *M_dFESpace->mesh(), M_dFESpace->dof(), bchSolid, M_dFESpace->feBd(), 1.,
              dataSolid().getTime() );

    // matrix is GlobalAssembled by  bcManage

    if ( !bchFluid.bdUpdateDone() )
        bchFluid.bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );

    vector_type rhsFull(rhsFullSolid);
    bcManage( *M_monolithicMatrix, rhsFull, *M_uFESpace->mesh(), M_uFESpace->dof(), bchFluid, M_uFESpace->feBd(), 1.,
              dataSolid().getTime() );

    M_solid->getDisplayer().leaderPrint("rhs norm = ", rhs->NormInf() );

    *rhs = rhsFull;

    evalResidual(sol,rhs, res, diagonalScaling);
}

void  Monolithic::setupBlockPrec(vector_type& rhs)
{
    boost::shared_ptr<IfpackComposedPrec>  ifpackCompPrec;


    if( !M_reusePrec || M_resetPrec )
    {
        switch(M_DDBlockPrec)
        {
        case 1:
            *M_precMatrPtr += *M_fluidBlock;
            *M_precMatrPtr += *M_solidBlock;
            this->M_fluid->updateStab( *M_precMatrPtr);//applies the stabilization terms
            if(!M_precMatrPtr->getMatrixPtr()->Filled())
            {
                couplingMatrix(M_precMatrPtr, 7);
                //addDiagonalEntries(M_entry,M_precMatrPtr, M_interfaceMap, M_solidAndFluidDim);
            }
            else
                this->M_solid->getDisplayer().leaderPrint("ERROR: probably the tolerance fixed for Newton is too low. Remember that this type of preconditioner can be used only in the geometry-explicit case \n");
            bcManageMatrix( *M_precMatrPtr, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_flux, M_uFESpace->feBd(), 1., dataSolid().getTime() );
            bcManageMatrix( *M_precMatrPtr, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime() );

            M_precMatrPtr->GlobalAssemble();
            bcManageMatrix( *M_precMatrPtr, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1.,
                            dataSolid().getTime() );
            bcManageMatrix( *M_precMatrPtr, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_u, M_uFESpace->feBd(), 1.,
                            dataFluid().getTime() );

            M_precMatrPtr->GlobalAssemble();
            break;

        case 2:
            *M_precMatrPtr += *M_monolithicMatrix;
            M_precMatrPtr->GlobalAssemble();
            break;

        case 3:
            if(false && !M_robinCoupling.get())
            {
                M_robinCoupling.reset( new matrix_type(*M_monolithicMap));
                robinCoupling(M_robinCoupling, M_alphaf, M_alphas);
                M_robinCoupling->GlobalAssemble();
            }
            *M_precMatrPtr += *M_monolithicMatrix;
            M_precMatrPtr->GlobalAssemble();

            {
                std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
                for(short i=0; i<nDimensions; ++i)
                {
                    zeroBlock(M_precMatrPtr/*tmpPrecPtr*/, *M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_dFESpace->dof().numTotalDof()+M_offset); //(2, 4)
                    zeroBlock(M_precMatrPtr/*tmpPrecPtr*/, *M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);//(2, 5)
                }
            }

            break;

        case 4:
            *M_precMatrPtr += *M_monolithicMatrix;
            M_precMatrPtr->GlobalAssemble();

            {
                std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
                for(short i=0; i<nDimensions; ++i)
                {
                    zeroBlock(M_precMatrPtr, *M_numerationInterface, locDofMap, M_offset + i*M_dFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);//(4, 5)
                }
            }
            break;

        case 5:
            break;
        case 6:
            {
            displayer().leaderPrint("Preconditioner type not yet implemented," );
            displayer().leaderPrint("change the entry DDBlockPrec in the data file");
            throw WRONG_PREC_EXCEPTION();
            break;
            }
        case 7:

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_precPtr);

            M_BCh_flux->setOffset(M_offset-M_fluxes);
            M_BCh_Robin->setOffset(M_offset);

            if(ifpackCompPrec->set())
            {
                ifpackCompPrec->replace(M_fluidOper, 1);
            }
            else
            {
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime());
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().getTime() );
                M_solidBlockPrec->GlobalAssemble();

                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));

                ifpackCompPrec->buildPreconditioner(M_solidOper);
                ifpackCompPrec->push_back(M_fluidOper);
            }

            //
            break;

        case 8:

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_precPtr);

            M_BCh_flux->setOffset(M_offset-M_fluxes);
            M_BCh_Robin->setOffset(M_offset);


            if(ifpackCompPrec->set())
            {
                ifpackCompPrec->replace(M_fluidOper, 0);
            }
            else
            {
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime() );
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().getTime() );
                M_solidBlockPrec->GlobalAssemble();

                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));

                ifpackCompPrec->buildPreconditioner(M_fluidOper);
                ifpackCompPrec->push_back(M_solidOper);
            }

            break;

        case 9:
        case 10:
        case 11:
        case 12:
            displayer().leaderPrint("Preconditioner type not yet implemented," );
            displayer().leaderPrint("change the entry DDBlockPrec in the data file");
            throw WRONG_PREC_EXCEPTION();
            break;
        default:
            {
                for(short i=0; i<nDimensions; ++i)
                    *M_precMatrPtr += *M_monolithicMatrix;
                M_precMatrPtr->GlobalAssemble();

            }
            break;
        }
    }
}

void  Monolithic::solveJac(vector_type         &_step,
                           const vector_type   &_res,
                           const Real         /*_linearRelTol*/)
{
    setupBlockPrec(*const_cast<vector_type*>(&_res));
    // #if OBSOLETE
    //     if(!M_dataFluid->useShapeDerivatives())
    // #endif

    //M_precMatrPtr->spy("p");


    M_solid->getDisplayer().leaderPrint("solveJac: NormInf res ", _res.NormInf());

    M_solid->getDisplayer().leaderPrint("Solving Jacobian system... \n" );

    switch(M_DDBlockPrec)
    {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
        this->iterateMonolithic(_res, _step, M_precMatrPtr, M_linearSolver);
        break;

    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        this->iterateMonolithic(_res, _step, M_precPtr,     M_linearSolver);
        break;

    default:
        this->iterateMonolithic(_res, _step, M_precMatrPtr, M_linearSolver);
        break;
    }

    this->M_solid->getDisplayer().leaderPrint("done.\n");
}

void
Monolithic::iterateMesh(const vector_type& disp)
{
    vector_type lambdaFluid(this->M_interfaceMap, Unique);

    monolithicToInterface(lambdaFluid, disp);

    lambdaFluid *= (M_dataFluid->getTimeStep()*(M_solid->rescaleFactor()));

    this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

    M_meshMotion->iterate(*M_BCh_mesh);

}

//void Monolithic::variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
void Monolithic::variablesInit(const std::string& dOrder)
{
    //    EpetraMap interfaceMap(*M_solidInterfaceMap);
    M_solidMeshPart.reset( new  partitionMesh< FSIOperator::mesh_type > (*M_dataSolid->mesh(), *M_epetraComm/*, M_solidInterfaceMap->getMap(Unique).get(), M_solidInterfaceMap->getMap(Repeated).get()*/));

    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_solidMeshPart,
                                                       dOrder,
                                                       //*refFE_struct,
                                                       //*qR_struct,
                                                       //*bdQr_struct,
                                                       3,
                                                       *M_epetraComm));
    // INITIALIZATION OF THE VARIABLES
    M_lambdaFluid.reset(new vector_type(*M_fluidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset(new vector_type(*M_fluidInterfaceMap, Repeated) );

}


void Monolithic::
applyPreconditioner( matrix_ptrtype robinCoupling, vector_type& rhs)
{
    applyPreconditioner(robinCoupling, M_monolithicMatrix);
    rhs=*robinCoupling*(rhs);

}

void Monolithic::
applyPreconditioner( matrix_ptrtype robinCoupling, matrix_ptrtype& prec )
{
    matrix_type tmpMatrix(*M_monolithicMap/* *prec*/);
    EpetraExt::MatrixMatrix::Multiply( *robinCoupling->getMatrixPtr(),
                                       false,
                                       *prec->getMatrixPtr(),
                                       false,
                                       *tmpMatrix.getMatrixPtr());
    prec->swapCrsMatrix(tmpMatrix);
}


void Monolithic::
updateMatrix(matrix_type & bigMatrixStokes)
{
    bigMatrixStokes += *M_couplingMatrix;
    bigMatrixStokes += *M_solidBlock;
}

void Monolithic::
updateSolidSystem( vector_ptrtype & rhsFluidCoupling )
{
    M_solid->updateSystem();
    *rhsFluidCoupling += *M_solid->rhsWithoutBC();
}

// void
// Monolithic::setupBDF( vector_type const& u0)
// {
//     M_bdf.reset(new BdfT<vector_type>(M_dataFluid->getBDF_order()));
//     std::cout<<"initial condition : "<<u0.NormInf()<<std::endl;
//     M_bdf->initialize_unk(u0);
// }

void
Monolithic::setupSystem( )
{
    M_fluid->setUp( M_dataFile );
    M_meshMotion->setUp( M_dataFile );
    setUp( M_dataFile );
}

void
Monolithic::setUp( const GetPot& dataFile )
{
    M_solid->getDisplayer().leaderPrint("\n S-  Displacement unknowns: ",  M_dFESpace->dof().numTotalDof() );
    M_solid->getDisplayer().leaderPrint(" S-  Computing mass and linear strain matrices ... \n");
    M_linearSolver.reset(new solver_type(*M_epetraComm));
    M_reusePrec     = dataFile( "solid/prec/reuse", true);
    //M_maxIterSolver = dataFile( "solid/solver/max_iter", -1);

    M_linearSolver->setDataFromGetPot( dataFile, "solid/solver" );
    M_linearSolver->setUpPrec(dataFile, "solid/prec");//to avoid if we have already a prec.

    // We should use a factory here, like reset( PRECFactory::instance().createObject( precType ) );
    M_precPtr.reset(new IfpackComposedPrec(M_epetraComm.get()));
    M_precPtr->setDataFromGetPot(dataFile, "solid/prec");//to avoid if we build the prec from a matrix.

    M_reusePrec     = dataFile( "solid/prec/reuse", true);
    //    M_maxIterForReuse = data_file( "solid/solver/max_iter_reuse", M_maxIterSolver*8/10);
    M_maxIterSolver = dataFile( "solid/solver/max_iter", -1);
}

void Monolithic::
diagonalScale(vector_type& rhs, matrix_ptrtype matrFull)
{
    Epetra_Vector diagonal(*rhs.getMap().getMap(Unique));
    //M_matrFull->getMatrixPtr()->InvRowSums(diagonal);
    //M_matrFull->getMatrixPtr()->InvRowMaxs(diagonal);
    //M_matrFull->getMatrixPtr()->InvColSums(diagonal);
    matrFull->getMatrixPtr()->InvColMaxs(diagonal);
    matrFull->getMatrixPtr()->LeftScale(diagonal);
    rhs.getEpetraVector().Multiply(1, rhs.getEpetraVector(), diagonal,0);
}


//void Monolithic::solidInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
void Monolithic::solidInit(std::string const& dOrder)
{   // Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataSolid->mesh(),
                                                       dOrder,
                                                       //*refFE_struct,
                                                       //*qR_struct,
                                                       //*bdQr_struct,
                                                       3,
                                                       *M_epetraComm));
}


void
Monolithic::setFluxBC             (fluid_bchandler_type bc_flux)
{
    if (isFluid())
    {
        M_BCh_flux = bc_flux;
    }
}


void
Monolithic::setRobinBC             (fluid_bchandler_type bc_Robin)
{
    if (isFluid())
    {
        M_BCh_Robin = bc_Robin;
    }
}


void Monolithic::initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                             FSIOperator::solid_type::value_type::Function const& p0,
                             FSIOperator::solid_type::value_type::Function const& d0,
                             FSIOperator::solid_type::value_type::Function const& w0,
                             FSIOperator::solid_type::value_type::Function const& /*w0*/ )
{
    vector_type u(M_uFESpace->map());
    M_uFESpace->interpolate(u0, u, M_dataFluid->getTime());

    vector_type p(M_pFESpace->map());
    M_pFESpace->interpolate(p0, p, M_dataFluid->getTime());

    vector_type d(M_dFESpace->map());
    M_dFESpace->interpolate(d0, d, M_dataSolid->getTime());

    //vector_type w(M_dFESpace->map());
    //M_dFESpace->interpolate(w0, w, M_dataSolid->getTime());

    initialize(u, p, d);
    //fluid().initialize(u0,p0);
    //solid().initialize(d0, w0);
}

void Monolithic::initialize( const vector_type& u0, const vector_type& p0, const vector_type& d0)
{
    *M_un=u0;
    M_un->add(p0, nDimensions*M_uFESpace->dof().numTotalDof());
    M_un->add(d0, M_offset);
    //  M_bdf->initialize_unk(*M_un);
}

#ifdef HAVE_TRILINOS_ANASAZI
void Monolithic::computeMaxSingularValue()
{
    typedef Epetra_Operator                                                operator_type;

    M_PAAP.reset(new ComposedPreconditioner<Epetra_Operator>(M_epetraComm.get()));

    boost::shared_ptr<operator_type>  ComposedPrecPtr(M_linearSolver->getPrec()->getPrec());

    M_monolithicMatrix->getMatrixPtr()->OptimizeStorage();
    boost::shared_ptr<Epetra_FECrsMatrix> matrCrsPtr(new Epetra_FECrsMatrix(*M_monolithicMatrix->getMatrixPtr()));

    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(ComposedPrecPtr/*matrCrsPtr*/));
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(/*ComposedPrecPtr*/matrCrsPtr),  true);
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(/*ComposedPrecPtr*/matrCrsPtr), true, true);
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(ComposedPrecPtr), false, true);

    std::vector<LifeV::Real> real;
    std::vector<LifeV::Real> imaginary;

    boost::shared_ptr<EigenSolver> eig;

    UInt nev = M_dataFile("eigensolver/nevec", 10);//number of eigenvectors
    if(nev)
    {
        eig.reset(new EigenSolver(M_PAAP, M_PAAP->OperatorDomainMap(), nev));
        eig->setDataFromGetPot(M_dataFile, "eigensolver/");
        eig->solve();
        eig->eigenvalues(real, imaginary);
    }
    else
    {
        throw UNDEF_EIGENSOLVER_EXCEPTION();
    }
    for (int i=0; i<real.size(); ++i)
    {
        displayer().leaderPrint("\n real part ", real[i]);
        displayer().leaderPrint("\n imaginary part ", imaginary[i]);
    }

}
#endif

namespace
{
static Real fZero(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{return 0.;}
static Real fOne(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{return 1.;}
}

void Monolithic::enableWssComputation(EntityFlag flag)
{
    M_BChWSS.reset(new solid_bchandler_raw_type());
    BCFunctionBase bcfZero(fZero);
    BCFunctionBase bcfOne(fOne);
    M_bcfWss.setFunctions_Mixte(bcfOne,bcfOne);

    M_BChWSS->addBC("WSS", (EntityFlag) flag, Mixte, Full, M_bcfWss, 3);
}

boost::shared_ptr<EpetraVector> Monolithic::computeWS()
{
    //M_BChWSS->setOffset(M_offset);
    M_bdMass.reset(new matrix_type(M_interfaceMap));
    if ( !M_BChWSS->bdUpdateDone() )
        M_BChWSS->bdUpdate(*M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    bcManageMatrix(*M_bdMass, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BChWSS, M_dFESpace->feBd(), 1., dataSolid().getTime() );
    M_bdMass->GlobalAssemble();

    vector_type lambda(M_monolithicInterfaceMap);
    lambda.subset(*M_un, M_solidAndFluidDim);

    solver_type solverMass(*M_epetraComm);
    solverMass.setDataFromGetPot( M_dataFile, "solid/solver" );
    solverMass.setUpPrec(M_dataFile, "solid/prec");//to avoid if we have already a prec.

    boost::shared_ptr<IfpackPreconditioner> P(new IfpackPreconditioner());

    vector_ptrtype sol(new vector_type(M_monolithicInterfaceMap));
    solverMass.setMatrix(*M_bdMass);
    solverMass.setReusePreconditioner(false);
    int numIter = solverMass.solveSystem( lambda, *sol, M_bdMass);

    EpetraExt::MultiVector_Reindex reindexMV(*M_interfaceMap.getMap(Unique));
    boost::shared_ptr<EpetraMap> newMap(new EpetraMap( M_interfaceMap ));
    M_wss.reset(new vector_type(reindexMV(sol->getEpetraVector()), newMap, Unique));
    return M_wss;
}

void Monolithic::initializeMesh(vector_ptrtype fluid_dispOld)
{
    meshMotion().setDisplacement(*fluid_dispOld);
}


namespace
{
FSIOperator* createM(){ return new Monolithic(); }
}
static bool reg = FSIFactory::instance().registerProduct( "monolithic", &createM );

}
