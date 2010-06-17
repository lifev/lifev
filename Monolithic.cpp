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
    M_rhsFull(),
    M_BCh_flux(new fluid_bchandler_raw_type),
    M_BCh_Robin(new fluid_bchandler_raw_type),
    M_fluxes(0),
    M_BChWSS(),
    M_bdMass(),
    M_robinCoupling(),
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
    M_listFluid(),
    M_listSolid(),
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
    M_maxIterSolver(-1),
    M_restarts(false)
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
    M_dFESpace.reset( new FESpace<mesh_type, EpetraMap>( M_data->dataSolid()->dataMesh()->mesh(),
                                                         M_data->dataSolid()->order(),
                                                         nDimensions,
                                                         *M_epetraComm));
}

void
Monolithic::setupDOF( void )
{
	M_dofStructureToHarmonicExtension->setup(   M_uFESpace->refFE(), M_uFESpace->dof(),
											    M_dFESpace->refFE(), M_dFESpace->dof() );
	M_dofStructureToHarmonicExtension->update( *M_uFESpace->mesh(),  M_data->structureInterfaceFlag(),
											   *M_dFESpace->mesh(),  M_data->harmonicInterfaceFlag(),
											    M_data->interfaceTolerance() );

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
Monolithic::setDataFile( const GetPot& dataFile )
{
    super::setDataFile( dataFile );

    M_fullMonolithic   = !(M_data->method().compare("fullMonolithic"));
    M_diagonalScale    = dataFile( "solid/prec/diagonalScaling",  false );
    M_entry            = dataFile( "solid/prec/entry",  0. );
    M_restarts         = dataFile( "exporter/start"  ,  0   );
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
    *M_solidBlock *= (M_data->dataSolid()->dataTime()->getTimeStep()*M_solid->rescaleFactor());
    M_solid->rescaleMatrices();
    M_couplingMatrix->GlobalAssemble();

    if(M_data->DDBlockPreconditioner()==7 || M_data->DDBlockPreconditioner()==8 || M_data->DDBlockPreconditioner()>=13)
    {
        EpetraMap fluidPressureMap(M_uFESpace->map());
        fluidPressureMap+= M_pFESpace->map();
        if(M_fluxes)
            fluidPressureMap+= M_fluxes;

        M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, 1));
        *M_solidBlockPrec += *M_solidBlock;
         if(M_data->DDBlockPreconditioner()>=13)
         {
             matrix_ptrtype swap;
             M_solidBlockPrec->GlobalAssemble();
             *M_solidBlockPrec *=0.5;
             swap=M_solidBlockPrec;
             M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, swap->getMeanNumEntries()));
             *M_solidBlockPrec += *swap;
         }
         addDiagonalEntries(1., M_solidBlockPrec, fluidPressureMap );
         if(!(M_data->DDBlockPreconditioner()==14))
             addDiagonalEntries(1., M_solidBlockPrec, M_monolithicInterfaceMap, M_solidAndFluidDim);
         if(M_data->DDBlockPreconditioner()==8)
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
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), (-1.*M_solid->rescaleFactor()/*/M_data->dataFluid()->timestep()*/));//low right

                bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim , (int)(*this->M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, 0.0);
            }
        }
    }
}


void
Monolithic::robinCoupling( matrix_ptrtype& IdentityMatrix, const Real& alphaf, const Real& alphas, int coupling ) // not working with non-matching grids
{
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //Real* entry(new Real(1.));

    //addDiagonalEntries(alphaf, IdentityMatrix, *M_solidInterfaceMap, M_offset, true);
    if(coupling-4 >= 0)
    {
            for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
            {
                if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    for(UInt dim = 0; dim < nDimensions; ++dim)
                    {
                        IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), alphas);//low right//to enable
                    }
                }
            }
            coupling -= 4;
    }
    if(coupling-2 >= 0)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                for(UInt dim = 0; dim < nDimensions; ++dim)
                {
                    IdentityMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, alphaf );//right up
                }
            }
        }
        coupling -= 2;
    }
    if(coupling-1 >= 0)
    {
            for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
            {
                if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    for(UInt dim = 0; dim < nDimensions; ++dim)
                    {
                    IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), alphas);//low left
                    }
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
Monolithic::updateSystem()
{

    vector_type solution(*this->M_monolithicMap);
    monolithicToX(*this->M_un, solution, M_uFESpace->map(), UInt(0));
    this->M_bdf->shift_right(solution);

    M_meshMotion->updateSystem();

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


    if((iter==0)|| !this->M_data->dataFluid()->isSemiImplicit())
    {
        setDispSolid(disp);

        vector_type lambdaFluid(this->M_interfaceMap, Unique);

        monolithicToInterface(lambdaFluid, disp);

        lambdaFluid *= (M_data->dataFluid()->dataTime()->getTimeStep()*(M_solid->rescaleFactor()));//because of the matrix scaling
        this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface
        M_meshMotion->iterate(*M_BCh_mesh);
        M_meshMotion->updateDispDiff();

        M_beta.reset(new vector_type(M_uFESpace->map()));
        vector_type meshDispDiff( M_meshMotion->disp(), Repeated );

        this->moveMesh(meshDispDiff);//initialize the mesh position with the total displacement

        meshDispDiff=M_meshMotion->dispDiff();//repeating the mesh dispDiff
        this->interpolateVelocity(meshDispDiff, *this->M_beta);

        double alpha = 1./M_data->dataFluid()->dataTime()->getTimeStep();//mesh velocity w

        *this->M_beta *= -alpha;
        vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
        fluid->subset(*M_un, (UInt)0);
        *this->M_beta += *fluid/*M_un*/;//relative velocity beta=un-w


        M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap));

        if(M_data->DDBlockPreconditioner()==1)
        {
            M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
            M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, M_fluidBlock );
            M_fluidBlock->GlobalAssemble();
            *M_monolithicMatrix += *M_fluidBlock;
        }
        else
            if( M_data->DDBlockPreconditioner()==7 || M_data->DDBlockPreconditioner()==8 || M_data->DDBlockPreconditioner()>=13 )
            {
                matrix_ptrtype newMatrix(new matrix_type(*M_monolithicMap));
                M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, newMatrix );
                this->M_fluid->updateStab(*newMatrix);
                newMatrix->GlobalAssemble();
                *M_monolithicMatrix += *newMatrix;

                M_fluidBlock.reset(new matrix_type(*M_monolithicMap));

                if( M_data->DDBlockPreconditioner() >= 13 )
                {
                    *newMatrix *= 0.5;
                    if( M_data->DDBlockPreconditioner() == 13 )
                        addDiagonalEntries(1., M_fluidBlock, M_monolithicInterfaceMap, M_solidAndFluidDim);
                }
                *M_fluidBlock += *newMatrix;

                addDiagonalEntries(1., M_fluidBlock, M_dFESpace->map(), M_offset);
                if(M_data->DDBlockPreconditioner()==7)
                    couplingMatrix(M_fluidBlock, 7);
                else
                    if(M_data->DDBlockPreconditioner()==8)
                        couplingMatrix(M_fluidBlock, 6);

                if ( !M_BCh_flux->bdUpdateDone() )
                    M_BCh_flux->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
                bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_flux, M_uFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
                M_fluidBlock->GlobalAssemble();
                M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            }
            else
            {
                M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs, M_monolithicMatrix );//here it assembles the fluid matrices
            }

        if(!( M_data->DDBlockPreconditioner()==7 || M_data->DDBlockPreconditioner()==8 || M_data->DDBlockPreconditioner()>=13 ))
            this->M_fluid->updateStab(*M_monolithicMatrix);



        updateMatrix(*M_monolithicMatrix);

        //M_monolithicMatrix->GlobalAssemble();
        if(iter==0)
        {
            M_nbEval = 0; // new time step
            M_resetPrec=true;
            *this->M_rhs               += M_fluid->matrMass()*M_bdf->time_der( M_data->dataFluid()->dataTime()->getTimeStep() );
            couplingRhs(this->M_rhs, M_un);

            if (!M_restarts)
            {
                this->M_solid->updateVel();
                M_restarts = false;
            }
            updateSolidSystem(this->M_rhs);
        }

        ASSERT(!(M_data->RobinNeumannCoupling() && (M_data->DDBlockPreconditioner() >= 13)), "disable robinNeumannCoupling!");
            if(M_data->RobinNeumannCoupling())
            {
                M_RNcoupling.reset(new matrix_type(*M_monolithicMap, 0));
                matrix_ptrtype tmpMatrix(new matrix_type(*M_monolithicMap, 0));
                Real zero(0.);
                robinCoupling(M_RNcoupling, zero, M_data->RobinNeumannSolidCoefficient(), 1);
                M_RNcoupling->GlobalAssemble();
                //*M_RNcoupling *= M_data->dataFluid()->getTimeStep();
                M_RNcoupling->spy("rn");
                M_monolithicMatrix->GlobalAssemble();
                int err = EpetraExt::MatrixMatrix::
                    Multiply( *M_RNcoupling->getMatrixPtr(),
                              false,
                              *M_monolithicMatrix->getMatrixPtr(),
                              false,
                              *tmpMatrix->getMatrixPtr()
                              );
                tmpMatrix->GlobalAssemble();
                tmpMatrix->spy("tmp");
                matrix_ptrtype swap;

                if(M_data->DDBlockPreconditioner()>6)
                {
                    *M_fluidBlock += *tmpMatrix;
                    M_fluidBlock->GlobalAssemble();
                    M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));
                    M_fluidBlock->spy("fluidblock");
                }
                swap=M_monolithicMatrix;

                M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap));
                *M_monolithicMatrix  += *swap;
                *M_monolithicMatrix += *tmpMatrix;

                M_RNcoupling->GlobalAssemble();

                this->M_rhs->spy("rb");
                *this->M_rhs += *M_RNcoupling*(*this->M_rhs);
                this->M_rhs->spy("ra");
            }
            else
            {
                if(M_fluidBlock.get())
                {
                    //to move away
                    M_fluidBlock->GlobalAssemble();
                    M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));
                }
            }

        M_rhsFull.reset(new vector_type(*M_rhs));

        M_BCh_flux->setOffset(M_offset-M_fluxes);
        if ( !M_BCh_flux->bdUpdateDone() )
            M_BCh_flux->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
        bcManage( *M_monolithicMatrix, *this->M_rhsFull, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_flux, M_uFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );

        M_BCh_Robin->setOffset(M_offset);
        if ( !M_BCh_Robin->bdUpdateDone() )
            M_BCh_Robin->bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
        bcManage( *M_monolithicMatrix, *this->M_rhsFull, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );

        evalResidual( *M_BCh_u, *M_BCh_d, disp, M_rhsFull, res, M_diagonalScale);

        if(M_data->DDBlockPreconditioner()>=2 && M_data->DDBlockPreconditioner()!=5 && M_data->DDBlockPreconditioner()!=7 && M_data->DDBlockPreconditioner()!=8 && M_data->DDBlockPreconditioner()<13)
        {
            if(!M_robinCoupling.get())
            {
                M_robinCoupling.reset( new matrix_type(*M_monolithicMap));
                addDiagonalEntries(1., M_robinCoupling, *M_monolithicMap);
                robinCoupling(M_robinCoupling, M_data->RobinNeumannFluidCoefficient(), M_data->RobinNeumannSolidCoefficient());
                M_robinCoupling->GlobalAssemble();
            }
            this->applyPreconditioner(M_robinCoupling, *M_rhsFull);
        }

        if(M_data->DDBlockPreconditioner()!=5 && M_data->DDBlockPreconditioner()!=7 &&M_data->DDBlockPreconditioner()!=8 && M_data->DDBlockPreconditioner()<13)
            M_precMatrPtr.reset(new matrix_type(*M_monolithicMap));
        M_nbEval++ ;
    }
    evalResidual( disp,  M_rhsFull, res, false);
}


void Monolithic::
evalResidual( const vector_type& sol,const vector_ptrtype& rhs, vector_type& res, bool diagonalScaling)
{
    M_monolithicMatrix->GlobalAssemble();
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

    if(M_solid->offset())
        bchSolid.setOffset(M_solid->offset());
    if ( !bchSolid.bdUpdateDone() )
        bchSolid.bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );

    bcManage( *M_monolithicMatrix, *rhs, *M_dFESpace->mesh(), M_dFESpace->dof(), bchSolid, M_dFESpace->feBd(), 1.,
              dataFluid().dataTime()->getTime() );

    // matrix is GlobalAssembled by  bcManage

    if ( !bchFluid.bdUpdateDone() )
        bchFluid.bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );

    bcManage( *M_monolithicMatrix, *rhs, *M_uFESpace->mesh(), M_uFESpace->dof(), bchFluid, M_uFESpace->feBd(), 1.,
              dataFluid().dataTime()->getTime() );

    //M_solid->getDisplayer().leaderPrint("rhs norm = ", rhs->NormInf() );

    evalResidual(sol,rhs, res, diagonalScaling);
}

int  Monolithic::setupBlockPrec(vector_type& rhs)
{
    boost::shared_ptr<IfpackComposedPrec>  ifpackCompPrec;
    boost::shared_ptr<ComposedPreconditioner<Epetra_Operator> >  firstCompPrec;
    boost::shared_ptr<ComposedPreconditioner<Epetra_Operator> >  secondCompPrec;

    if( !M_reusePrec || M_resetPrec )
    {
        switch(M_data->DDBlockPreconditioner())
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
            bcManageMatrix( *M_precMatrPtr, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_flux, M_uFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
            bcManageMatrix( *M_precMatrPtr, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );

            M_precMatrPtr->GlobalAssemble();
            bcManageMatrix( *M_precMatrPtr, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1.,
                            dataFluid().dataTime()->getTime() );
            bcManageMatrix( *M_precMatrPtr, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_u, M_uFESpace->feBd(), 1.,
                            dataFluid().dataTime()->getTime() );

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
                addDiagonalEntries(1., M_robinCoupling, *M_monolithicMap);
                robinCoupling(M_robinCoupling, M_data->RobinNeumannFluidCoefficient(), M_data->RobinNeumannSolidCoefficient());
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

            M_BCh_Robin->setOffset(M_offset);

            if(ifpackCompPrec->set())
            {
                ifpackCompPrec->replace(M_fluidOper, 1);
            }
            else
            {
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime());
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
                M_solidBlockPrec->GlobalAssemble();

                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));

                ifpackCompPrec->buildPreconditioner(M_solidOper);
                ifpackCompPrec->push_back(M_fluidOper);
            }

            //
            break;

        case 8:

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_precPtr);

            M_BCh_Robin->setOffset(M_offset);


            if(ifpackCompPrec->set())
            {
                ifpackCompPrec->replace(M_fluidOper, 0);
            }
            else
            {
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
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
        case 13:
            {
            firstCompPrec.reset(new ComposedPreconditioner<Epetra_Operator>(M_epetraComm.get()));
            M_CompPrec.reset(new ComposedPreconditioner<Epetra_Operator>(M_epetraComm.get()));
            secondCompPrec.reset(new ComposedPreconditioner<Epetra_Operator>(M_epetraComm.get()));

            M_BCh_Robin->setOffset(M_offset);

            M_couplingFluid1.reset(new matrix_type(*M_monolithicMap, 1));
            M_couplingFluid2.reset(new matrix_type(*M_monolithicMap, 1));
            M_couplingSolid1.reset(new matrix_type(*M_monolithicMap, 1));
            M_couplingSolid2.reset(new matrix_type(*M_monolithicMap, 1));

            matrix_ptrtype tmpMatrix(new matrix_type(*M_monolithicMap, 0));
            addDiagonalEntries(1., M_couplingFluid1, *M_monolithicMap);
            addDiagonalEntries(1., M_couplingFluid2, *M_monolithicMap);
            addDiagonalEntries(1., M_couplingSolid1, *M_monolithicMap);
            EpetraMap fluidPressureMap(M_uFESpace->map());
            fluidPressureMap+= M_pFESpace->map();
            if(M_fluxes)
                fluidPressureMap+= M_fluxes;
            addDiagonalEntries(1., M_couplingSolid2, fluidPressureMap );
            addDiagonalEntries(1., M_couplingSolid2, M_dFESpace->map(), M_offset);
            addDiagonalEntries(-1., M_couplingSolid2, M_monolithicInterfaceMap, M_solidAndFluidDim);

            // fluid part of the preconditioner
            couplingMatrix(tmpMatrix, 4);
            tmpMatrix->GlobalAssemble();
            *tmpMatrix *= -1.;//-0.5;
            tmpMatrix->spy("upright");
            *M_couplingFluid1 += *tmpMatrix;
            M_couplingFluid1->GlobalAssemble();

            tmpMatrix.reset(new matrix_type(*M_monolithicMap, 0));
            couplingMatrix(tmpMatrix, 2);
            tmpMatrix->GlobalAssemble();
            *tmpMatrix *= -1.;// -0.5;
            tmpMatrix->spy("lowleft");
            *M_couplingFluid2 += *tmpMatrix;
            M_couplingFluid2->GlobalAssemble();

            //TODO: diagonalize the coupling matrices for the Dirichlet b.c.? no

            bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
            M_solidBlockPrec->GlobalAssemble();
            bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
            M_solidBlockPrec->GlobalAssemble();

            // solid part of the preconditioner
            tmpMatrix.reset(new matrix_type(*M_monolithicMap, 0));
            couplingMatrix(tmpMatrix, 8);
            tmpMatrix->GlobalAssemble();
            *tmpMatrix *= -1.;
            tmpMatrix->spy("rightlow");
            *M_couplingSolid1 += *tmpMatrix;
            M_couplingSolid1->GlobalAssemble();

            tmpMatrix.reset(new matrix_type(*M_monolithicMap, 0));
            couplingMatrix(tmpMatrix, 1);
            tmpMatrix->GlobalAssemble();
            tmpMatrix->spy("lowright");
            *M_couplingSolid2 += *tmpMatrix;
            M_couplingSolid2->GlobalAssemble();

            M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));

            Chrono chrono;

            boost::shared_ptr<Ifpack_Preconditioner> precFluid;
            this->displayer().leaderPrint("  M-  Computing prec. factorization ...        ");
            chrono.start();
            createIfpackList(M_dataFile, "/solid/prec", M_listFluid);
            int overlapLevel = M_listFluid.get("overlap level", 2);
            std::string precType     = M_listFluid.get("prectype", "Amesos");
            Ifpack factory;
            precFluid.reset(factory.Create(precType, M_fluidOper->getMatrixPtr().get(), overlapLevel));

            if ( !precFluid.get() )
            {
                ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
            }

            precFluid->SetParameters(M_listFluid);
            precFluid->Initialize();
            precFluid->Compute();

            chrono.stop();
            this->displayer().leaderPrintMax("done in ", chrono.diff());

            boost::shared_ptr<Ifpack_Preconditioner> precSolid;
            this->displayer().leaderPrint("  M-  Computing prec. factorization ...        ");
            chrono.start();
            createIfpackList(M_dataFile, "/solid/prec",M_listSolid);
            overlapLevel = M_listSolid.get("overlap level", 2);
            precType     = M_listSolid.get("prectype", "Amesos");

            precSolid.reset(factory.Create(precType, M_solidOper->getMatrixPtr().get(), overlapLevel));
            if ( !precSolid.get() )
            {
                ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
            }
            precSolid->SetParameters( M_listSolid);
            precSolid->Initialize();
            precSolid->Compute();
            chrono.stop();
            this->displayer().leaderPrintMax("done in ", chrono.diff());


            //fluid part of the preconditioner
            boost::shared_ptr<Epetra_FECrsMatrix> matrCrsPtr1(new Epetra_FECrsMatrix(*M_couplingFluid1->getMatrixPtr()));
            boost::shared_ptr<Epetra_FECrsMatrix> matrCrsPtr2(new Epetra_FECrsMatrix(*M_couplingSolid1->getMatrixPtr()));

                firstCompPrec->push_back(precFluid, false, false);
                firstCompPrec->push_back(precSolid, false, false);
                matrCrsPtr1.reset (new Epetra_FECrsMatrix(*M_couplingFluid2->getMatrixPtr()));
                matrCrsPtr2.reset (new Epetra_FECrsMatrix(*M_couplingSolid2->getMatrixPtr()));


                secondCompPrec->push_back(boost::dynamic_pointer_cast<Epetra_Operator>(matrCrsPtr1), true, false/*already inverted*/);
                secondCompPrec->push_back(precFluid, false, false);
                secondCompPrec->push_back(boost::dynamic_pointer_cast<Epetra_Operator>(matrCrsPtr2), true, false/*already inverted*/);
                secondCompPrec->push_back(precSolid, false, false);

                M_CompPrec->push_back(firstCompPrec, false, false, false);
            }
                break;
        case 14:
            {
            firstCompPrec.reset(new ComposedPreconditioner<Epetra_Operator>(M_epetraComm.get()));
            M_CompPrec.reset(new ComposedPreconditioner<Epetra_Operator>);
            secondCompPrec.reset(new ComposedPreconditioner<Epetra_Operator>(M_epetraComm.get()));

            M_BCh_Robin->setOffset(M_offset);

            M_couplingFluid1.reset(new matrix_type(*M_monolithicMap, 1));
            M_couplingFluid2.reset(new matrix_type(*M_monolithicMap, 1));
            M_couplingSolid1.reset(new matrix_type(*M_monolithicMap, 1));
            M_couplingSolid2.reset(new matrix_type(*M_monolithicMap, 1));

            bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
            M_solidBlockPrec->GlobalAssemble();
            bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );
            M_solidBlockPrec->GlobalAssemble();

            if ( !M_BCh_u->bdUpdateDone() )
                M_BCh_u->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataFluid().dataTime()->getTime() );


            // fluid part of the preconditioner
            matrix_ptrtype tmpMatrix(new matrix_type(*M_monolithicMap, 0));
            couplingMatrix(tmpMatrix, 4);
            tmpMatrix->GlobalAssemble();
            *tmpMatrix *= 0.5;
            tmpMatrix->spy("rightup");
            *M_couplingFluid1 += *tmpMatrix;
            couplingMatrix(M_couplingFluid1, 2);
            *M_couplingFluid1 += *M_fluidBlock;
            M_couplingFluid1->GlobalAssemble();


            tmpMatrix.reset(new matrix_type(*M_monolithicMap, 0));
            couplingMatrix(tmpMatrix, 8);
            tmpMatrix->GlobalAssemble();
            *tmpMatrix *= 0.5;
            tmpMatrix->spy("rightlow");
            addDiagonalEntries(1., M_couplingSolid1, M_monolithicInterfaceMap, M_solidAndFluidDim);
            *M_couplingSolid1 += *tmpMatrix;
            *M_couplingSolid1 += *M_solidBlockPrec;
            M_couplingSolid1->GlobalAssemble();

            // solid part of the preconditioner
            tmpMatrix.reset(new matrix_type(*M_monolithicMap, 0));
            couplingMatrix(tmpMatrix, 8);
            tmpMatrix->GlobalAssemble();
            *tmpMatrix *= 0.5;
            tmpMatrix->spy("rightlow");
            couplingMatrix(M_couplingSolid2, 1);
            *M_couplingSolid2 += *tmpMatrix;
            *M_couplingSolid2 += *M_solidBlockPrec;

            M_couplingSolid2->GlobalAssemble();

            tmpMatrix.reset(new matrix_type(*M_monolithicMap, 0));
            couplingMatrix(tmpMatrix, 4);
            tmpMatrix->GlobalAssemble();
            *tmpMatrix *= 0.5;
            tmpMatrix->spy("rightup");
            addDiagonalEntries(1., M_couplingFluid2, M_monolithicInterfaceMap, M_solidAndFluidDim);
            *M_couplingFluid2 += *tmpMatrix;
            *M_couplingFluid2 += *M_fluidBlock;
            M_couplingFluid2->GlobalAssemble();


            Chrono chrono;

            boost::shared_ptr<Ifpack_Preconditioner> precFluid1;
            this->displayer().leaderPrint("  M-  Computing prec. factorization ...        ");
            chrono.start();
            createIfpackList(M_dataFile, "/solid/prec", M_listFluid);
            int overlapLevel = M_listFluid.get("overlap level", 4);
            std::string precType     = M_listFluid.get("prectype", "Amesos");
            Ifpack factory;
            precFluid1.reset(factory.Create(precType, M_couplingFluid1->getMatrixPtr().get(), overlapLevel));
            if ( !precFluid1.get() )
            {
                ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
            }
            IFPACK_CHK_ERR(precFluid1->SetParameters(M_listFluid));
            IFPACK_CHK_ERR(precFluid1->Initialize());
            IFPACK_CHK_ERR(precFluid1->Compute());
            chrono.stop();
            this->displayer().leaderPrintMax("done in ", chrono.diff());

            boost::shared_ptr<Ifpack_Preconditioner> precFluid2;
            this->displayer().leaderPrint("  M-  Computing prec. factorization ...        ");
            chrono.start();
            createIfpackList(M_dataFile, "/solid/prec", M_listFluid);
            overlapLevel = M_listFluid.get("overlap level", 4);
            precType     = M_listFluid.get("prectype", "Amesos");
            precFluid2.reset(factory.Create(precType, M_couplingFluid2->getMatrixPtr().get(), overlapLevel));
            if ( !precFluid2.get() )
            {
                ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
            }
            IFPACK_CHK_ERR(precFluid2->SetParameters(M_listFluid));
            IFPACK_CHK_ERR(precFluid2->Initialize());
            IFPACK_CHK_ERR(precFluid2->Compute());
            chrono.stop();
            this->displayer().leaderPrint(" overlap: ", M_listFluid.get("overlap level", 4));
            this->displayer().leaderPrintMax(" done in ", chrono.diff());


            boost::shared_ptr<Ifpack_Preconditioner> precSolid1;
            this->displayer().leaderPrint("  M-  Computing prec. factorization ...        ");
            chrono.start();
            createIfpackList(M_dataFile, "/solid/prec",M_listSolid);
            overlapLevel = M_listSolid.get("overlap level", 2);
            precType     = M_listSolid.get("prectype", "Amesos");
            precSolid1.reset(factory.Create(precType, M_couplingSolid1->getMatrixPtr().get(), overlapLevel));
            if ( !precSolid1.get() )
            {
                ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
            }
            IFPACK_CHK_ERR(precSolid1->SetParameters( M_listSolid));
            IFPACK_CHK_ERR(precSolid1->Initialize());
            IFPACK_CHK_ERR(precSolid1->Compute());
            chrono.stop();
            this->displayer().leaderPrintMax("done in ", chrono.diff());

            boost::shared_ptr<Ifpack_Preconditioner> precSolid2;
            this->displayer().leaderPrint("  M-  Computing prec. factorization ...        ");
            chrono.start();
            //Teuchos::ParameterList listSolid;
            createIfpackList(M_dataFile, "/solid/prec",M_listSolid);
            overlapLevel = M_listSolid.get("overlap level", 2);
            precType     = M_listSolid.get("prectype", "Amesos");
            precSolid2.reset(factory.Create(precType, M_couplingSolid2->getMatrixPtr().get(), overlapLevel));
            if ( !precSolid2.get() )
            {
                ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
            }
            IFPACK_CHK_ERR(precSolid2->SetParameters( M_listSolid));
            IFPACK_CHK_ERR(precSolid2->Initialize());
            IFPACK_CHK_ERR(precSolid2->Compute());
            chrono.stop();
            this->displayer().leaderPrintMax("done in ", chrono.diff());


            firstCompPrec->push_back(precFluid1, false, false);
            firstCompPrec->push_back(precSolid1, false, false);

            secondCompPrec->push_back(precFluid2, false, false);
            secondCompPrec->push_back(precSolid2, false, false);

            M_CompPrec->push_back(firstCompPrec, false, false, false);
            M_CompPrec->push_back(secondCompPrec, false, false, true /*sum*/);
        }
            break;

        default:
            {
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
    //     if(!M_data->dataFluid()->useShapeDerivatives())
    // #endif

    //M_precMatrPtr->spy("p");

    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _res.NormInf(), "\n");
    M_solid->getDisplayer().leaderPrint("  M-  Solving Jacobian system ...              \n" );

    switch(M_data->DDBlockPreconditioner())
    {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
        M_linearSolver->buildPreconditioner(M_precMatrPtr);
        this->iterateMonolithic(_res, _step, M_precMatrPtr->getMatrixPtr(), M_linearSolver);
        break;

    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        this->iterateMonolithic(_res, _step, M_precPtr->getPrecPtr(),     M_linearSolver);
        break;

    case 13:
    case 14:
        this->iterateMonolithic(_res, _step, M_CompPrec,     M_linearSolver);
        break;

    default:
        if(!M_precMatrPtr.get())
        {
            displayer().leaderPrint("Preconditioner type not implemented," );
            displayer().leaderPrint("change the entry DDBlockPrec in the data file");
            throw WRONG_PREC_EXCEPTION();
        }
        else
        {
            displayer().leaderPrint("default (Amesos) preconditioner choice" );
        }
        this->iterateMonolithic(_res, _step, M_precMatrPtr->getMatrixPtr(), M_linearSolver);
        break;
    }

    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _step.NormInf(), "\n");
}

void
Monolithic::iterateMesh(const vector_type& disp)
{
    vector_type lambdaFluid(this->M_interfaceMap, Unique);

    monolithicToInterface(lambdaFluid, disp);

    lambdaFluid *= (M_data->dataFluid()->dataTime()->getTimeStep()*(M_solid->rescaleFactor()));

    this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

    M_meshMotion->iterate(*M_BCh_mesh);

}

//void Monolithic::variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
void Monolithic::variablesInit(const std::string& dOrder)
{
    //    EpetraMap interfaceMap(*M_solidInterfaceMap);
    M_solidMeshPart.reset( new  partitionMesh< FSIOperator::mesh_type > (*M_data->dataSolid()->dataMesh()->mesh(), *M_epetraComm/*, M_solidInterfaceMap->getMap(Unique).get(), M_solidInterfaceMap->getMap(Repeated).get()*/));

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
//     M_bdf.reset(new BdfT<vector_type>(M_data->dataFluid()->getBDF_order()));
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
    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_data->dataSolid()->dataMesh()->mesh(),
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
    M_uFESpace->interpolate(u0, u, M_data->dataFluid()->dataTime()->getTime());

    vector_type p(M_pFESpace->map());
    M_pFESpace->interpolate(p0, p, M_data->dataFluid()->dataTime()->getTime());

    vector_type d(M_dFESpace->map());
    M_dFESpace->interpolate(d0, d, M_data->dataSolid()->dataTime()->getTime());

    //vector_type w(M_dFESpace->map());
    //M_dFESpace->interpolate(w0, w, M_data->dataSolid()->getTime());

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
Real& Monolithic::computeMaxSingularValue( )
{
    typedef Epetra_Operator                                                operator_type;

    M_PAAP.reset(new ComposedPreconditioner<Epetra_Operator>(M_epetraComm.get()));

    boost::shared_ptr<operator_type>  ComposedPrecPtr(M_linearSolver->getPrec()->getPrec());

    //M_monolithicMatrix->getMatrixPtr()->OptimizeStorage();
    boost::shared_ptr<Epetra_FECrsMatrix> matrCrsPtr(new Epetra_FECrsMatrix(*M_monolithicMatrix->getMatrixPtr()));

    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(/*ComposedPrecPtr*/matrCrsPtr));
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(ComposedPrecPtr/*matrCrsPtr*/),  true);
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(ComposedPrecPtr/*matrCrsPtr*/), true, true);
    M_PAAP->push_back(boost::dynamic_pointer_cast<operator_type>(/*ComposedPrecPtr*/matrCrsPtr), false, true);

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
    return real[0];
}
#endif


void Monolithic::computeFNormals( vector_type& normals)
{
    BCNormalManager<mesh_type, matrix_type> normalManager(*M_uFESpace->mesh());
    if ( !M_BChWSS->bdUpdateDone() )//possibly to avoid
        M_BChWSS->bdUpdate(*M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    normalManager.init((*M_BChWSS)[0], 0.);
    normalManager.computeIntegratedNormals(M_uFESpace->dof(), M_uFESpace->feBd(), normals, *M_uFESpace->mesh());
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
