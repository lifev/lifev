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
  @brief Contains an extension to the elemOper for including
  @a the fibers in the stiffness term for the Bidomain, Monodomain, and
  @ Electromechanical problems

  @date 12-2009
  @author Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>

  @contributors Simone Rossi <simone.rossi@epfl.ch>
  @mantainer Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>

 */




namespace LifeV
{

template <class vector_type>
void stiff(const Real sigma_l, const Real sigma_t, const vector_type& cos, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock);

template <class reduced_sigma, class vector_type>
void stiff( const reduced_sigma& red_sigma, const Real sigma_l, const Real sigma_t, const vector_type& cos, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock, ID id=0);

template <class reduced_sigma>
void stiff( reduced_sigma red_sigma, const Real D, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock, ID id=0);

template <class vector_type>
void stiffNL(vector_type& U, Real coef, ElemMat& elmat, const CurrentFE& fe,
             const Dof& dof, UInt iblock, UInt jblock, const Real beta);

//template <class vector_type>
//void stiffNL(vector_type& U, Real coef, ElemMat& elmat, const CurrentFE& fe,
//             const Dof& dof, UInt iblock, UInt jblock, UInt nb, const Real beta);

template <class vector_type>
void stiffNL(const vector_type& U, const Real sigma_l, const Real sigma_t,
                 const vector_type& cos, ElemMat& elmat, const CurrentFE& fe,
             const Dof& dof, UInt iblock, UInt jblock, const Real beta);

template <class vector_type>
void stiff( const Real sigma_l, const Real sigma_t, const vector_type& cos, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock)
{
	//! Assembling the righthand side
    //ASSERT_PRE( fe.hasFirstDeriv(),
    //            "Stiffness matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Int iloc, jloc;
    UInt i, icoor, jcoor, ig;
    Real s;
    ID eleId=fe.currentLocalId();
    UInt dim = dof.numTotalDof();

    Vector a_l(fe.nbCoor());
    Vector u_x(fe.nbQuadPt());
    Vector u_y(fe.nbQuadPt());
    Vector u_z(fe.nbQuadPt());

    for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
        u_x[ig] = u_y[ig] = u_z[ig] = 0;
        for (i=0;i<fe.nbFEDof();i++){
           u_x[ig]+=cos[dof.localToGlobal(eleId,i+1)]*fe.phi(i,ig);    //(one component)
           u_y[ig]+=cos[dof.localToGlobal(eleId,i+1)+dim]*fe.phi(i,ig);
           u_z[ig]+=cos[dof.localToGlobal(eleId,i+1)+2*dim]*fe.phi(i,ig);
        }
    }
    //
    // diagonal
    //

    for ( i = 0;i < fe.nbDiag();i++ )
       {
           iloc = fe.patternFirst( i );
           s = 0;
           for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
               a_l[0] = u_x[ig];
               a_l[1] = u_y[ig];
               a_l[2] = u_z[ig];
               Real norm = sqrt(a_l[0]*a_l[0]+a_l[1]*a_l[1]+a_l[2]*a_l[2]);
               a_l[0] = a_l[0]/norm; a_l[1] = a_l[1]/norm; a_l[2] = a_l[2]/norm;

             //  std::cout<< a_l[0] << a_l[1] << a_l[2] << " ";            //  D = sigma_t * I + (sigma_l-sigma_t) * a_l * a_l^T
               for ( icoor = 0;icoor < fe.nbCoor();icoor++ ){
                   s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig ) *
                           fe.weightDet( ig )* sigma_t;
                   for ( jcoor = 0;jcoor < fe.nbCoor(); jcoor++ )
                       s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, jcoor, ig ) *
                           fe.weightDet( ig )* (sigma_l-sigma_t)*a_l[icoor] * a_l[jcoor];
                    }

           }
           mat( iloc, iloc ) += s;
       }



    //
    // extra diagonal
    //
    for ( i = fe.nbDiag();i < fe.nbDiag() + fe.nbUpper();i++ )
    {
        s = 0;
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
            a_l[0] = u_x[ig];
            a_l[1] = u_y[ig];
            a_l[2] = u_z[ig];
            Real norm = sqrt(a_l[0]*a_l[0]+a_l[1]*a_l[1]+a_l[2]*a_l[2]);
            a_l[0] = a_l[0]/norm; a_l[1] = a_l[1]/norm; a_l[2] = a_l[2]/norm;
            //  D = sigma_t * I + (sigma_l-sigma_t) * a_l * a_l^T
            for ( icoor = 0;icoor < fe.nbCoor(); icoor++ ){
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                        fe.weightDet( ig )* sigma_t;  //diagonal
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                    s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, jcoor, ig ) *
                        fe.weightDet( ig )* (sigma_l-sigma_t)*a_l[icoor] * a_l[jcoor];
                    }
            }

        mat( iloc, jloc ) += s;
        mat( jloc, iloc ) += s;
    }
}


template <class reduced_sigma, class vector_type>
void stiff( const reduced_sigma& red_sigma, const Real sigma_l, const Real sigma_t, const vector_type& cos, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock, ID id)
{
//     ASSERT_PRE( fe.hasFirstDeriv(),
//                 "Stiffness matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Int iloc, jloc;
    UInt i, icoor, jcoor, ig;
    Real s;
    ID eleId=fe.currentLocalId();
    Real x,y,z;


    Vector a_l(fe.nbCoor());
    UInt dim = dof.numTotalDof();
    Vector u_x(fe.nbQuadPt());
    Vector u_y(fe.nbQuadPt());
    Vector u_z(fe.nbQuadPt());

    for ( ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        u_x[ig] = u_y[ig] = u_z[ig] = 0;
        for (i=0;i<fe.nbFEDof();i++)
	{
           u_x[ig]+=cos[dof.localToGlobal(eleId,i+1)]*fe.phi(i,ig);    //(one component)
           u_y[ig]+=cos[dim+dof.localToGlobal(eleId,i+1)]*fe.phi(i,ig);
           u_z[ig]+=cos[2*dim+dof.localToGlobal(eleId,i+1)]*fe.phi(i,ig);
        }
    }

    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag();i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
       for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
            a_l[0] = u_x[ig];
            a_l[1] = u_y[ig];
            a_l[2] = u_z[ig];
            Real norm = sqrt(a_l[0]*a_l[0]+a_l[1]*a_l[1]+a_l[2]*a_l[2]);
            a_l[0] = a_l[0]/norm; a_l[1] = a_l[1]/norm; a_l[2] = a_l[2]/norm;

            fe.coorQuadPt(x,y,z,ig);
           //  D = sigma_t * I + (sigma_l-sigma_t) * a_l * a_l^T
            for ( icoor = 0;icoor < fe.nbCoor();icoor++ )
	    {
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig ) *
                        fe.weightDet( ig )* red_sigma(x,y,z,0,id,sigma_t);
                for ( jcoor = 0;jcoor < fe.nbCoor(); jcoor++ )
                    s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, jcoor, ig ) *
                        fe.weightDet( ig )* (red_sigma(x,y,z,0,id, sigma_l)-red_sigma(x,y,z,0,id, sigma_t))*a_l[icoor] * a_l[jcoor];
                 }

        }
        mat( iloc, iloc ) += s;
    }
    //
    // extra diagonal
    //
    for ( i = fe.nbDiag();i < fe.nbDiag() + fe.nbUpper();i++ )
    {
        s = 0;
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        for ( ig = 0;ig < fe.nbQuadPt();ig++ )
	{
            a_l[0] = u_x[ig];
            a_l[1] = u_y[ig];
            a_l[2] = u_z[ig];

            Real norm = sqrt(a_l[0]*a_l[0]+a_l[1]*a_l[1]+a_l[2]*a_l[2]);
            a_l[0] = a_l[0]/norm; a_l[1] = a_l[1]/norm; a_l[2] = a_l[2]/norm;

            //  D = sigma_t * I + (sigma_l-sigma_t) * a_l * a_l^T
            for ( icoor = 0;icoor < fe.nbCoor(); icoor++ )
	    {
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                        fe.weightDet( ig )* red_sigma(x,y,z,0,id,sigma_t);  //diagonal
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                    s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, jcoor, ig ) *
                        fe.weightDet( ig )* (red_sigma(x,y,z,0,id,sigma_l)-red_sigma(x,y,z,0,id,sigma_t))*a_l[icoor] * a_l[jcoor];
            }
         }

        mat( iloc, jloc ) += s;
        mat( jloc, iloc ) += s;
    }
}


template <class reduced_sigma>
void stiff( reduced_sigma red_sigma, const Real D, ElemMat& elmat, const CurrentFE& fe, const Dof& /*dof*/, UInt iblock, UInt jblock, ID id)
{

//     ASSERT_PRE( fe.hasFirstDeriv(),
//                 "Stiffness matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Int iloc, jloc;
    UInt i, icoor, ig;
    Real s;

    Real x,y,z;
    for ( i = 0;i < fe.nbDiag();i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt();ig++ )
        {
            fe.coorQuadPt(x,y,z,ig);
            for ( icoor = 0;icoor < fe.nbCoor();icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig )
                    * fe.weightDet( ig )*red_sigma(x,y,z,0,id, D);
        }
        mat( iloc, iloc ) += s;
    }


    //
    // extra diagonal
    //


    for ( i = fe.nbDiag();i < fe.nbDiag() + fe.nbUpper();i++ )
    {
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt();ig++ )
        {
            fe.coorQuadPt(x,y,z,ig);
            for ( icoor = 0;icoor < fe.nbCoor();icoor++ )
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                    fe.weightDet( ig )*red_sigma(x,y,z,0,id,D);
        }
        mat( iloc, jloc ) += s;
        mat( jloc, iloc ) += s;
    }
}

//Without fibers
template <class vector_type>
void stiffNL(vector_type& U, Real coef, ElemMat& elmat, const CurrentFE& fe,
             const Dof& dof, UInt iblock, UInt jblock, const Real beta)
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Int iloc, jloc, i, icoor, jcoor, ig;
    Real s, coef_s;
    ID eleId=fe.currentLocalId();
    UInt dim = dof.numTotalDof();
    Int iu;
    std::vector<Real> locU(fe.nbFEDof());
    std::vector<Real> MM_l(fe.nbCoor());
    Real bPt;

    for (i=0;i<fe.nbFEDof();i++){
        locU[i]=U[dof.localToGlobal(eleId,i+1)];    //(one component)
    }

    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag();i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt();ig++ )
        {
            bPt = 0.0;
            for(iu=0;iu<fe.nbFEDof();iu++){
                bPt+=locU[iu]*fe.phi(iu,ig);
            }
            MM_l[0] = 1.0/(1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt);
            MM_l[1] = 1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;
            MM_l[2] = 1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;

            for ( icoor = 0;icoor < fe.nbCoor();icoor++ ){
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig ) *
                    fe.weightDet( ig )*MM_l[icoor];
            }
        }
        mat( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //

    for ( i = fe.nbDiag();i < fe.nbDiag() + fe.nbUpper();i++ )
    {
        s = 0;
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
            bPt = 0.0;
            for(iu=0;iu<fe.nbFEDof();iu++){
                bPt+=locU[iu]*fe.phi(iu,ig);
            }
            MM_l[0] = 1.0/(1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt);
            MM_l[1]=1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;
            MM_l[2]=1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;

            for ( icoor = 0;icoor < fe.nbCoor();icoor++ ){
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                    fe.weightDet( ig )*MM_l[icoor];
            }
        }
        coef_s= coef * s;
        mat( iloc, jloc ) += coef_s;
        mat( jloc, iloc ) += coef_s;
    }
}


//----------------------------------------------------------------------------

/* This one is for the vectorial case, only 'UInt nb' has been added to the argument
template <class vector_type>
void stiffNL(vector_type& U, Real coef, ElemMat& elmat, const CurrentFE& fe,
             const Dof& dof, UInt iblock, UInt jblock, UInt nb, const Real beta)
{
    Tab2d mat_tmp( fe.nbFEDof(), fe.nbFEDof() );
    mat_tmp = ZeroMatrix( fe.nbFEDof(), fe.nbFEDof() );

    Int iloc, jloc, i, icoor, jcoor, ig;
    Real s, coef_s;

    ID eleId=fe.currentLocalId();
    UInt dim = dof.numTotalDof();
    Int iu;
    std::vector<Real> locU(fe.nbFEDof());
    std::vector<Real> MM_l(fe.nbCoor());
    Real bPt;

    for (i=0;i<fe.nbFEDof();i++){
        locU[i]=U[dof.localToGlobal(eleId,i+1)];    //(one component)
    }

    //
    // diagonal
    //
    for ( i = 0;i < fe.nbDiag();i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt();ig++ )
        {
            bPt = 0.0;
            for(iu=0;iu<fe.nbFEDof();iu++){
                bPt+=locU[iu]*fe.phi(iu,ig);
            }
            MM_l[0] = 1.0/(1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt);
            MM_l[1]=1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;
            MM_l[2]=1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;

            for ( icoor = 0;icoor < fe.nbCoor();icoor++ ){
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig ) *
                    fe.weightDet( ig )*MM_l[icoor];
            }
        }
        mat_tmp( iloc, iloc ) += coef * s;
    }
    //
    // extra diagonal
    //

    for ( i = fe.nbDiag();i < fe.nbDiag() + fe.nbUpper();i++ )
    {
        s = 0;
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
            bPt = 0.0;
            for(iu=0;iu<fe.nbFEDof();iu++){
                bPt+=locU[iu]*fe.phi(iu,ig);
            }
            MM_l[0] = 1.0/(1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt);
            MM_l[1]=1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;
            MM_l[2]=1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;

            for ( icoor = 0;icoor < fe.nbCoor();icoor++ ){
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                    fe.weightDet( ig )*MM_l[icoor];
            }
        }
        coef_s= coef * s;
        mat_tmp( iloc, jloc ) += coef_s;
        mat_tmp( jloc, iloc ) += coef_s;
    }
    // copy on the components
    for ( Int icomp = 0;icomp < nb;icomp++ )
    {ElemMat::matrix_view mat_icomp = elmat.block( iblock + icomp, jblock + icomp );
        for ( i = 0;i < fe.nbDiag();i++ )
        {
            iloc = fe.patternFirst( i );
            mat_icomp( iloc, iloc ) += mat_tmp( iloc, iloc );
        }
        for ( i = fe.nbDiag();i < fe.nbDiag() + fe.nbUpper();i++ )
        {
            iloc = fe.patternFirst( i );
            jloc = fe.patternSecond( i );
            mat_icomp( iloc, jloc ) += mat_tmp( iloc, jloc );
            mat_icomp( jloc, iloc ) += mat_tmp( jloc, iloc );
        }
    }
}

*/

//----------------------------------------------------------------------------



//with the fibers
template <class vector_type>
void stiffNL(const vector_type& U, const Real sigma_l, const Real sigma_t,
             const vector_type& cos, ElemMat& elmat, const CurrentFE& fe,
             const Dof& dof, UInt iblock, UInt jblock, const Real beta)
{
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    Int iloc, jloc;
    UInt i, icoor, jcoor, ig;
    Real s;
    ID eleId=fe.currentLocalId();
    UInt dim = dof.numTotalDof();
    Int iu;
    Vector locU(fe.nbFEDof());
    Vector MM_l(fe.nbCoor());
    Real bPt;
    Vector a_l(fe.nbCoor());
    Vector u_x(fe.nbQuadPt());
    Vector u_y(fe.nbQuadPt());
    Vector u_z(fe.nbQuadPt());

    for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
        u_x[ig] = u_y[ig] = u_z[ig] = 0;
        for (i=0;i<fe.nbFEDof();i++){
            locU[i]=U[dof.localToGlobal(eleId,i+1)]; // ojo!!! quizas es -1
            u_x[ig]+=cos[dof.localToGlobal(eleId,i+1)]*fe.phi(i,ig);    //(one component)
            u_y[ig]+=cos[dof.localToGlobal(eleId,i+1)+dim]*fe.phi(i,ig);
            u_z[ig]+=cos[dof.localToGlobal(eleId,i+1)+2*dim]*fe.phi(i,ig);
        }
    }
    //
    // diagonal
    //

    for ( i = 0;i < fe.nbDiag();i++ )
    {
        iloc = fe.patternFirst( i );
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
            bPt = 0.0;
            for(iu=0;iu<fe.nbFEDof();iu++){
                bPt+=locU[iu]*fe.phi(iu,ig);
            }
            MM_l[0]= 1.0/(1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt);
            MM_l[1]= 1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;
            MM_l[2]= 1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;
            a_l[0] = u_x[ig];
            a_l[1] = u_y[ig];
            a_l[2] = u_z[ig];
            Real norm = sqrt(a_l[0]*a_l[0]+a_l[1]*a_l[1]+a_l[2]*a_l[2]);
            a_l[0] = a_l[0]/norm; a_l[1] = a_l[1]/norm; a_l[2] = a_l[2]/norm;

            //  D = sigma_t * I + (sigma_l-sigma_t) * a_l * a_l^T
            for ( icoor = 0;icoor < fe.nbCoor();icoor++ ){
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, icoor, ig ) *
                    fe.weightDet( ig )* sigma_t*MM_l[icoor];
                for ( jcoor = 0;jcoor < fe.nbCoor(); jcoor++ )
                    s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( iloc, jcoor, ig ) *
                        fe.weightDet( ig )* (sigma_l-sigma_t)*a_l[icoor] * a_l[jcoor];
            }

        }
        mat( iloc, iloc ) += s;
    }



    //
    // extra diagonal
    //
    for ( i = fe.nbDiag();i < fe.nbDiag() + fe.nbUpper();i++ )
    {
        s = 0;
        iloc = fe.patternFirst( i );
        jloc = fe.patternSecond( i );
        for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
            bPt = 0.0;
            for(iu=0;iu<fe.nbFEDof();iu++){
                bPt+=locU[iu]*fe.phi(iu,ig);
            }
            MM_l[0]= 1.0/(1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt);
            MM_l[1]= 1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;
        MM_l[2]= 1.0+beta*(bPt+84.0)/(184.0+bPt)+0.001*bPt;
            a_l[0] = u_x[ig];
            a_l[1] = u_y[ig];
            a_l[2] = u_z[ig];
            Real norm = sqrt(a_l[0]*a_l[0]+a_l[1]*a_l[1]+a_l[2]*a_l[2]);
            a_l[0] = a_l[0]/norm; a_l[1] = a_l[1]/norm; a_l[2] = a_l[2]/norm;
            //  D = sigma_t * I + (sigma_l-sigma_t) * a_l * a_l^T
            for ( icoor = 0;icoor < fe.nbCoor(); icoor++ ){
                s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, icoor, ig ) *
                    fe.weightDet( ig )* sigma_t*MM_l[icoor];  //diagonal
                for ( jcoor = 0; jcoor < fe.nbCoor(); jcoor++ )
                    s += fe.phiDer( iloc, icoor, ig ) * fe.phiDer( jloc, jcoor, ig ) *
                        fe.weightDet( ig )* (sigma_l-sigma_t)*a_l[icoor] * a_l[jcoor];
            }
        }

        mat( iloc, jloc ) += s;
        mat( jloc, iloc ) += s;
    }
}


/*-------------------------------------------------------------------------------*/



}
