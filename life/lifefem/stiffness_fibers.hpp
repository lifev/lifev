namespace LifeV
{

template <class vector_type>
void stiff(const Real sigma_l, const Real sigma_t, const vector_type& cos, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock);

template <class reduced_sigma, class vector_type>
void stiff( const reduced_sigma& red_sigma, const Real sigma_l, const Real sigma_t, const vector_type& cos, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock, ID id=0);

template <class reduced_sigma>
void stiff( reduced_sigma red_sigma, const Real D, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock, ID id=0);




template <class vector_type>
void stiff( const Real sigma_l, const Real sigma_t, const vector_type& cos, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock)
{
	//! Assembling the righthand side
    //ASSERT_PRE( fe.hasFirstDeriv(),
    //            "Stiffness matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int iloc, jloc;
    int i, icoor, jcoor, ig;
    Real s;
    ID eleId=fe.currentLocalId();
    UInt dim = dof.numTotalDof();

    Vector a_l(fe.nbCoor());
    Vector u_x(fe.nbQuadPt());
    Vector u_y(fe.nbQuadPt());
    Vector u_z(fe.nbQuadPt());

    for ( ig = 0;ig < fe.nbQuadPt();ig++ ){
        u_x[ig] = u_y[ig] = u_z[ig] = 0;
        for (i=0;i<fe.nbFENode();i++){
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
    int iloc, jloc;
    int i, icoor, jcoor, ig;
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
        for (i=0;i<fe.nbFENode();i++)
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
void stiff( reduced_sigma red_sigma, const Real D, ElemMat& elmat, const CurrentFE& fe, const Dof& dof, UInt iblock, UInt jblock, ID id)
{

//     ASSERT_PRE( fe.hasFirstDeriv(),
//                 "Stiffness matrix needs at least the first derivatives" );
    ElemMat::matrix_view mat = elmat.block( iblock, jblock );
    int iloc, jloc;
    int i, icoor, ig;
    double s, coef_s;
    int iu;
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


}
