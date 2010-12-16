//  cbessjy.cpp -- complex Bessel functions.
//  Algorithms and coefficient values from "Computation of Special
//  Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
#include <complex>
#include <testsuite/functions/bessel/bessel.hpp>
double gamma(double);
using namespace std;

static complex<double> cii(0.0,1.0);
static complex<double> cone(1.0,0.0);
static complex<double> czero(0.0,0.0);

int cbessjy01(complex<double> z,complex<double> &cj0,complex<double> &cj1,
    complex<double> &cy0,complex<double> &cy1,complex<double> &cj0p,
    complex<double> &cj1p,complex<double> &cy0p,complex<double> &cy1p)
{
    complex<double> z1,z2,cr,cp,cs,cp0,cq0,cp1,cq1,ct1,ct2,cu;
    double a0,w0,w1;
    int k,kz;

    static double a[] = {
        -7.03125e-2,
         0.112152099609375,
        -0.5725014209747314,
         6.074042001273483,
        -1.100171402692467e2,
         3.038090510922384e3,
        -1.188384262567832e5,
         6.252951493434797e6,
        -4.259392165047669e8,
         3.646840080706556e10,
        -3.833534661393944e12,
         4.854014686852901e14,
        -7.286857349377656e16,
         1.279721941975975e19};
    static double b[] = {
         7.32421875e-2,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1,
         5.513358961220206e2,
        -1.825775547429318e4,
         8.328593040162893e5,
        -5.006958953198893e7,
         3.836255180230433e9,
        -3.649010818849833e11,
         4.218971570284096e13,
        -5.827244631566907e15,
         9.476288099260110e17,
        -1.792162323051699e20};
    static double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947,
         1.215978918765359e2,
        -3.302272294480852e3,
         1.276412726461746e5,
        -6.656367718817688e6,
         4.502786003050393e8,
        -3.833857520742790e10,
         4.011838599133198e12,
        -5.060568503314727e14,
         7.572616461117958e16,
        -1.326257285320556e19};
    static double b1[] = {
        -0.1025390625,
         0.2775764465332031,
        -1.993531733751297,
         2.724882731126854e1,
        -6.038440767050702e2,
         1.971837591223663e4,
        -8.902978767070678e5,
         5.310411010968522e7,
        -4.043620325107754e9,
         3.827011346598605e11,
        -4.406481417852278e13,
         6.065091351222699e15,
        -9.833883876590679e17,
         1.855045211579828e20};

    a0 = abs(z);
    z2 = z*z;
    z1 = z;
    if (a0 == 0.0) {
        cj0 = cone;
        cj1 = czero;
        cy0 = complex<double>(-1e308,0);
        cy1 = complex<double>(-1e308,0);
        cj0p = czero;
        cj1p = complex<double>(0.5,0.0);
        cy0p = complex<double>(1e308,0);
        cy1p = complex<double>(1e308,0);
        return 0;
    }
    if (real(z) < 0.0) z1 = -z;
    if (a0 <= 12.0) {
        cj0 = cone;
        cr = cone;
        for (k=1;k<=40;k++) {
            cr *= -0.25*z2/(double)(k*k);
            cj0 += cr;
            if (abs(cr) < abs(cj0)*eps) break;
        }
        cj1 = cone;
        cr = cone;
        for (k=1;k<=40;k++) {
            cr *= -0.25*z2/(k*(k+1.0));
            cj1 += cr;
            if (abs(cr) < abs(cj1)*eps) break;
        }
        cj1 *= 0.5*z1;
        w0 = 0.0;
        cr = cone;
        cs = czero;
        for (k=1;k<=40;k++) {
            w0 += 1.0/k;
            cr *= -0.25*z2/(double)(k*k);
            cp = cr*w0;
            cs += cp;
            if (abs(cp) < abs(cs)*eps) break;
        }
        cy0 = M_2_PI*((log(0.5*z1)+el)*cj0-cs);
        w1 = 0.0;
        cr = cone;
        cs = cone;
        for (k=1;k<=40;k++) {
            w1 += 1.0/k;
            cr *= -0.25*z2/(k*(k+1.0));
            cp = cr*(2.0*w1+1.0/(k+1.0));
            cs += cp;
            if (abs(cp) < abs(cs)*eps) break;
        }
        cy1 = M_2_PI*((log(0.5*z1)+el)*cj1-1.0/z1-0.25*z1*cs);
    }
    else {
        if (a0 >= 50.0) kz = 8;         // can be changed to 10
        else if (a0 >= 35.0) kz = 10;   //   "      "     "  12
        else kz = 12;                   //   "      "     "  14
        ct1 = z1 - M_PI_4;
        cp0 = cone;
        for (k=0;k<kz;k++) {
            cp0 += a[k]*pow(z1,-2.0*k-2.0);
        }
        cq0 = -0.125/z1;
        for (k=0;k<kz;k++) {
            cq0 += b[k]*pow(z1,-2.0*k-3.0);
        }
        cu = sqrt(M_2_PI/z1);
        cj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1));
        cy0 = cu*(cp0*sin(ct1)+cq0*cos(ct1));
        ct2 = z1 - 0.75*M_PI;
        cp1 = cone;
        for (k=0;k<kz;k++) {
            cp1 += a1[k]*pow(z1,-2.0*k-2.0);
        }
        cq1 = 0.375/z1;
        for (k=0;k<kz;k++) {
            cq1 += b1[k]*pow(z1,-2.0*k-3.0);
        }
        cj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2));
        cy1 = cu*(cp1*sin(ct2)+cq1*cos(ct2));
    }
    if (real(z) < 0.0) {
        if (imag(z) < 0.0) {
            cy0 -= 2.0*cii*cj0;
            cy1 = -(cy1-2.0*cii*cj1);
        }
        else if (imag(z) > 0.0) {
            cy0 += 2.0*cii*cj0;
            cy1 = -(cy1+2.0*cii*cj1);
        }
        cj1 = -cj1;
    }
    cj0p = -cj1;
    cj1p = cj0-cj1/z;
    cy0p = -cy1;
    cy1p = cy0-cy1/z;
    return 0;
}

int cbessjyna(int n,complex<double> z,int &nm,complex<double> *cj,
    complex<double> *cy,complex<double> *cjp,complex<double> *cyp)
{
    complex<double> cbj0,cbj1,cby0,cby1,cj0,cjk,cj1,cf,cf1,cf2;
    complex<double> cs,cg0,cg1,cyk,cyl1,cyl2,cylk,cp11,cp12,cp21,cp22;
    complex<double> ch0,ch1,ch2;
    double a0,yak,ya1,ya0,wa;
    int m,k,lb,lb0;

    if (n < 0) return 1;
    a0 = abs(z);
    nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            cj[k] = czero;
            cy[k] = complex<double> (-1e308,0);
            cjp[k] = czero;
            cyp[k] = complex<double>(1e308,0);
        }
        cj[0] = cone;
        cjp[1] = complex<double>(0.5,0.0);
        return 0;
    }
    cbessjy01(z,cj[0],cj[1],cy[0],cy[1],cjp[0],cjp[1],cyp[0],cyp[1]);
    cbj0 = cj[0];
    cbj1 = cj[1];
    cby0 = cy[0];
    cby1 = cy[1];
    if (n <= 1) return 0;
    if (n < (int)0.25*a0) {
        cj0 = cbj0;
        cj1 = cbj1;
        for (k=2;k<=n;k++) {
            cjk = 2.0*(k-1.0)*cj1/z-cj0;
            cj[k] = cjk;
            cj0 = cj1;
            cj1 = cjk;
        }
    }
    else {
        m = msta1(a0,200);
        if (m < n) nm = m;
        else m = msta2(a0,n,15);
        cf2 = czero;
        cf1 = complex<double> (1.0e-100,0.0);
        for (k=m;k>=0;k--) {
            cf = 2.0*(k+1.0)*cf1/z-cf2;
            if (k <=nm) cj[k] = cf;
            cf2 = cf1;
            cf1 = cf;
        }
        if (abs(cbj0) > abs(cbj1)) cs = cbj0/cf;
        else cs = cbj1/cf2;
        for (k=0;k<=nm;k++) {
            cj[k] *= cs;
        }
    }
    for (k=2;k<=nm;k++) {
        cjp[k] = cj[k-1]-(double)k*cj[k]/z;
    }
    ya0 = abs(cby0);
    lb = 0;
    cg0 = cby0;
    cg1 = cby1;
    for (k=2;k<=nm;k++) {
        cyk = 2.0*(k-1.0)*cg1/z-cg0;
        yak = abs(cyk);
        ya1 = abs(cg0);
        if ((yak < ya0) && (yak < ya1)) lb = k;
        cy[k] = cyk;
        cg0 = cg1;
        cg1 = cyk;
    }
    lb0 = 0;
    if ((lb > 4) && (imag(z) != 0.0)) {
        while (lb != lb0) {
            ch2 = cone;
            ch1 = czero;
            lb0 = lb;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*k*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp12 = ch0;
            cp22 = ch2;
            ch2 = czero;
            ch1 = cone;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*k*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp11 = ch0;
            cp21 = ch2;
            if (lb == nm)
                cj[lb+1] = 2.0*lb*cj[lb]/z-cj[lb-1];
            if (abs(cj[0]) > abs(cj[1])) {
                cy[lb+1] = (cj[lb+1]*cby0-2.0*cp11/(M_PI*z))/cj[0];
                cy[lb] = (cj[lb]*cby0+2.0*cp12/(M_PI*z))/cj[0];
            }
            else {
                cy[lb+1] = (cj[lb+1]*cby1-2.0*cp21/(M_PI*z))/cj[1];
                cy[lb] = (cj[lb]*cby1+2.0*cp22/(M_PI*z))/cj[1];
            }
            cyl2 = cy[lb+1];
            cyl1 = cy[lb];
            for (k=lb-1;k>=0;k--) {
                cylk = 2.0*(k+1.0)*cyl1/z-cyl2;
                cy[k] = cylk;
                cyl2 = cyl1;
                cyl1 = cylk;
            }
            cyl1 = cy[lb];
            cyl2 = cy[lb+1];
            for (k=lb+1;k<n;k++) {
                cylk = 2.0*k*cyl2/z-cyl1;
                cy[k+1] = cylk;
                cyl1 = cyl2;
                cyl2 = cylk;
            }
            for (k=2;k<=nm;k++) {
                wa = abs(cy[k]);
                if (wa < abs(cy[k-1])) lb = k;
            }
        }
    }
    for (k=2;k<=nm;k++) {
        cyp[k] = cy[k-1]-(double)k*cy[k]/z;
    }
    return 0;
}

int cbessjynb(int n,complex<double> z,int &nm,complex<double> *cj,
    complex<double> *cy,complex<double> *cjp,complex<double> *cyp)
{
    complex<double> cf,cf0,cf1,cf2,cbs,csu,csv,cs0,ce;
    complex<double> ct1,cp0,cq0,cp1,cq1,cu,cbj0,cby0,cbj1,cby1;
    complex<double> cyy,cbjk,ct2;
    double a0,y0;
    int k,m;
    static double a[] = {
        -0.7031250000000000e-1,
         0.1121520996093750,
        -0.5725014209747314,
         6.074042001273483};
    static double b[] = {
         0.7324218750000000e-1,
        -0.2271080017089844,
         1.727727502584457,
        -2.438052969955606e1};
    static double a1[] = {
         0.1171875,
        -0.1441955566406250,
         0.6765925884246826,
        -6.883914268109947};
    static double b1[] = {
       -0.1025390625,
        0.2775764465332031,
       -1.993531733751297,
        2.724882731126854e1};

    y0 = abs(imag(z));
    a0 = abs(z);
    nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            cj[k] = czero;
            cy[k] = complex<double> (-1e308,0);
            cjp[k] = czero;
            cyp[k] = complex<double>(1e308,0);
        }
        cj[0] = cone;
        cjp[1] = complex<double>(0.5,0.0);
        return 0;
    }
    if ((a0 <= 300.0) || (n > (int)(0.25*a0))) {
        if (n == 0) nm = 1;
        m = msta1(a0,200);
        if (m < nm) nm = m;
        else m = msta2(a0,nm,15);
        cbs = czero;
        csu = czero;
        csv = czero;
        cf2 = czero;
        cf1 = complex<double> (1.0e-100,0.0);
        for (k=m;k>=0;k--) {
            cf = 2.0*(k+1.0)*cf1/z-cf2;
            if (k <= nm) cj[k] = cf;
            if (((k & 1) == 0) && (k != 0)) {
                if (y0 <= 1.0) {
                    cbs += 2.0*cf;
                }
                else {
                    cbs += (-1)*((k & 2)-1)*2.0*cf;
                }
                csu += (double)((-1)*((k & 2)-1))*cf/(double)k;
            }
            else if (k > 1) {
                csv += (double)((-1)*((k & 2)-1)*k)*cf/(double)(k*k-1.0);
            }
            cf2 = cf1;
            cf1 = cf;
        }
        if (y0 <= 1.0) cs0 = cbs+cf;
        else cs0 = (cbs+cf)/cos(z);
        for (k=0;k<=nm;k++) {
            cj[k] /= cs0;
        }
        ce = log(0.5*z)+el;
        cy[0] = M_2_PI*(ce*cj[0]-4.0*csu/cs0);
        cy[1] = M_2_PI*(-cj[0]/z+(ce-1.0)*cj[1]-4.0*csv/cs0);
    }
    else {
        ct1 = z-M_PI_4;
        cp0 = cone;
        for (k=0;k<4;k++) {
            cp0 += a[k]*pow(z,-2.0*k-2.0);
        }
        cq0 = -0.125/z;
        for (k=0;k<4;k++) {
            cq0 += b[k] *pow(z,-2.0*k-3.0);
        }
        cu = sqrt(M_2_PI/z);
        cbj0 = cu*(cp0*cos(ct1)-cq0*sin(ct1));
        cby0 = cu*(cp0*sin(ct1)+cq0*cos(ct1));
        cj[0] = cbj0;
        cy[0] = cby0;
        ct2 = z-0.75*M_PI;
        cp1 = cone;
        for (k=0;k<4;k++) {
            cp1 += a1[k]*pow(z,-2.0*k-2.0);
        }
        cq1 = 0.375/z;
        for (k=0;k<4;k++) {
            cq1 += b1[k]*pow(z,-2.0*k-3.0);
        }
        cbj1 = cu*(cp1*cos(ct2)-cq1*sin(ct2));
        cby1 = cu*(cp1*sin(ct2)+cq1*cos(ct2));
        cj[1] = cbj1;
        cy[1] = cby1;
        for (k=2;k<=n;k++) {
            cbjk = 2.0*(k-1.0)*cbj1/z-cbj0;
            cj[k] = cbjk;
            cbj0 = cbj1;
            cbj1 = cbjk;
        }
    }
    cjp[0] = -cj[1];
    for (k=1;k<=nm;k++) {
        cjp[k] = cj[k-1]-(double)k*cj[k]/z;
    }
    if (abs(cj[0]) > 1.0)
        cy[1] = (cj[1]*cy[0]-2.0/(M_PI*z))/cj[0];
    for (k=2;k<=nm;k++) {
        if (abs(cj[k-1]) >= abs(cj[k-2]))
            cyy = (cj[k]*cy[k-1]-2.0/(M_PI*z))/cj[k-1];
        else
            cyy = (cj[k]*cy[k-2]-4.0*(k-1.0)/(M_PI*z*z))/cj[k-2];
        cy[k] = cyy;
    }
    cyp[0] = -cy[1];
    for (k=1;k<=nm;k++) {
        cyp[k] = cy[k-1]-(double)k*cy[k]/z;
    }

    return 0;
}

int cbessjyva(double v,complex<double> z,double &vm,complex<double>*cjv,
    complex<double>*cyv,complex<double>*cjvp,complex<double>*cyvp)
{
    complex<double> z1,z2,zk,cjvl,cr,ca,cjv0,cjv1,cpz,crp;
    complex<double> cqz,crq,ca0,cck,csk,cyv0,cyv1,cju0,cju1,cb;
    complex<double> cs,cs0,cr0,cs1,cr1,cec,cf,cf0,cf1,cf2;
    complex<double> cfac0,cfac1,cg0,cg1,cyk,cp11,cp12,cp21,cp22;
    complex<double> ch0,ch1,ch2,cyl1,cyl2,cylk;

    double a0,v0,pv0,pv1,vl,ga,gb,vg,vv,w0,w1,ya0,yak,ya1,wa;
    int j,n,k,kz,l,lb,lb0,m;

    a0 = abs(z);
    z1 = z;
    z2 = z*z;
    n = (int)v;


    v0 = v-n;

    pv0 = M_PI*v0;
    pv1 = M_PI*(1.0+v0);
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            cjv[k] = czero;
            cyv[k] = complex<double> (-1e308,0);
            cjvp[k] = czero;
            cyvp[k] = complex<double> (1e308,0);

        }
        if (v0 == 0.0) {
            cjv[0] = cone;
            cjvp[1] = complex<double> (0.5,0.0);
        }
        else {
            cjvp[0] = complex<double> (1e308,0);
        }
        vm = v;
        return 0;
    }
    if (real(z1) < 0.0) z1 = -z;
    if (a0 <= 12.0) {
        for (l=0;l<2;l++) {
            vl = v0+l;
            cjvl = cone;
            cr = cone;
            for (k=1;k<=40;k++) {
                cr *= -0.25*z2/(k*(k+vl));
                cjvl += cr;
                if (abs(cr) < abs(cjvl)*eps) break;
            }
           vg = 1.0 + vl;
           ga = gamma(vg);
           ca = pow(0.5*z1,vl)/ga;
           if (l == 0) cjv0 = cjvl*ca;
           else cjv1 = cjvl*ca;
        }
    }
    else {
        if (a0 >= 50.0) kz = 8;
        else if (a0 >= 35.0) kz = 10;
        else kz = 11;
        for (j=0;j<2;j++) {
            vv = 4.0*(j+v0)*(j+v0);
            cpz = cone;
            crp = cone;
            for (k=1;k<=kz;k++) {
                crp = -0.78125e-2*crp*(vv-pow(4.0*k-3.0,2.0))*
                    (vv-pow(4.0*k-1.0,2.0))/(k*(2.0*k-1.0)*z2);
                cpz += crp;
            }
            cqz = cone;
            crq = cone;
            for (k=1;k<=kz;k++) {
                crq = -0.78125e-2*crq*(vv-pow(4.0*k-1.0,2.0))*
                    (vv-pow(4.0*k+1.0,2.0))/(k*(2.0*k+1.0)*z2);
                cqz += crq;
            }
            cqz *= 0.125*(vv-1.0)/z1;
            zk = z1-(0.5*(j+v0)+0.25)*M_PI;
            ca0 = sqrt(M_2_PI/z1);
            cck = cos(zk);
            csk = sin(zk);
            if (j == 0) {
                cjv0 = ca0*(cpz*cck-cqz*csk);
                cyv0 = ca0*(cpz*csk+cqz+cck);
            }
            else {
                cjv1 = ca0*(cpz*cck-cqz*csk);
                cyv1 = ca0*(cpz*csk+cqz*cck);
            }
        }
    }
    if (a0 <= 12.0) {
        if (v0 != 0.0) {
            for (l=0;l<2;l++) {
                vl = v0+l;
                cjvl = cone;
                cr = cone;
                for (k=1;k<=40;k++) {
                    cr *= -0.25*z2/(k*(k-vl));
                    cjvl += cr;
                    if (abs(cr) < abs(cjvl)*eps) break;
                }
                vg = 1.0-vl;
                gb = gamma(vg);
                cb = pow(2.0/z1,vl)/gb;
                if (l == 0) cju0 = cjvl*cb;
                else cju1 = cjvl*cb;
            }
            cyv0 = (cjv0*cos(pv0)-cju0)/sin(pv0);
            cyv1 = (cjv1*cos(pv1)-cju1)/sin(pv1);
        }
        else {
            cec = log(0.5*z1)+el;
            cs0 = czero;
            w0 = 0.0;
            cr0 = cone;
            for (k=1;k<=30;k++) {
                w0 += 1.0/k;
                cr0 *= -0.25*z2/(double)(k*k);
                cs0 += cr0*w0;
            }
            cyv0 = M_2_PI*(cec*cjv0-cs0);
            cs1 = cone;
            w1 = 0.0;
            cr1 = cone;
            for (k=1;k<=30;k++) {
                w1 += 1.0/k;
                cr1 *= -0.25*z2/(k*(k+1.0));
                cs1 += cr1*(2.0*w1+1.0/(k+1.0));
            }
            cyv1 = M_2_PI*(cec*cjv1-1.0/z1-0.25*z1*cs1);
        }
    }
    if (real(z) < 0.0) {
        cfac0 = exp(pv0*cii);
        cfac1 = exp(pv1*cii);
        if (imag(z) < 0.0) {
            cyv0 = cfac0*cyv0-2.0*cii*cos(pv0)*cjv0;
            cyv1 = cfac1*cyv1-2.0*cii*cos(pv1)*cjv1;
            cjv0 /= cfac0;
            cjv1 /= cfac1;
        }
        else if (imag(z) > 0.0) {
            cyv0 = cyv0/cfac0+2.0*cii*cos(pv0)*cjv0;
            cyv1 = cyv1/cfac1+2.0*cii*cos(pv1)*cjv1;
            cjv0 *= cfac0;
            cjv1 *= cfac1;
        }
    }
    cjv[0] = cjv0;
    cjv[1] = cjv1;
    if ((n >= 2) && (n <= (int)(0.25*a0))) {
        cf0 = cjv0;
        cf1 = cjv1;
        for (k=2;k<= n;k++) {
            cf = 2.0*(k+v0-1.0)*cf1/z-cf0;
            cjv[k] = cf;
            cf0 = cf1;
            cf1 = cf;
        }
    }
    else if (n >= 2) {
        m = msta1(a0,200);
        if (m < n) n = m;
        else  m = msta2(a0,n,15);
        cf2 = czero;
        cf1 = complex<double>(1.0e-100,0.0);
        for (k=m;k>=0;k--) {
            cf = 2.0*(v0+k+1.0)*cf1/z-cf2;
            if (k <= n) cjv[k] = cf;
            cf2 = cf1;
            cf1 = cf;
        }
        if (abs(cjv0) > abs(cjv1)) cs = cjv0/cf;
        else cs = cjv1/cf2;
        for (k=0;k<=n;k++) {
            cjv[k] *= cs;
        }
    }
    cjvp[0] = v0*cjv[0]/z-cjv[1];
    for (k=1;k<=n;k++) {
        cjvp[k] = -(k+v0)*cjv[k]/z+cjv[k-1];
    }
    cyv[0] = cyv0;
    cyv[1] = cyv1;
    ya0 = abs(cyv0);
    lb = 0;
    cg0 = cyv0;
    cg1 = cyv1;
    for (k=2;k<=n;k++) {
        cyk = 2.0*(v0+k-1.0)*cg1/z-cg0;
        yak = abs(cyk);
        ya1 = abs(cg0);
        if ((yak < ya0) && (yak< ya1)) lb = k;
        cyv[k] = cyk;
        cg0 = cg1;
        cg1 = cyk;
    }
    lb0 = 0;
    if ((lb > 4) && (imag(z) != 0.0)) {
        while(lb != lb0) {
            ch2 = cone;
            ch1 = czero;
            lb0 = lb;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*(k+v0)*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp12 = ch0;
            cp22 = ch2;
            ch2 = czero;
            ch1 = cone;
            for (k=lb;k>=1;k--) {
                ch0 = 2.0*(k+v0)*ch1/z-ch2;
                ch2 = ch1;
                ch1 = ch0;
            }
            cp11 = ch0;
            cp21 = ch2;
            if (lb == n)
                cjv[lb+1] = 2.0*(lb+v0)*cjv[lb]/z-cjv[lb-1];
            if (abs(cjv[0]) > abs(cjv[1])) {
                cyv[lb+1] = (cjv[lb+1]*cyv0-2.0*cp11/(M_PI*z))/cjv[0];
                cyv[lb] = (cjv[lb]*cyv0+2.0*cp12/(M_PI*z))/cjv[0];
            }
            else {
                cyv[lb+1] = (cjv[lb+1]*cyv1-2.0*cp21/(M_PI*z))/cjv[1];
                cyv[lb] = (cjv[lb]*cyv1+2.0*cp22/(M_PI*z))/cjv[1];
            }
            cyl2 = cyv[lb+1];
            cyl1 = cyv[lb];
            for (k=lb-1;k>=0;k--) {
                cylk = 2.0*(k+v0+1.0)*cyl1/z-cyl2;
                cyv[k] = cylk;
                cyl2 = cyl1;
                cyl1 = cylk;
            }
            cyl1 = cyv[lb];
            cyl2 = cyv[lb+1];
            for (k=lb+1;k<n;k++) {
                cylk = 2.0*(k+v0)*cyl2/z-cyl1;
                cyv[k+1] = cylk;
                cyl1 = cyl2;
                cyl2 = cylk;
            }
            for (k=2;k<=n;k++) {
                wa = abs(cyv[k]);
                if (wa < abs(cyv[k-1])) lb = k;
            }
        }
    }
    cyvp[0] = v0*cyv[0]/z-cyv[1];
    for (k=1;k<=n;k++) {
        cyvp[k] = cyv[k-1]-(k+v0)*cyv[k]/z;
    }
    vm = n+v0;
    return 0;
}


