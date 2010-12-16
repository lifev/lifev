//  cbessik.cpp -- complex modified Bessel functions.
//  Algorithms and coefficient values from "Computation of Special
//  Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
#include <complex>
#include <life/lifefunctions/bessel/bessel.hpp>

using namespace std;

static complex<double> cii(0.0,1.0);
static complex<double> czero(0.0,0.0);
static complex<double> cone(1.0,0.0);

double gamma(double x);

int cbessik01(complex<double>z,complex<double>&ci0,complex<double>&ci1,
    complex<double>&ck0,complex<double>&ck1,complex<double>&ci0p,
    complex<double>&ci1p,complex<double>&ck0p,complex<double>&ck1p)
{
    complex<double> z1,z2,zr,zr2,cr,ca,cb,cs,ct,cw;
    double a0,w0;
    int k,kz;
    static double a[] = {
        0.125,
        7.03125e-2,
        7.32421875e-2,
        1.1215209960938e-1,
        2.2710800170898e-1,
        5.7250142097473e-1,
        1.7277275025845,
        6.0740420012735,
        2.4380529699556e1,
        1.1001714026925e2,
        5.5133589612202e2,
        3.0380905109224e3};
    static double b[] = {
        -0.375,
        -1.171875e-1,
        -1.025390625e-1,
        -1.4419555664063e-1,
        -2.7757644653320e-1,
        -6.7659258842468e-1,
        -1.9935317337513,
        -6.8839142681099,
        -2.7248827311269e1,
        -1.2159789187654e2,
        -6.0384407670507e2,
        -3.3022722944809e3};
    static double a1[] = {
        0.125,
        0.2109375,
        1.0986328125,
        1.1775970458984e1,
        2.1461706161499e2,
        5.9511522710323e3,
        2.3347645606175e5,
        1.2312234987631e7,
        8.401390346421e08,
        7.2031420482627e10};

    a0 = abs(z);
    z2 = z*z;
    z1 = z;
    if (a0 == 0.0) {
        ci0 = cone;
        ci1 = czero;
        ck0 = complex<double> (1e308,0);
        ck1 = complex<double> (1e308,0);
        ci0p = czero;
        ci1p = complex<double>(0.5,0.0);
        ck0p = complex<double>(-1e308,0);
        ck1p = complex<double>(-1e308,0);
        return 0;
    }
    if (real(z) < 0.0) z1 = -z;
    if (a0 <= 18.0) {
        ci0 = cone;
        cr = cone;
        for (k=1;k<=50;k++) {
            cr *= 0.25*z2/(double)(k*k);
            ci0 += cr;
            if (abs(cr/ci0) < eps) break;
        }
        ci1 = cone;
        cr = cone;
        for (k=1;k<=50;k++) {
            cr *= 0.25*z2/(double)(k*(k+1.0));
            ci1 += cr;
            if (abs(cr/ci1) < eps) break;
        }
        ci1 *= 0.5*z1;
    }
    else {
        if (a0 >= 50.0) kz = 7;
        else if (a0 >= 35.0) kz = 9;
        else kz = 12;
        ca = exp(z1)/sqrt(2.0*M_PI*z1);
        ci0 = cone;
        zr = 1.0/z1;
        for (k=0;k<kz;k++) {
            ci0 += a[k]*pow(zr,k+1.0);
        }
        ci0 *= ca;
        ci1 = cone;
        for (k=0;k<kz;k++) {
            ci1 += b[k]*pow(zr,k+1.0);
        }
        ci1 *= ca;
    }
    if (a0 <= 9.0) {
        cs = czero;
        ct = -log(0.5*z1)-el;
        w0 = 0.0;
        cr = cone;
        for (k=1;k<=50;k++) {
            w0 += 1.0/k;
            cr *= 0.25*z2/(double)(k*k);
            cs += cr*(w0+ct);
            if (abs((cs-cw)/cs) < eps) break;
            cw = cs;
        }
        ck0 = ct+cs;
    }
    else {
        cb = 0.5/z1;
        zr2 = 1.0/z2;
        ck0 = cone;
        for (k=0;k<10;k++) {
            ck0 += a1[k]*pow(zr2,k+1.0);
        }
        ck0 *= cb/ci0;
    }
    ck1 = (1.0/z1 - ci1*ck0)/ci0;
    if (real(z) < 0.0) {
        if (imag(z) < 0.0) {
            ck0 += cii*M_PI*ci0;
            ck1 = -ck1+cii*M_PI*ci1;
        }
        else if (imag(z) > 0.0) {
            ck0 -= cii*M_PI*ci0;
            ck1 = -ck1-cii*M_PI*ci1;
        }
        ci1 = -ci1;
    }
    ci0p = ci1;
    ci1p = ci0-1.0*ci1/z;
    ck0p = -ck1;
    ck1p = -ck0-1.0*ck1/z;
    return 0;
}
int cbessikna(int n,complex<double> z,int &nm,complex<double> *ci,
    complex<double> *ck,complex<double> *cip,complex<double> *ckp)
{
    complex<double> ci0,ci1,ck0,ck1,ckk,cf,cf1,cf2,cs;
    double a0;
    int k,m,ecode;
    a0 = abs(z);
    nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            ci[k] = czero;
            ck[k] = complex<double>(-1e308,0);
            cip[k] = czero;
            ckp[k] = complex<double>(1e308,0);
        }
        ci[0] = cone;
        cip[1] = complex<double>(0.5,0.0);
        return 0;
    }
    ecode = cbessik01(z,ci[0],ci[1],ck[0],ck[1],cip[0],cip[1],ckp[0],ckp[1]);
    if (n < 2) return 0;
    ci0 = ci[0];
    ci1 = ci[1];
    ck0 = ck[0];
    ck1 = ck[1];
    m = msta1(a0,200);
    if (m < n) nm = m;
    else m = msta2(a0,n,15);
    cf2 = czero;
    cf1 = complex<double>(1.0e-100,0.0);
    for (k=m;k>=0;k--) {
        cf = 2.0*(k+1.0)*cf1/z+cf2;
        if (k <= nm) ci[k] = cf;
        cf2 = cf1;
        cf1 = cf;
    }
    cs = ci0/cf;
    for (k=0;k<=nm;k++) {
        ci[k] *= cs;
    }
    for (k=2;k<=nm;k++) {
        if (abs(ci[k-1]) > abs(ci[k-2])) {
            ckk = (1.0/z-ci[k]*ck[k-1])/ci[k-1];
        }
        else {
            ckk = (ci[k]*ck[k-2]+2.0*(k-1.0)/(z*z))/ci[k-2];
        }
        ck[k] = ckk;
    }
    for (k=2;k<=nm;k++) {
        cip[k] = ci[k-1]-(double)k*ci[k]/z;
        ckp[k] = -ck[k-1]-(double)k*ck[k]/z;
    }
    return 0;
}
int cbessiknb(int n,complex<double> z,int &nm,complex<double> *ci,
    complex<double> *ck,complex<double> *cip,complex<double> *ckp)
{
    complex<double> z1,cbs,csk0,cf,cf0,cf1,ca0,cbkl;
    complex<double> cg,cg0,cg1,cs0,cs,cr;
    double a0,vt,fac;
    int k,kz,l,m;

    a0 = abs(z);
    nm = n;
    if (a0 < 1.0e-100) {
        for (k=0;k<=n;k++) {
            ci[k] = czero;
            ck[k] = complex<double>(1e308,0);
            cip[k] = czero;
            ckp[k] = complex<double>(-1e308,0);
        }
        ci[0] = complex<double>(1.0,0.0);
        cip[1] = complex<double>(0.5,0.0);
        return 0;
    }
    z1 = z;
    if (real(z) < 0.0) z1 = -z;
    if (n == 0) nm = 1;
    m = msta1(a0,200);
    if (m < nm) nm = m;
    else m = msta2(a0,nm,15);
    cbs = czero;
    csk0 = czero;
    cf0 = czero;
    cf1 = complex<double>(1.0e-100,0.0);
    for (k=m;k>=0;k--) {
        cf = 2.0*(k+1.0)*cf1/z1+cf0;
        if (k <=nm) ci[k] = cf;
        if ((k != 0) && (k == 2*(k>>1)))  csk0 += 4.0*cf/(double)k;
        cbs += 2.0*cf;
        cf0 = cf1;
        cf1 = cf;
    }
    cs0 = exp(z1)/(cbs-cf);
    for (k=0;k<=nm;k++) {
        ci[k] *= cs0;
    }
    if (a0 <= 9.0) {
        ck[0] = -(log(0.5*z1)+el)*ci[0]+cs0*csk0;
        ck[1] = (1.0/z1-ci[1]*ck[0])/ci[0];
    }
    else {
        ca0 = sqrt(M_PI_2/z1)*exp(-z1);
        if (a0 >= 200.0) kz = 6;
        else if (a0 >= 80.0) kz = 8;
        else if (a0 >= 25.0) kz = 10;
        else kz = 16;
        for (l=0;l<2;l++) {
            cbkl = cone;
            vt = 4.0*l;
            cr = cone;
            for (k=1;k<=kz;k++) {
                cr *= 0.125*(vt-pow(2.0*k-1.0,2.0))/((double)k*z1);
                cbkl += cr;
            }
            ck[l] = ca0*cbkl;
        }
    }
    cg0 = ck[0];
    cg1 = ck[1];
    for (k=2;k<=nm;k++) {
        cg = 2.0*(k-1.0)*cg1/z1+cg0;
        ck[k] = cg;
        cg0 = cg1;
        cg1 = cg;
    }
    if (real(z) < 0.0) {
        fac = 1.0;
        for (k=0;k<=nm;k++) {
            if (imag(z) < 0.0) {
                ck[k] = fac*ck[k]+cii*M_PI*ci[k];
            }
            else {
                ck[k] = fac*ck[k]-cii*M_PI*ci[k];
            }
            ci[k] *= fac;
            fac = -fac;
        }
    }
    cip[0] = ci[1];
    ckp[0] = -ck[1];
    for (k=1;k<=nm;k++) {
        cip[k] = ci[k-1]-(double)k*ci[k]/z;
        ckp[k] = -ck[k-1]-(double)k*ck[k]/z;
    }
    return 0;
}
int cbessikv(double v,complex<double>z,double &vm,complex<double> *civ,
    complex<double> *ckv,complex<double> *civp,complex<double> *ckvp)
{
    complex<double> z1,z2,ca1,ca,cs,cr,ci0,cbi0,cf,cf1,cf2;
    complex<double> ct,cp,cbk0,ca2,cr1,cr2,csu,cws,cb;
    complex<double> cg0,cg1,cgk,cbk1,cvk;
    double a0,v0,v0p,v0n,vt,w0,piv,gap(0),gan;
    int m,n,k,kz;

    a0 = abs(z);
    z1 = z;
    z2 = z*z;
    n = (int)v;
    v0 = v-n;
    piv = M_PI*v0;
    vt = 4.0*v0*v0;
    if (n == 0) n = 1;
    if (a0 < 1e-100) {
        for (k=0;k<=n;k++) {
            civ[k] = czero;
            ckv[k] = complex<double>(-1e308,0);
            civp[k] = czero;
            ckvp[k] = complex<double>(1e308,0);
        }
        if (v0 == 0.0) {
            civ[0] = cone;
            civp[1] = complex<double> (0.5,0.0);
        }
        vm = v;
        return 0;
    }
    if (a0 >= 50.0) kz = 8;
    else if (a0 >= 35.0) kz = 10;
    else kz = 14;
    if (real(z) <= 0.0) z1 = -z;
    if (a0 < 18.0) {
        if (v0 == 0.0) {
            ca1 = cone;
        }
        else {
            v0p = 1.0+v0;
            gap = gamma(v0p);
            ca1 = pow(0.5*z1,v0)/gap;
        }
        ci0 = cone;
        cr = cone;
        for (k=1;k<=50;k++) {
            cr *= 0.25*z2/(k*(k+v0));
            ci0 += cr;
            if (abs(cr/ci0) < eps) break;
        }
        cbi0 = ci0*ca1;
    }
    else {
        ca = exp(z1)/sqrt(2.0*M_PI*z1);
        cs = cone;
        cr = cone;
        for (k=1;k<=kz;k++) {
            cr *= -0.125*(vt-pow(2.0*k-1.0,2.0))/((double)k*z1);
            cs += cr;
        }
        cbi0 = ca*cs;
    }
    m = msta1(a0,200);
    if (m < n) n = m;
    else m = msta2(a0,n,15);
    cf2 = czero;
    cf1 = complex<double>(1.0e-100,0.0);
    for (k=m;k>=0;k--) {
        cf = 2.0*(v0+k+1.0)*cf1/z1+cf2;
        if (k <= n) civ[k] = cf;
        cf2 = cf1;
        cf1 = cf;
    }
    cs = cbi0/cf;
    for (k=0;k<=n;k++) {
        civ[k] *= cs;
    }
    if (a0 <= 9.0) {
        if (v0 == 0.0) {
            ct = -log(0.5*z1)-el;
            cs = czero;
            w0 = 0.0;
            cr = cone;
            for (k=1;k<=50;k++) {
                w0 += 1.0/k;
                cr *= 0.25*z2/(double)(k*k);
                cp = cr*(w0+ct);
                cs += cp;
                if ((k >= 10) && (abs(cp/cs) < eps)) break;
            }
            cbk0 = ct+cs;
        }
        else {
            v0n = 1.0-v0;
            gan = gamma(v0n);
            ca2 = 1.0/(gan*pow(0.5*z1,v0));
            ca1 = pow(0.5*z1,v0)/gap;
            csu = ca2-ca1;
            cr1 = cone;
            cr2 = cone;
            cws = czero;
            for (k=1;k<=50;k++) {
                cr1 *= 0.25*z2/(k*(k-v0));
                cr2 *= 0.25*z2/(k*(k+v0));
                csu += ca2*cr1-ca1*cr2;
                if ((k >= 10) && (abs((cws-csu)/csu) < eps)) break;
                cws = csu;
            }
            cbk0 = csu*M_PI_2/sin(piv);
        }
    }
    else {
        cb = exp(-z1)*sqrt(M_PI_2/z1);
        cs = cone;
        cr = cone;
        for (k=1;k<=kz;k++) {
            cr *= 0.125*(vt-pow(2.0*k-1.0,2.0))/((double)k*z1);
            cs += cr;
        }
        cbk0 = cb*cs;
    }
    cbk1 = (1.0/z1-civ[1]*cbk0)/civ[0];
    ckv[0] = cbk0;
    ckv[1] = cbk1;
    cg0 = cbk0;
    cg1 = cbk1;
    for (k=2;k<=n;k++) {
        cgk = 2.0*(v0+k-1.0)*cg1/z1+cg0;
        ckv[k] = cgk;
        cg0 = cg1;
        cg1 = cgk;
    }
    if (real(z) < 0.0) {
        for (k=0;k<=n;k++) {
            cvk = exp((k+v0)*M_PI*cii);
            if (imag(z) < 0.0) {
                ckv[k] = cvk*ckv[k]+M_PI*cii*civ[k];
                civ[k] /= cvk;
            }
            else if (imag(z) > 0.0) {
                ckv[k] = ckv[k]/cvk-M_PI*cii*civ[k];
                civ[k] *= cvk;
            }
        }
    }
    civp[0] = v0*civ[0]/z+civ[1];
    ckvp[0] = v0*ckv[0]/z-ckv[1];
    for (k=1;k<=n;k++) {
        civp[k] = -(k+v0)*civ[k]/z+civ[k-1];
        ckvp[k] = -(k+v0)*ckv[k]/z-ckv[k-1];
    }
    vm = n+v0;
    return 0;
}
