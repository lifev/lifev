namespace LifeV {

Real uexact(const Real& t, const Real& x, const Real& y, const Real& z,
            const ID& i);

Real f(const Real& t, const Real& x, const Real& y, const Real& z,
       const ID& i) {
    Real nu=0.0001;
    Real sigma=0.0;
    Real a=0.75;
    Real b=0.75;
    switch(i) {
        case 1:
            return
                -2.0*nu*exp(a*x-a*z+b*y-b*z)*b*b*b
                +2.0*nu*exp(a*z-a*y+b*x-b*y)*a*b*b
                +2.0*nu*exp(a*z-a*y+b*x-b*y)*a*a*b
                -2.0*nu*exp(a*x-a*z+b*y-b*z)*b*a*a
                +2.0*nu*exp(a*z-a*y+b*x-b*y)*a*a*a
                -2.0*nu*exp(a*x-a*z+b*y-b*z)*b*b*a
                + sigma*uexact(t,x,y,z,i);
            break;
        case 2:
            return
                -2.0*nu*exp(a*y-a*x+b*z-b*x)*b*b*b
                +2.0*nu*exp(a*x-a*z+b*y-b*z)*b*a*a
                +2.0*nu*exp(a*x-a*z+b*y-b*z)*a*a*a
                +2.0*nu*exp(a*x-a*z+b*y-b*z)*b*b*a
                -2.0*nu*exp(a*y-a*x+b*z-b*x)*b*b*a
                -2.0*nu*exp(a*y-a*x+b*z-b*x)*b*a*a
                + sigma*uexact(t,x,y,z,i);
            break;
        case 3:
            return
                -2.0*nu*exp(a*z-a*y+b*x-b*y)*b*b*b
                +2.0*nu*exp(a*y-a*x+b*z-b*x)*a*a*a
                -2.0*nu*exp(a*z-a*y+b*x-b*y)*a*b*b
                -2.0*nu*exp(a*z-a*y+b*x-b*y)*a*a*b
                +2.0*nu*exp(a*y-a*x+b*z-b*x)*b*b*a
                +2.0*nu*exp(a*y-a*x+b*z-b*x)*b*a*a
                + sigma*uexact(t,x,y,z,i);
            break;
    }
    exit(1);
}

Real fZero(const Real& t, const Real& x, const Real& y, const Real& z,
           const ID& i) {
    return 0;
}

Real uexact(const Real& t, const Real& x, const Real& y, const Real& z,
            const ID& i) {
    Real a=0.75;
    Real b=0.75;

    switch(i) {
        case 1:
            return
                b*exp(a*(x-z)+b*(y-z))-
                a*exp(a*(z-y)+b*(x-y));
            break;
        case 2:
            return
                b*exp(a*(y-x)+b*(z-x))-
                a*exp(a*(x-z)+b*(y-z));
            break;
        case 3:
            return
                b*exp(a*(z-y)+b*(x-y))-
                a*exp(a*(y-x)+b*(z-x));
            break;
    }

    exit(1);
}


Real pexact(const Real& t, const Real& x, const Real& y, const Real& z,
            const ID& i) {
    Real a=0.75;
    Real b=0.75;
    return (a*a+b*b+a*b)*(exp(a*(x-y)+b*(x-z))+
                          exp(a*(y-z)+b*(y-x))+
                          exp(a*(z-x)+b*(z-y)));

}

Real xexact(const Real& t, const Real& x, const Real& y, const Real& z,
            const ID& i) {

    switch(i) {
        case 1:
        case 2:
        case 3:
            return uexact(t, x, y, z, i);
            break;
        case 4:
            return pexact(t, x, y, z, 1);
            break;
        default:
            exit(1);
    }
}


// Initial velocity
Real u0(const Real& t, const Real& x, const Real& y, const Real& z,
        const ID& i) {
    // return 0.0;
    return uexact(t,x,y,z,i);
}

} // namespace LifeV
