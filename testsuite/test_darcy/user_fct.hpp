Real g1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
Real g2(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
Real g3(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

class SourceFct
{
public:
  inline Real operator()(Real x,Real y,Real z,int ic=0) const {
    return 0.;
  }
};
