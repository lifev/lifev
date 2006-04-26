#include <ensight.hpp>
namespace LifeV
{
  EnsightData::EnsightData(const  EnsightData::Type type, const std::string prefix, const VectorRange& vr, const UInt dim):
    M_prefix(prefix),
    M_vr(vr),
    M_dim(dim),
    M_type(type)
  {}

  std::string EnsightData::prefix() const 
  {
    return M_prefix;
  }
  
  Real EnsightData::operator()(const UInt i) const {
    return M_vr[i];
  }
  
  UInt EnsightData::dim() const {
    return M_dim;
  }

  EnsightData::Type EnsightData::type() const {
    return M_type;
  }

}
