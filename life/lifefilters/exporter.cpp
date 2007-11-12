#include <life/lifefilters/exporter.hpp>
namespace LifeV
{
ExporterData::ExporterData(const  ExporterData::Type type, const std::string prefix, const vector_ptrtype& vr, UInt start, UInt dim):
    M_prefix(prefix),
    M_vr(vr),
    M_dim(dim),
    M_start( start ),
    M_type(type)
{}

std::string ExporterData::prefix() const
{
    return M_prefix;
}

Real ExporterData::operator()(const UInt i) const
{
    return (*M_vr)[i];
}

UInt ExporterData::dim() const {
    return M_dim;
}

ExporterData::Type ExporterData::type() const {
    return M_type;
}

}
