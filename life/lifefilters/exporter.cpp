#include <life/lifefilters/exporter.hpp>
namespace LifeV
{
ExporterData::ExporterData(const  ExporterData::Type type, const std::string prefix, const vector_ptrtype& vr, UInt start, UInt dim, UInt steady):
    M_prefix(prefix),
    M_vr(vr),
    M_dim(dim),
    M_start( start ),
    M_type(type),
    M_steady(steady)
{};

std::string ExporterData::prefix() const
{
    return M_prefix;
}

Real ExporterData::operator()(const UInt i) const
{
    return (*M_vr)[i];
}

Real& ExporterData::operator()(const UInt i)
{
    return (*M_vr)[i];
}

UInt ExporterData::dim() const {
    return M_dim;
}

ExporterData::Type ExporterData::type() const {
    return M_type;
}

//! returns Scalar or Vector strings
std::string ExporterData::typeName() const
{
    switch (M_type) {
    case Scalar:
        return "Scalar";
    case Vector:
        return "Vector";
    }

    return "ERROR string";
}

//! returns 1 (if Scalar) or 3 (if Vector)
UInt ExporterData::typeDim() const
{
    switch (M_type) {
    case Scalar:
        return 1;
    case Vector:
        return 3;
    }

    return 0;

}


}
