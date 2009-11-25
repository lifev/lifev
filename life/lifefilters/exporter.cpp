#include <life/lifefilters/exporter.hpp>
namespace LifeV
{
ExporterData::ExporterData(const  ExporterData::Type type, const std::string variableName, const vector_ptrtype& vr, UInt start, UInt size, UInt steady):
    M_variableName(variableName),
    M_vr(vr),
    M_size(size),
    M_start( start ),
    M_type(type),
    M_steady(steady)
{};

const std::string& ExporterData::variableName() const
{
    return M_variableName;
}

Real ExporterData::operator()(const UInt i) const
{
    return (*M_vr)[i];
}

Real& ExporterData::operator()(const UInt i)
{
    return (*M_vr)[i];
}

const UInt& ExporterData::size() const {
    return M_size;
}

const ExporterData::Type& ExporterData::type() const {
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
