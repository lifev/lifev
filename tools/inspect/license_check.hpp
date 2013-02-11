//  license_check header  ----------------------------------------------------//

//  Copyright Beman Dawes 2002, 2003.
//  Copyright Rene Rivera 2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_LICENSE_CHECK_HPP
#define BOOST_LICENSE_CHECK_HPP

#include "inspector.hpp"

namespace boost
{
namespace inspect
{
class license_check : public source_inspector
{
    long m_files_with_errors;
public:

    license_check();
    virtual const char* name() const
    {
        return "license-check";
    }
    virtual const char* desc() const
    {
        return "missing license info (LGPL or GPL)";
    }

    virtual void inspect (
        const std::string& library_name,
        const path& full_path,
        const std::string& contents );

    virtual ~license_check()
    {
        std::cout << "  "
                  << m_files_with_errors << " files missing license info (LGPL or GPL)\n";
    }
};
}
}

#endif // BOOST_LICENSE_CHECK_HPP
