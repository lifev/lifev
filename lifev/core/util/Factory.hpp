//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
  @file
  @brief Factory class

  @date 4-10-2004
  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef FACTORY_H
#define FACTORY_H 1


#include <stdexcept>

#include <boost/bind.hpp>


#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/LifeDebug.hpp>
#include <lifev/core/util/FactoryTypeInfo.hpp>

namespace LifeV
{

//! @struct FactoryDefaultError
/*!
  Manages the "Unknown Type" error in an object factory.
*/
template <class AbstractProduct>
struct FactoryDefaultError
{
    class Exception :  public std::exception
    {
    public:
        //! @name Constructor and destructor
        //@{

        Exception ( const std::string& id ) : std::exception(), M_exception()
        {
            std::ostringstream ex_str;
            ex_str << "[factory] Unknown Type : " + id;
            M_exception = ex_str.str();
        }

        ~Exception() throw() {}

        //@}


        //! @name  Methods
        //@{

        const char* what() const throw ()
        {
            return M_exception.c_str();
        }

        //@}

    private:
        std::string M_exception;
    };

    static AbstractProduct* onUnknownType (const std::string& id )
    {
        throw Exception ( id );
    }
};





/*!
  @class factory
  @brief Implements a generic object factory

  @sa factoryDefaultError, factoryClone, FactoryTypeInfo
*/
template < class AbstractProduct, typename IdentifierType,
         typename ProductCreator = boost::function<AbstractProduct*() >,
         template<class> class FactoryErrorPolicy = FactoryDefaultError >
class Factory : public FactoryErrorPolicy<AbstractProduct>
{
public:
    //! @name Public Typedefs
    //@{

    typedef IdentifierType identifier_Type;
    typedef AbstractProduct product_Type;
    typedef ProductCreator creator_Type;
    typedef FactoryErrorPolicy<product_Type> super;

    //@}


    //! @name Constructor and destructor
    //@{

    Factory() {}

    virtual ~Factory() {}

    //@}


    //! @name  Methods
    //@{

    /**
     * Register a product.
     *
     * A product is composed of an identifier (typically a
     * std::string) and a functor that will create the associated
     * object.
     *
     * @param id identifier for the object to be registered
     * @param creator the functor that will create the registered
     * object
     *
     * @return true if registration went fine, false otherwise
     */
    bool registerProduct ( const identifier_Type& id, creator_Type creator )
    {
#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 2200 ) << "Registered type with id : " << id << "\n";
#endif
        return M_associations.insert ( typename productId_Type::value_type ( id, creator ) ).second;
    }

    /**
     * Unregister a product
     *
     * @param id
     * @sa registerProduct
     * @return true if unregistration went fine, false otherwise
     */
    bool unregisterProduct ( const identifier_Type& id )
    {
#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 2200 ) << "Unregistered type with id : " << id << "\n";
#endif
        return M_associations.erase ( id ) == 1;
    }

    /**
     * Create an object from a product registered in the factory using
     * identifier \c id
     *
     * @param id identifier of the product to instantiate
     *
     * @return the object associate with \c id
     */
    product_Type* createObject ( const identifier_Type& id )
    {
        typename productId_Type::const_iterator i = M_associations.find ( id );
        if (i != M_associations.end() )
        {
#ifdef HAVE_LIFEV_DEBUG
            debugStream ( 2200 ) << "Creating type with id : " << id << "\n";
#endif
            return (i->second) ();
        }
#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 2200 ) << "Unknown type with id : " << id << "\n";
#endif
        return super::onUnknownType ( id );
    }

    /**
     * Create an object from a product registered in the factory using
     * identifier \c id (of type enum) and a map to catch the exception.
     *
     * @param id identifier of the product to instantiate
     *
     * @return the object associate with \c id
     */
    template< typename map_Type >
    product_Type* createObject ( const identifier_Type& id, const map_Type& map )
    {
        typename productId_Type::const_iterator i = M_associations.find ( id );
        if ( i != M_associations.end() )
        {
#ifdef HAVE_LIFEV_DEBUG
            debugStream ( 2200 ) << "Creating type with id : " << enum2String ( id, map ) << "\n";
#endif
            return (i->second) ();
        }
#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 2200 ) << "Unknown type with id : " << enum2String ( id, map ) << "\n";
#endif
        return super::onUnknownType ( enum2String ( id, map ) );
    }

    //@}

private:

    //! @name Private typedefs
    //@{

    typedef std::map<identifier_Type, creator_Type> productId_Type;

    //@}

    productId_Type M_associations;
};


}
#endif // FACTORY_H
