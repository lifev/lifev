/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Daniele A. di Pietro <dipietro@unibg.it>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2004-10-12

 Copyright (C) 2004 EPFL

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file connectivity.hpp
   \author Daniele A. di Pietro <dipietro@unibg.it>
           Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-10-12
*/
#ifndef _CONNECTIVITY_HPP_
#define _CONNECTIVITY_HPP_

#include <set>

namespace LifeV
{

/*! \class Connectivity
 *  Defines the interactions between matrix entries
 */
class Connectivity {
public:

    typedef std::pair<ID, ID> Connection;
    
    class RowIterator {
    public:
        RowIterator() {};
        const Connection& operator*() { return *M_setIterator; }
        const Connection* operator->(){ return &(*M_setIterator); }
        RowIterator& operator++() { M_setIterator++; return *this; }
        bool operator==( const RowIterator& rhs ) const
            {
                return rhs.M_setIterator == M_setIterator;
            }
        bool operator!=( const RowIterator& rhs ) const
            {
                return !(rhs==*this);
            }
        friend class Connectivity;
    private:
        RowIterator(const std::set<Connection>::iterator& setIterator)
            : M_setIterator(setIterator) {}
        std::set<Connection>::iterator M_setIterator;
    };
    
    //! returns a RowIterator pointing to the first Connection
    RowIterator rowBegin() const
        {
            return RowIterator(M_connectionSet.begin());
        }
    
    //! returns a RowIterator pointing to the end of the Connection container
    RowIterator rowEnd() const
        {
            return RowIterator(M_connectionSet.end());
        }

protected:

    /*! Insert a new Connection. This method is protected, only derived classes
     *  should insert new Connections
     */
    void insert(Connection newConnection)
        { 
            M_connectionSet.insert(newConnection); 
        }
    /*! Insert a new Connection defined by the IDs id2 and id2. This method is
     *  protected, only derived classes should insert new Connections
     */
    void insert(ID id1, ID id2)
        {
            M_connectionSet.insert(make_pair(id1, id2));
        }

private:
    //! the set of connections
    std::set<Connection> M_connectionSet;
};

//! defines connectivity of standard continuous fem
template<typename DOF>
class ContinuousFemConnectivity : public Connectivity
{
public:
    ContinuousFemConnectivity(const DOF& dof);
};

class TestConnectivity : public Connectivity
{
public:
    TestConnectivity()
        {
            this->insert(Connection(3,4));
            this->insert(1,2);
            this->insert(1,1);
        }
};

}

#endif /* _CONNECTIVITY_HPP_ */
