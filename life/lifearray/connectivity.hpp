#ifndef _CONNECTIVITY_HPP_
#define _CONNECTIVITY_HPP_

#include <set>

namespace LifeV
{

/*! \class Connectivity
 *  Defines the interactions between matrix entries
 *  \author Daniele A. di Pietro, Christoph Winkelmann
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
