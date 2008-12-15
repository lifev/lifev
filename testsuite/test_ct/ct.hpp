#ifndef __CT_H
#define __CT_H 1

#include <life/lifecore/application.hpp>

enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

/*!
 * \class CT
 * \brief Test for Chorin Temam Method
 *
 *  @author MG
 *  @see
 */
class CT
{
  public:


    /** @name Typedefs
     */
    //@{

    //@}

    /** @name Constructors, destructor
     */
    //@{

    CT( int argc,
              char** argv,
              LifeV::AboutData const& ad,
              LifeV::po::options_description const& od );

    ~CT()
        {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    void run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> d;
};

#endif /* __CT_H */
