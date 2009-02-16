#ifndef __CT_H
#define __CT_H 1

#include <life/lifecore/application.hpp>

/*!
 * \class CT
 * \brief Class for running Chorin-Temam / Projection methods. 
 *        Uses opaque pointer to hide implementation of specific case study stuff
 *        from the run itself. 
 *
 * @author
 * @see
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
    Epetra_Comm *M_comm;
    struct CTcase;
    boost::shared_ptr<CTcase> C_case;

};

#endif /* __CT_H */
