#ifndef __CT_H
#define __CT_H 1

/*!
 * \class CT
 * \brief Class for running Chorin-Temam / Projection methods.
 *        Uses opaque pointer to hide implementation of specific case study stuff
 *        from the run itself.
 *
 * @author
 * @see
 */

struct CTcaseBase;

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
        char** argv );

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
    boost::shared_ptr<CTcaseBase> C_case;

};

#endif /* __CT_H */
