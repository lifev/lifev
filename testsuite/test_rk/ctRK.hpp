#ifndef __CTRK_H
#define __CTRK_H 1


/*!
 * \class CTRK
 * \brief Class for running Chorin-Temam methods with RK2 time stepping.
 *        Uses opaque pointer to hide implementation of specific case study stuff
 *        from the run itself.
 *
 * @author
 * @see
 */

struct CTRKcaseBase;

class CTRK
{
  public:


    /** @name Typedefs
     */
    //@{

    //@}

    /** @name Constructors, destructor
     */
    //@{

    CTRK( int argc,
          char** argv );

    ~CTRK()
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
    boost::shared_ptr<CTRKcaseBase> C_case;

};

#endif /* __CTRK_H */
