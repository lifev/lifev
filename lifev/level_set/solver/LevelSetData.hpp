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
    @brief Contains the data container for the LevelSetSolver

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 10-11-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */

#ifndef DATALEVELSET_H
#define DATALEVELSET_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/fem/TimeAdvanceData.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

//! dataLevelSet - Container for the data for the level set solver
/*!

  There are two ways of filling this data container: either by reading a file
  and using the setup method, or by setting the data one by one using the
  setters defined in the class. Mixture of the two solutions are also possible.

  <b> Reading from a data file </b>

  If one read the data from a file, the following informations
  can be provided:

  \verbatim
  [section]

      [./time_discretization]

          #everything for the TimeData class

      [../]

      stabilization = ip  #options: none, ip

      [./ip]

          coefficient = 0.5
          treatment = implicit  #options: implicit, semi-implicit, explicit

      [../]

  [../]
  \endverbatim

  The code that one need then (usually in the main file) is simply

  \code

  GetPot dataFile( ... );

  DataLevelSet data_level_set;
  data_level_set.setup(dataFile,"section");

  \endcode

  <b> Using the setters </b>

  This way allows you to define all the members
  of the class. Simply use the setters defined.

  For the stabilization and its treatment, prefer
  using the setters with the enumerated types, this makes
  code more robust with respect to typos.

  <b> Mixing the two solutions </b>

  Mixing the two solutions might be a good idea is some cases,
  for example to write the time discretization only once in a
  data file for a problem that consists in coupling smaller
  problems. Here, the idea would be to initialize this data
  container with the dataFile (this will put default parameters
  for the time discretization), create a TimeData (outside
  DataLevelSet) and set this new TimeData in DataLevelSet. As
  DataLevelSet will retain a pointer to the TimeData, any change
  to it will be repercuted in the DataLevelSet. Then, only one
  TimeData needs to be managed.

  Beware that, when mixing the two solutions, any call to the
  setup will overwrite ALL the previous informations stored
  in the DataLevelSet.

    @author Samuel Quinodoz
 */
class DataLevelSet
{
public:

    //! @name Public Types
    //@{

    typedef TimeData                               time_Type;
    typedef std::shared_ptr<time_Type>           timePtr_Type;

    typedef TimeAdvanceData                        timeAdvance_Type;
    typedef std::shared_ptr<timeAdvance_Type>    timeAdvancePtr_Type;

    //! \enum Enumerated type for the stabilization
    enum stabilization_type
    {
        NONE, //!< No stabilization will be added
        IP //!< The IP stabiilzation will be used
    };

    //! \enum Enumerated type for the treatment of the stabilization (if IP is chosen)
    enum IPTreatment_type
    {
        IMPLICIT, //!< Fully implicit stabilization (standard procedure)
        SEMI_IMPLICIT, //!< Semi-implicit, best choice for faster computations
        EXPLICIT //!< Completly explicit, not very usefull but available.
    };

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    DataLevelSet();

    //! Destructor
    virtual ~DataLevelSet() {};

    //@}


    //! @name Methods
    //@{

    //! Fill the DataLevelSet with the informations of the GetPot object
    /*!
      Using this method overrides all the previously stored informations!
     */
    void setup ( const GetPot& dataFile, const std::string& section = "level-set");

    //! ShowMe method
    void showMe (std::ostream& out = std::cout) const ;

    //@}


    //! @name Set Methods
    //@{

    //! Set data time container
    /*!
     * @param TimeData shared_ptr to TimeData container
     */
    void setTimeData ( const timePtr_Type timeData )
    {
        M_time = timeData;
    }

    //! Set data time advance container
    /*!
     * @param timeAdvanceData shared_ptr to TimeAdvanceData container
     */
    void setTimeAdvanceData ( const timeAdvancePtr_Type timeAdvanceData )
    {
        M_timeAdvance = timeAdvanceData;
    }

    //! Set the stabilization type
    /*!
      This set the stabilization to the one given in the string format.
      Accepted possibilities are "none" and "ip".
      @Warning: prefer using the other version of setStabilization
      that is not sensible to typos.
    */
    void setStabilization (const std::string& stab);

    //! Set the stabilization type
    inline void setStabilization (const stabilization_type& stab)
    {
        M_stabilization = stab;
    }

    //! Set the treatment for the IP stabilization
    /*!
      This set the treatment of the ip stabilization to the one given in the string format.
      Accepted possibilities are "implicit", "semi-implicit" and "explicit".
      @Warning: prefer using the other version of setIPTreatment
      that is not sensible to typos.
    */
    void setIPTreatment (const std::string& treat);

    //! Set the IP treatment
    void setIPTreatment (const IPTreatment_type& treat)
    {
        M_IPTreatment = treat;
    }

    //! Set the IP coefficient
    inline void setIPCoef (const Real& coef)
    {
        M_IPCoef = coef;
    };

    //@}


    //! @name Get Methods
    //@{

    //! Get data time container
    /*!
     * @return shared_ptr to TimeData container
     */
    timePtr_Type dataTime() const
    {
        return M_time;
    }

    //! Get data time advance container
    /*!
     * @return shared_ptr to TimeAdvanceData container
     */
    timeAdvancePtr_Type dataTimeAdvance() const
    {
        return M_timeAdvance;
    }

    //! Getter for the stabilization type
    inline stabilization_type stabilization() const
    {
        return M_stabilization;
    };

    //! Getter for the IP treatment
    inline IPTreatment_type IPTreatment() const
    {
        return M_IPTreatment;
    };

    //! Getter for the IP coefficient
    inline Real IPCoef() const
    {
        return M_IPCoef;
    };

    //@}

private:

    // No copy
    DataLevelSet (const DataLevelSet&);

    // Data for the time
    timePtr_Type        M_time;
    timeAdvancePtr_Type M_timeAdvance;

    // Stabilization type
    stabilization_type M_stabilization;

    // IP Stabilization treatment
    IPTreatment_type M_IPTreatment;

    // Coefficient for the IP
    Real M_IPCoef;

};


} // Namespace LifeV

#endif /* DATALEVELSET_H */
