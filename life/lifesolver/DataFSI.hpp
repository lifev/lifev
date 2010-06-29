//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2009-2010 EPFL, Politecnico di Milano

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief DataFSI - File containing a data container for FSI problems
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 10-06-2010
 *
 *  version 1 + rand(0,1), 3 months earlier
 *  @author Gilles fourestey <gilles.fourestey@epfl.ch>
 */

#ifndef DATAFSI_H
#define DATAFSI_H

#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifesolver/dataElasticStructure.hpp>

namespace LifeV {

enum Preconditioner
{
    NO_PRECONDITIONER = -1, //!< No preconditioner
    NEUMANN_DIRICHLET,      //!< Neuman Dirichlet
    DIRICHLET_NEUMANN,      //!< Dirichlet Neumann
    NEUMANN_NEUMANN,        //!< Neumann Neumann
    NEWTON                  //!< Newton
};

enum DDNPreconditioner
{
    DDN_NO_PRECONDITIONER = -1, //!< No preconditioner
    DDN_NEUMANN_DIRICHLET,      //!< Neuman Dirichlet
    DDN_DIRICHLET_NEUMANN,      //!< Dirichlet Neumann
    DDN_NEUMANN_NEUMANN         //!< Neumann Neumann
};

//! DataFSI - Data container for FSI problems
/*!
 *  @author Cristiano Malossi
 */

class DataFSI
{
public:

    //! @name Type definitions
    //@{

    typedef DataNavierStokes                  dataFluid_Type;
    typedef boost::shared_ptr< dataFluid_Type >     dataFluid_PtrType;

    typedef DataElasticStructure              dataSolid_Type;
    typedef boost::shared_ptr< dataSolid_Type >     dataSolid_PtrType;

    //@}


	//! @name Constructors & Destructor
	//@{

    //! Empty Constructor
    DataFSI();

	//! Copy constructor
	/*!
	 * @param DataFSI - DataFSI
	 */
	DataFSI( const DataFSI& DataFSI );

	//@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param DataFSI - DataFSI
     */
    DataFSI& operator=( const DataFSI& DataFSI );

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile - data file
     */
    void setup( const GetPot& dataFile, const std::string& section = "problem" );

    //! Display the values
    void showMe( std::ostream& output = std::cout );

    bool isMonolithic();

    //@}


    //! @name Set methods
    //@{

    //! Set data fluid container
    /*!
     * @param dataFluid shared_ptr to dataFluid container
     */
    inline void setDataFluid( const dataFluid_PtrType dataFluid ) { M_dataFluid = dataFluid; }

    //! Set data solid container
    /*!
     * @param dataFluid shared_ptr to dataSolid container
     */
    inline void setDataSolid( const dataSolid_PtrType dataSolid ) { M_dataSolid = dataSolid; }

    //! Set preconditioner type
    /*!
     * @param preconditioner preconditioner type
     */
    inline void setPreconditioner( const Preconditioner& preconditioner ) {M_preconditioner = preconditioner; }

    //@}


    //! @name Get methods
    //@{

    //! Get data fluid container
    /*!
     * @return shared_ptr to dataFluid container
     */
    inline dataFluid_PtrType dataFluid() const { return M_dataFluid; }

    //! Get data solid container
    /*!
     * @return shared_ptr to dataSolid container
     */
    inline dataSolid_PtrType dataSolid() const { return M_dataSolid; }

    //! Get maximum number of subiterations
    /*!
     * @return maximum number of subiterations
     */
    inline UInt maxSubIterationNumber() const { return M_maxSubIterationNumber; }

    //! Get absolute tolerance
    /*!
     * @return absolute tolerance
     */
    inline Real absoluteTolerance() const { return M_absoluteTolerance; }

    //! Get relative tolerance
    /*!
     * @return relative tolerance
     */
    inline Real relativeTolerance() const { return M_relativeTolerance; }

    //! Get error tolerance
    /*!
     * @return error tolerance
     */
    inline Real errorTolerance() const { return M_errorTolerance; }

    //! Get linesearch
    /*!
     * @return linesearch
     */
    inline Int linesearch() const { return M_linesearch; }

    //! Get preconditioner type
    /*!
     * @return preconditioner type
     */
    inline Preconditioner preconditioner() const { return M_preconditioner; }

    //! Get DDNpreconditioner type
    /*!
     * @return DDNpreconditioner type
     */
    inline DDNPreconditioner DDNpreconditioner() const { return M_DDNpreconditioner; }

    //! Get DDBlockPreconditioner type
    /*!
     * @return DDBlockPreconditioner type
     */
    inline UInt DDBlockPreconditioner() const { return M_DDBlockPreconditioner; }

    //! Get method type
    /*!
     * @return method type
     */
    inline std::string method() const { return M_method; }

    //! Get algorithm type
    /*!
     * @return algorithm type
     */
    inline std::string algorithm() const { return M_algorithm; }

    //! Get default omega for Aitken iterations
    /*!
     * @return default omega for Aitken iterations
     */
    inline Real defaultOmega() const { return M_defaultOmega; }

    //! Get update every
    /*!
     * If M_updateEvery == 1, normal fixedPoint algorithm
     * If M_updateEvery  > 1, recompute computational domain every M_updateEvery iterations (transpiration)
     * If M_updateEvery <= 0, recompute computational domain and matrices only at first subiteration (semi-implicit)
     *
     * @return updateEvery value
     */
    inline int updateEvery() const { return M_updateEvery; }

    //! Get algorithm type
    /*!
     * @return algorithm type
     */
    inline int fluidInterfaceFlag() const { return M_fluidInterfaceFlag; }

    //! Get algorithm type
    /*!
     * @return algorithm type
     */
    inline int solidInterfaceFlag() const { return M_solidInterfaceFlag; }

    //! Get algorithm type
    /*!
     * @return algorithm type
     */
    inline int structureInterfaceFlag() const { return M_structureInterfaceFlag; }

    //! Get algorithm type
    /*!
     * @return algorithm type
     */
    inline int harmonicInterfaceFlag() const { return M_harmonicInterfaceFlag; }

    //! Get algorithm type
    /*!
     * @return algorithm type
     */
    inline Real interfaceTolerance() const { return M_interfaceTolerance; }

    //! Get Robin-Neumann coupling flag
    /*!
     * @return Robin-Neumann coupling flag
     */
    inline bool RobinNeumannCoupling() const { return M_RobinNeumannCoupling; }

    //! Get Robin-Neumann fluid coefficient
    /*!
     * @return Robin-Neumann fluid coefficient
     */
    const Real& RobinNeumannFluidCoefficient() const { return M_RobinNeumannFluidCoefficient; }

    //! Get Robin-Neumann solid coefficient
    /*!
     * @return Robin-Neumann solid coefficient
     */
    const Real& RobinNeumannSolidCoefficient() const { return M_RobinNeumannSolidCoefficient; }

    //@}

private:

    dataFluid_PtrType             M_dataFluid;
    dataSolid_PtrType             M_dataSolid;

    // Problem - Non Linear Richardson parameters
    UInt                          M_maxSubIterationNumber;
    Real                          M_absoluteTolerance;
    Real                          M_relativeTolerance;
    Real                          M_errorTolerance;
    Int                           M_linesearch;

    // Problem - Preconditioner
    Preconditioner                M_preconditioner;
    DDNPreconditioner             M_DDNpreconditioner;
    UInt                          M_DDBlockPreconditioner;

    // Problem - Methods
    std::string                   M_method;
    std::string                   M_algorithm;

    // Problem - FixPoint / EJ
    Real                          M_defaultOmega;
    int                           M_updateEvery;

    // Interface
    int                           M_fluidInterfaceFlag;
    int                           M_solidInterfaceFlag;
    int                           M_structureInterfaceFlag;
    int                           M_harmonicInterfaceFlag;
    Real                          M_interfaceTolerance;

    bool                          M_RobinNeumannCoupling;
    Real                          M_RobinNeumannFluidCoefficient;
    Real                          M_RobinNeumannSolidCoefficient;
};



} // end namespace LifeV

#endif // end DATAFSI_H
