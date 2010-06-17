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
template <typename Mesh>
class DataFSI
{
public:

    //! @name Type definitions
    //@{

    typedef DataNavierStokes<Mesh>                  dataFluid_Type;
    typedef boost::shared_ptr< dataFluid_Type >     dataFluid_PtrType;

    typedef DataElasticStructure<Mesh>              dataSolid_Type;
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



// ===================================================
// Constructors
// ===================================================
template <typename Mesh>
DataFSI<Mesh>::DataFSI( ) :
    M_dataFluid                     ( new dataFluid_Type() ),
    M_dataSolid                     ( new dataSolid_Type() ),
    M_maxSubIterationNumber         (),
    M_absoluteTolerance             (),
    M_relativeTolerance             (),
    M_errorTolerance                (),
    M_linesearch                    (),
    M_preconditioner                (),
    M_DDNpreconditioner             (),
    M_DDBlockPreconditioner         (),
    M_method                        (),
    M_algorithm                     (),
    M_defaultOmega                  (),
    M_updateEvery                   (),
    M_fluidInterfaceFlag            (),
    M_solidInterfaceFlag            (),
    M_structureInterfaceFlag        (),
    M_harmonicInterfaceFlag         (),
    M_interfaceTolerance            (),
    M_RobinNeumannCoupling          (),
    M_RobinNeumannFluidCoefficient  (),
    M_RobinNeumannSolidCoefficient  ()
{
}

template <typename Mesh>
DataFSI<Mesh>::DataFSI( const DataFSI& DataFSI ) :
    M_dataFluid                     ( DataFSI.M_dataFluid ),
    M_dataSolid                     ( DataFSI.M_dataSolid ),
    M_maxSubIterationNumber         ( DataFSI.M_maxSubIterationNumber ),
    M_absoluteTolerance             ( DataFSI.M_absoluteTolerance ),
    M_relativeTolerance             ( DataFSI.M_relativeTolerance ),
    M_errorTolerance                ( DataFSI.M_errorTolerance ),
    M_linesearch                    ( DataFSI.M_linesearch ),
    M_preconditioner                ( DataFSI.M_preconditioner ),
    M_DDNpreconditioner             ( DataFSI.M_DDNpreconditioner ),
    M_DDBlockPreconditioner         ( DataFSI.M_DDBlockPreconditioner ),
    M_method                        ( DataFSI.M_method ),
    M_algorithm                     ( DataFSI.M_algorithm ),
    M_defaultOmega                  ( DataFSI.M_defaultOmega ),
    M_updateEvery                   ( DataFSI.M_updateEvery ),
    M_fluidInterfaceFlag            ( DataFSI.M_fluidInterfaceFlag ),
    M_solidInterfaceFlag            ( DataFSI.M_solidInterfaceFlag ),
    M_structureInterfaceFlag        ( DataFSI.M_structureInterfaceFlag ),
    M_harmonicInterfaceFlag         ( DataFSI.M_harmonicInterfaceFlag ),
    M_interfaceTolerance            ( DataFSI.M_interfaceTolerance ),
    M_RobinNeumannCoupling          ( DataFSI.M_RobinNeumannCoupling ),
    M_RobinNeumannFluidCoefficient  ( DataFSI.M_RobinNeumannFluidCoefficient ),
    M_RobinNeumannSolidCoefficient  ( DataFSI.M_RobinNeumannSolidCoefficient )
{
}

// ===================================================
// Methods
// ===================================================
template <typename Mesh>
DataFSI<Mesh>&
DataFSI<Mesh>::operator=( const DataFSI& DataFSI )
{
    if ( this != &DataFSI )
    {
        M_dataFluid                     = DataFSI.M_dataFluid;
        M_dataSolid                     = DataFSI.M_dataSolid;
        M_maxSubIterationNumber         = DataFSI.M_maxSubIterationNumber;
        M_absoluteTolerance             = DataFSI.M_absoluteTolerance;
        M_relativeTolerance             = DataFSI.M_relativeTolerance;
        M_errorTolerance                = DataFSI.M_errorTolerance;
        M_linesearch                    = DataFSI.M_linesearch;
        M_preconditioner                = DataFSI.M_preconditioner;
        M_DDNpreconditioner             = DataFSI.M_DDNpreconditioner;
        M_DDBlockPreconditioner         = DataFSI.M_DDBlockPreconditioner;
        M_method                        = DataFSI.M_method;
        M_algorithm                     = DataFSI.M_algorithm;
        M_defaultOmega                  = DataFSI.M_defaultOmega;
        M_updateEvery                   = DataFSI.M_updateEvery;
        M_fluidInterfaceFlag            = DataFSI.M_fluidInterfaceFlag;
        M_solidInterfaceFlag            = DataFSI.M_solidInterfaceFlag;
        M_structureInterfaceFlag        = DataFSI.M_structureInterfaceFlag;
        M_harmonicInterfaceFlag         = DataFSI.M_harmonicInterfaceFlag;
        M_interfaceTolerance            = DataFSI.M_interfaceTolerance;
        M_RobinNeumannCoupling          = DataFSI.M_RobinNeumannCoupling;
        M_RobinNeumannFluidCoefficient  = DataFSI.M_RobinNeumannFluidCoefficient;
        M_RobinNeumannSolidCoefficient  = DataFSI.M_RobinNeumannSolidCoefficient;
    }

	return *this;
}

template <typename Mesh>
void
DataFSI<Mesh>::setup( const GetPot& dataFile, const std::string& section )
{
    M_dataFluid->setup( dataFile );
    M_dataSolid->setup( dataFile );

    // Problem - Non Linear Richardson Parameters
    M_maxSubIterationNumber = dataFile( ( section + "/maxSubIter" ).data(), 300 );
    M_absoluteTolerance = dataFile( ( section + "/abstol" ).data(), 1.e-07 );
    M_relativeTolerance = dataFile( ( section + "/reltol" ).data(), 1.e-04 );
    M_errorTolerance = dataFile( ( section + "/etamax" ).data(), 1.e-03 );
    M_linesearch = static_cast<Int> ( dataFile( ( section + "/linesearch" ).data(), 0 ) );

    // Problem - Preconditioner
    M_preconditioner = static_cast<Preconditioner> ( dataFile( ( section + "/precond" ).data(), DIRICHLET_NEUMANN ) );
    M_DDNpreconditioner = static_cast<DDNPreconditioner> ( dataFile( ( section + "/DDNprecond" ).data(), DDN_DIRICHLET_NEUMANN ) );

    // Problem - Methods
    M_method = dataFile( ( section + "/method" ).data(), "steklovPoincare" );
    M_algorithm = dataFile( ( section + "/algorithm" ).data(), "DirichletNeumann" );

    // Problem - FixPoint / EJ
    M_defaultOmega = dataFile( ( section + "/defOmega" ).data(), 0.001);
    M_updateEvery = dataFile( ( section + "/updateEvery" ).data(), 1);

    // Interface
    M_fluidInterfaceFlag = dataFile( "interface/fluid_flag",     1 );
    M_solidInterfaceFlag = dataFile( "interface/solid_flag",     M_fluidInterfaceFlag );
    M_structureInterfaceFlag = dataFile( "interface/structure_flag", M_fluidInterfaceFlag );
    M_harmonicInterfaceFlag = dataFile( "interface/harmonic_flag",  M_fluidInterfaceFlag );
    M_interfaceTolerance = dataFile( "interface/tolerance",      0. );

    // Interface - Monolithic
    M_DDBlockPreconditioner = dataFile( "interface/DDBlockPrec", 0 );
    M_RobinNeumannCoupling  = dataFile( "interface/robinNeumannCoupling", false );
    M_RobinNeumannFluidCoefficient = dataFile( "interface/alphaf", 0.5 );
    M_RobinNeumannSolidCoefficient = dataFile( "interface/alphas", 0.5 );
}

template <typename Mesh>
bool
DataFSI<Mesh>::isMonolithic()
{
    return !( M_method.compare( "monolithic" ) && M_method.compare( "fullMonolithic" ) );
}

template <typename Mesh>
void
DataFSI<Mesh>::showMe( std::ostream& output )
{
    output << "\n*** Values for data fluid\n\n";
    M_dataFluid->showMe();

    output << "\n*** Values for data solid\n\n";
    M_dataSolid->showMe();

    output << "\n*** Values for problem\n\n";
    output << "Max subiteration number          = " << M_maxSubIterationNumber << std::endl;
    output << "Absolute tolerance               = " << M_absoluteTolerance << std::endl;
    output << "Relative tolerance               = " << M_relativeTolerance << std::endl;
    output << "Max error tolerance              = " << M_errorTolerance << std::endl;
    output << "Linesearch                       = " << M_linesearch << std::endl;

    output << "Preconditioner                   = " << M_preconditioner << std::endl;
    output << "DDNPreconditioner                = " << M_DDNpreconditioner << std::endl;
    output << "DDBlockPreconditioner            = " << M_DDBlockPreconditioner << std::endl;

    output << "Method                           = " << M_method << std::endl;
    output << "Algorithm                        = " << M_algorithm << std::endl;

    output << "Default Omega                    = " << M_defaultOmega << std::endl;
    output << "Update every                     = " << M_updateEvery << std::endl;

    output << "\n*** Values for interface\n\n";
    output << "Interface fluid                  = " << M_fluidInterfaceFlag << std::endl;
    output << "Interface solid                  = " << M_solidInterfaceFlag << std::endl;
    output << "Interface structure              = " << M_structureInterfaceFlag << std::endl;
    output << "Interface harmonic               = " << M_harmonicInterfaceFlag << std::endl;
    output << "Interface tolerance              = " << M_interfaceTolerance << std::endl;

    output << "Robin-Neumann coupling           = " << M_RobinNeumannCoupling << std::endl;
    output << "Robin-Neumann fluid coefficient  = " << M_RobinNeumannFluidCoefficient << std::endl;
    output << "Robin-Neumann solid coefficient  = " << M_RobinNeumannSolidCoefficient << std::endl;
}

} // end namespace LifeV

#endif // end DATAFSI_H
