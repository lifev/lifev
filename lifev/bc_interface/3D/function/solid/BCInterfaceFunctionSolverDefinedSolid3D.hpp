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
 *  @file
 *  @brief File containing the BCInterfaceFunctionSolverDefined class and specializations
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionSolverDefinedSolid3D_H
#define BCInterfaceFunctionSolverDefinedSolid3D_H 1

// Structural solver includes
#include <lifev/structure/solver/StructuralOperator.hpp>

// BCInterface includes
#include <lifev/bc_interface/3D/bc/BCInterfaceData3D.hpp>
#include <lifev/bc_interface/core/function/BCInterfaceFunctionSolverDefined.hpp>

namespace LifeV
{

//! BCInterfaceFunctionSolverDefined - Template specialization of \c BCInterfaceFunctionSolverDefined for Solid 3D problems
/*!
 *  @author Cristiano Malossi
 *
 *
 *  The BCInterfaceFunctionSolverDefined class provides the interface between the
 *  \c BCInterface3D and the solver defined boundary conditions of the \c StructuralOperator.
 *
 *  <b>DETAILS:</b> <BR>
 *  The constructor of the class takes a string contains the ID of the boundary condition to impose.
 *  The list of available conditions is the \c Solid3DFunction enum. These are:
 *
 *  <ol>
 *      <li> RobinWall
 *  </ol>
 */
template< >
class BCInterfaceFunctionSolverDefined< BCHandler, StructuralOperator<RegionMesh <LinearTetra> > >
{
public:

    //! @name Type definitions
    //@{

    typedef BCHandler                                                              bcHandler_Type;
    typedef boost::shared_ptr< bcHandler_Type >                                    bcHandlerPtr_Type;

    typedef StructuralOperator< RegionMesh <LinearTetra> >                         physicalSolver_Type;
    typedef boost::shared_ptr< physicalSolver_Type >                               physicalSolverPtr_Type;

    typedef BCInterfaceFactory< bcHandler_Type, physicalSolver_Type >              factory_Type;
    typedef BCInterfaceFunction< bcHandler_Type, physicalSolver_Type >             bcFunction_Type;
    typedef boost::shared_ptr< bcFunction_Type >                                   bcFunctionPtr_Type;
    typedef std::vector< bcFunctionPtr_Type >                                      vectorFunction_Type;

    typedef BCInterfaceFunctionParserSolver< bcHandler_Type, physicalSolver_Type > functionParserSolver_Type;
    typedef boost::shared_ptr< functionParserSolver_Type >                         functionParserSolverPtr_Type;

    typedef BCInterfaceData3D                                                      data_Type;
    typedef boost::shared_ptr< data_Type >                                         dataPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFunctionSolverDefined();

    //! Destructor
    virtual ~BCInterfaceFunctionSolverDefined() {}

    //@}


    //! @name Methods
    //@{

    //! Copy the stored parameters in the data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData ( dataPtr_Type& data );

    //! Assign a boundary function to the boundary condition vector base
    /*!
     * @param physicalSolver FSI physical solver,
     * @param base boundary condition base
     */
    template< class BCBaseType >
    void assignFunction ( BCBaseType& base )
    {
        checkFunction ( base );
    }

    //! Update the solver variables
    void updatePhysicalSolverVariables();

    //@}


    //! @name Set methods
    //@{

    //! Set data
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData ( const dataPtr_Type& data );

    //! Set the physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver ( const physicalSolverPtr_Type& physicalSolver )
    {
        M_physicalSolver = physicalSolver;
    }

    //@}


    //! @name Get methods
    //@{

    //! Detect the correct base type
    /*!
     * @param bcBaseType the type of the base
     */
    baseContainer_Type baseType() const;

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionSolverDefined ( const BCInterfaceFunctionSolverDefined& function );

    BCInterfaceFunctionSolverDefined& operator= ( const BCInterfaceFunctionSolverDefined& function );

    //@}


    //! @name Private Methods
    //@{

    void checkFunction ( BCVector& base );

    void checkFunction ( BCFunctionBase& base );

    void checkFunction ( BCVectorInterface& base );

    //@}

    enum Solid3DFunction
    {
        RobinWall
    };

    Solid3DFunction                                M_solid3DFunction;

    physicalSolverPtr_Type                         M_physicalSolver;

    // The following members are required since the Solid BC are applied
    // a posteriori, when setPhysicalSolver() is called.

    // Classical parameters
    bcName_Type                                    M_name;
    bcFlag_Type                                    M_flag;
    bcType_Type                                    M_type;
    bcMode_Type                                    M_mode;
    bcComponentsVec_Type                           M_componentsVector;

    // RobinViscoelastic
    vectorFunction_Type                            M_vectorFunctionRobin;
    physicalSolver_Type::vectorPtr_Type            M_robinRHS;
    physicalSolver_Type::vectorPtr_Type            M_robinAlphaCoefficient;
    physicalSolver_Type::vectorPtr_Type            M_robinBetaCoefficient;
};

} // Namespace LifeV

#endif /* BCInterfaceFunctionSolverDefinedSolid3D_H */
