//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 *  @brief BCInterface_Data
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 17-07-2009
 */

#ifndef BCInterface_Data_H
#define BCInterface_Data_H 1

#include <lifemc/lifesolver/BCInterface_Definitions.hpp>

namespace LifeV
{

//! BCInterface_Data - The BCInterface data container
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface_Data class provides a general container to pass information
 *  to all the BCInterface functions.
 */
class BCInterface_Data
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface_BaseList                                                       BCBaseList_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface_Data();

    //! Copy constructor
    /*!
     * @param data BCInterface_Data
     */
    BCInterface_Data( const BCInterface_Data& data );

    //! Destructor
    ~BCInterface_Data() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param data BCInterface_Data
     * @return reference to a copy of the class
     */
    BCInterface_Data& operator=( const BCInterface_Data& data );

    //@}


    //! @name Methods
    //@{
    /*!
     * @param FileName Name of the data file.
     * @param dataSection BC section
     * @param name name of the boundary condition
     */
    void ReadBC( const std::string& FileName, const std::string& dataSection, const BCName& name );

    //@}


    //! @name Set Methods
    //@{

    //! Set the name of the boundary condition
    /*!
     * @param name Boundary condition name
     */
    void SetName( const BCName& name );

    //! Set the flag of the boundary condition
    /*!
     * @param flag Boundary condition flag
     */
    void SetFlag( const BCFlag& flag );

    //! Set the type of the boundary condition
    /*!
     * @param type Boundary condition type
     */
    void SetType( const BCType& type );

    //! Set the mode of the boundary condition
    /*!
     * @param mode Boundary condition mode
     */
    void SetMode( const BCMode& mode );

    //! Set the components vector of the boundary condition
    /*!
     * @param comV Boundary condition components vector
     */
    void SetComV( const BCComV& comV );

    //! Add a component to the component vector of the boundary condition
    /*!
     * @param comV Boundary condition component
     */
    void AddComV( const UInt& comV );

    //! Set the i-component of the components vector of the boundary condition
    /*!
     * @param comV Boundary condition component
     * @param index Index value
     */
    void SetComV( const UInt& comV, const UInt& index );

    //! Set the direction string of the boundary condition
    /*!
     * @param direction Boundary condition direction string
     */
    void SetDirection( const std::string& direction );

    //! Set the base string of the boundary condition
    /*!
     * @param baseString Boundary condition base string
     */
    void SetBaseString( const std::string& baseString );

    //! Set the base type of the boundary condition
    /*!
     * @param base Boundary condition base type
     */
    void SetBase( const std::pair< std::string, BCBaseList_Type >& base );

    //@}


    //! @name Get Methods
    //@{

    //! Get the name of the boundary condition
    const BCName& GetName() const;

    //! Get the flag of the boundary condition
    const BCFlag& GetFlag() const;

    //! Get the type of the boundary condition
    const BCType& GetType() const;

    //! Get the mode of the boundary condition
    const BCMode& GetMode() const;

    //! Get the vector of components of the boundary condition
    const BCComV& GetComV() const;

    //! Get the number of components of the boundary condition
    const ID& GetComN() const;

    //! Get the direction string of the boundary condition
    const std::string& GetDirection() const;

    //! Get the base string of the boundary condition
    const std::string& GetBaseString() const;

    //! Get the base type of the boundary condition
    const std::pair< std::string, BCBaseList_Type >& GetBase() const;

    //@}

private:

    //! @name Private Methods
    //@{

    inline void ReadFlag( const std::string& FileName, const char* flag );

    inline void ReadType( const std::string& FileName, const char* type );

    inline void ReadMode( const std::string& FileName, const char* mode );

    inline void ReadComV( const std::string& FileName, const char* component );

    inline void ReadDirection( const std::string& FileName, const char* direction );

    inline void ReadBase( const std::string& FileName, const std::string& base );

    inline bool IsBase( const std::string& FileName, const char* base );

    //@}

    BCName                                    M_name;
    BCFlag                                    M_flag;
    BCType                                    M_type;
    BCMode                                    M_mode;
    BCComV                                    M_comV;
    std::string                               M_direction;
    std::string                               M_baseString;
    std::pair< std::string, BCBaseList_Type > M_base;

    // Maps
    std::map< std::string, BCType >           M_mapType;
    std::map< std::string, BCMode >           M_mapMode;
    std::map< std::string, BCBaseList_Type >  M_mapBase;
};

} // Namespace LifeV

#endif /* BCInterface_Data_H */
