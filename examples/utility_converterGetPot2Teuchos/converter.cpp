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
    @brief GetPot2Teuchos converter

    @author Umberto Villa <uvilla@emory.edu>
    @contributor
    @maintainer

    @date 04-10-2010

	Implementation.
 */
#include<boost/algorithm/string/split.hpp>
#include<boost/algorithm/string/classification.hpp>
#include<ctype.h>

#include"converter.hpp"


namespace LifeV
{

bool fillFromGetPot( const std::string & dataFileName, Teuchos::ParameterList & _pList)
{
    typedef STRING_VECTOR::iterator Iterator;

    // Two possible scenarios are possible:
    // 1) All processor load the dataFile, and we don't communicate nothing
    // 2) Only processor 0 load the dataFile, generate _pList and finally broadcast it.
    // I will implement 1 for simplicity
    GetPot dataFile(dataFileName);

    STRING_VECTOR varNames(dataFile.get_variable_names());
    STRING_VECTOR sectNames(dataFile.get_section_names());

    for (Iterator it=varNames.begin(); it != varNames.end(); ++it)
    {
        STRING_VECTOR subdirs;
        boost::split( subdirs, *it, boost::is_any_of( "/" ) );
        Teuchos::ParameterList* pList_ptr(&_pList);
        for (STRING_VECTOR::size_type i(0); i<subdirs.size() - 1; ++i)
            pList_ptr = &(pList_ptr->sublist(subdirs[i]));
        std::string stringValue = dataFile(it->c_str(), "WILL_NEVER_HAPPEN");
        GetPotVariable mytype(guessMyType(stringValue));
        switch (mytype)
        {
        case IntegerVariable:
        {
            int var;
            from_string(var, stringValue);
            pList_ptr->set(subdirs[subdirs.size()-1], var);
            break;
        }
        case RealVariable:
        {
            Real var;
            from_string(var, stringValue);
            pList_ptr->set(subdirs[subdirs.size()-1], var);
            break;
        }
        case BooleanVariable:
        {
            bool var;

            if (strcmp(stringValue.c_str(),"false")==0)
                var = false;
            else
            {
                if (strcmp(stringValue.c_str(),"true")==0)
                    var = true;
                else
                    exit(1);
            }
            pList_ptr->set(subdirs[subdirs.size()-1], var);
            break;
        }
        case StringVariable:
            pList_ptr->set(subdirs[subdirs.size()-1], stringValue);
            break;
        }
    }

    return true;
}

GetPotVariable guessMyType(std::string & stringValue)
{

    //remove semicolon if present
    if (*stringValue.rbegin() == ';')
        stringValue.erase(stringValue.end() - 1);

    typedef std::string::iterator StringIterator;
    typedef std::vector<StringIterator> Container;

    Container notDigitIndex;

    //store the non numeric digits
    for (std::string::iterator it=stringValue.begin(); it!=stringValue.end(); ++it)
        if (!isdigit(*it))
            notDigitIndex.push_back(it);

    //IntegerVariable: n, -n
    if (notDigitIndex.size()==0)
        return IntegerVariable;

    if (notDigitIndex.size() == 1 &&  *(notDigitIndex[0]) =='-')
        return IntegerVariable;

    //RealVariable: 1.0, -1.0, 1e5, -1e5, 1.e5, -1.e5, 1e-5, -1e-5, 1.e-5, -1.e-5
    if (notDigitIndex.size()<=4)
    {
        int nmenosplus(0), npoints(0), ne(0);
        bool isReal(true);

        for (Container::iterator it=notDigitIndex.begin(); it != notDigitIndex.end(); ++it)
        {
            switch (**it)
            {
            case '-':
                ++nmenosplus;
                break;
            case '+':
                ++nmenosplus;
                break;
            case 'e':
                ++npoints;
                break;
            case '.':
                ++ne;
                break;
            default:
                isReal = false;
                break;
            }
        }
        if (isReal && nmenosplus < 3 && npoints+ne > 0 && npoints < 2 && ne < 2)
            return RealVariable;
    }


    //It's a string
    //1) put it to lower case:
    for (std::string::iterator it=stringValue.begin(); it!=stringValue.end(); ++it)
        *it = tolower(*it);

    // check if is a "false"
    if (strcmp("false", stringValue.c_str())==0 || strcmp("f", stringValue.c_str())==0)
    {
        stringValue = "false";
        return BooleanVariable;
    }

    // check if is a "true"
    if (strcmp("true", stringValue.c_str())==0 || strcmp("t", stringValue.c_str())==0)
    {
        stringValue = "true";
        return BooleanVariable;
    }

    // it's a normal string
    return StringVariable;

}

}
