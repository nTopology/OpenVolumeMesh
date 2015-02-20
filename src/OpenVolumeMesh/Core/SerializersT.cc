/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 236 $                                                         *
 *   $Date: 2013-02-19 12:32:33 +0100 (Tue, 19 Feb 2013) $                    *
 *   $LastChangedBy: kremer $                                                *
 *                                                                           *
\*===========================================================================*/

#define SERIALIZERST_CC

#include "Serializers.hh"


namespace OpenVolumeMesh
{


template <typename ValueT>
std::ostream& serialize(std::ostream& _ostr, const ValueT& _rhs)
{
    _ostr << _rhs;
    return _ostr;
}


template <typename ValueT>
std::istream& deserialize(std::istream& _istr, ValueT& _rhs)
{
    _istr >> _rhs;
    return _istr;
}


template <typename KeyT, typename ValueT>
std::ostream& operator<<(std::ostream& os, const std::map< KeyT, ValueT >& rhs)
{
    os << rhs.size() << std::endl;
    for (typename std::map< KeyT, ValueT >::const_iterator it = rhs.begin();
         it != rhs.end();
         ++it)
    {
        serialize(os,it->first) << std::endl;
        serialize(os, it->second) << std::endl;
    }

    return os;
}

template <typename KeyT, typename ValueT>
std::istream& operator>>(std::istream& is, std::map< KeyT, ValueT >& rhs)
{

    size_t size;
    is >> size;
    rhs.clear();
    for (size_t i=0; i<size; i++)
    {
        KeyT key;
        ValueT value;
        deserialize(is, key);
        deserialize(is, value);
        rhs[key] = value;
    }

    return is;
}

template <typename ValueT>
std::ostream& operator<<(std::ostream& _ostr, const std::vector< ValueT >& _rhs)
{
    _ostr << _rhs.size() << std::endl;
    for (size_t i = 0; i < _rhs.size(); ++i)
        serialize(_ostr, _rhs[i]) << std::endl;
    return _ostr;
}

template <typename ValueT>
std::istream& operator>>(std::istream& _istr, std::vector< ValueT >& _rhs)
{
    size_t size;
    _istr >> size;
    _rhs.resize(size);
    for (size_t i=0; i<size; i++)
        deserialize(_istr,_rhs[i]);

    return _istr;
}


}
