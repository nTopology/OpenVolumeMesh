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
 *   $Revision$                                                         *
 *   $Date$                    *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef BASEPROPERTY_HH_
#define BASEPROPERTY_HH_

#include <string>

#include "OpenVolumeMeshHandle.hh"

namespace OpenVolumeMesh {

class ResourceManager;

class BaseProperty {
public:
    friend class ResourceManager;

    BaseProperty(ResourceManager& _resMan) : resMan_(_resMan) {}

    BaseProperty(const BaseProperty& _cpy) : resMan_(_cpy.resMan_) {}

    virtual ~BaseProperty() {}

    virtual BaseProperty& operator= (const BaseProperty& _rhs);

    virtual BaseProperty& operator= (BaseProperty& _rhs);

    virtual void resize(unsigned int /*_size*/) = 0;

    virtual const std::string& name() const = 0;

    virtual void delete_element(size_t _idx) = 0;

    virtual std::ostream& serialize(std::ostream& _ostr) const = 0;

    virtual std::istream& deserialize(std::istream& _istr) = 0;

protected:

    virtual void set_handle(const OpenVolumeMeshHandle& /*_handle*/) = 0;

    virtual OpenVolumeMeshHandle handle() const = 0;

    ResourceManager& resMan_;
};

} // Namespace OpenVolumeMesh

#endif /* BASEPROPERTY_HH_ */
