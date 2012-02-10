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

#ifndef PROPERTYPTR_HH_
#define PROPERTYPTR_HH_

#include <iostream>
#include <string>
#include <cassert>

#include "BaseProperty.hh"
#include "PropertyHandles.hh"

namespace OpenVolumeMesh {

class ResourceManager;

/**
 * \class PropertyPtr
 *
 * A smart-pointer-like class that counts the encapsulated
 * object's references and automatically deletes the memory
 * as soon as the object is not used anymore.
 */

template <class PropT, class HandleT>
class PropertyPtr : public BaseProperty {
protected:

    /**
     * \class Holder
     * A helper class that encapsulates a pointer
     * to an object and counts its number of references.
     */
    class Holder {
    public:
        Holder(PropT* _ptr) : ptr_(_ptr), count_(1u) {}
        ~Holder() {
            if(ptr_ != NULL) {
                delete ptr_;
                ptr_ = NULL;
            }
        }

        PropT* ptr_;
        unsigned int count_;
    };

    Holder* h_;

public:

    friend class ResourceManager;

    typedef HandleT Handle;

    /// Constructor
    PropertyPtr(PropT* _ptr, ResourceManager& _resMan, Handle _handle);

    /// Copy Constructor
    PropertyPtr(const PropertyPtr<PropT,HandleT>& _cpy);

    /// Destructor
    virtual ~PropertyPtr();

    /// Assignment operator with const reference
    PropertyPtr<PropT,HandleT>& operator= (const PropertyPtr<PropT,HandleT>& _rhs);

    /// Assignment operator with non-const reference
    PropertyPtr<PropT,HandleT>& operator= (PropertyPtr<PropT,HandleT>& _rhs);

    const typename PropT::const_reference operator[](size_t _idx) const;

    typename PropT::reference operator[](size_t _idx);

    virtual void resize(unsigned int _size);

    virtual const std::string& name() const;

    /// Access the encapsulated pointer
    PropT* operator->();

    /// Access the encapsulated dereferenced pointer
    PropT& operator*();

    virtual void set_persistent();

    virtual bool persistent() const;

    virtual void delete_element(size_t _idx);

protected:

    virtual void set_non_persistent();

    virtual void set_ref_count(unsigned int _c);

    virtual unsigned int ref_count() const;

    virtual void set_handle(const OpenVolumeMeshHandle& _handle);

    virtual OpenVolumeMeshHandle handle() const;

private:

     Handle handle_;
};

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(PROPERTYPTRT_CC)
#include "PropertyPtrT.cc"
#endif

#endif /* PROPERTYPTR_HH_ */
