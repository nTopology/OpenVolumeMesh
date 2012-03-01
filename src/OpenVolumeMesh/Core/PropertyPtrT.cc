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

#define PROPERTYPTRT_CC

#include "PropertyPtr.hh"
#include "ResourceManager.hh"
#include "PropertyDefines.hh"

namespace OpenVolumeMesh {

template <class PropT, class HandleT>
PropertyPtr<PropT,HandleT>::PropertyPtr(PropT* _ptr, ResourceManager& _resMan, Handle _handle) :
    BaseProperty(_resMan), h_(new Holder(_ptr)), handle_(_handle) {
}

template <class PropT, class HandleT>
PropertyPtr<PropT,HandleT>::PropertyPtr(const PropertyPtr<PropT,HandleT>& _cpy) :
    BaseProperty(_cpy), h_(_cpy.h_), handle_(_cpy.handle_) {
    ++h_->count_;
}

template <class PropT, class HandleT>
PropertyPtr<PropT,HandleT>::~PropertyPtr() {
    if(--h_->count_ == 0) {
        resMan_.released_property(handle_);
        delete h_;
    }
}

template <class PropT, class HandleT>
PropertyPtr<PropT,HandleT>&
PropertyPtr<PropT,HandleT>::operator= (const PropertyPtr<PropT,HandleT>& _rhs) {

    if(--h_->count_ == 0) delete h_;

    h_ = _rhs.h_;
    ++h_->count_;

    return *this;
}

template <class PropT, class HandleT>
PropertyPtr<PropT,HandleT>&
PropertyPtr<PropT,HandleT>::operator= (PropertyPtr<PropT,HandleT>& _rhs) {

    if(--h_->count_ == 0) delete h_;

    h_ = _rhs.h_;
    ++h_->count_;

    return *this;
}

template <class PropT, class HandleT>
const typename PropT::const_reference
PropertyPtr<PropT,HandleT>::operator[](size_t _idx) const {
    assert(h_->ptr_->n_elements() > _idx);
    return (*h_->ptr_)[_idx];
}

template <class PropT, class HandleT>
typename PropT::reference
PropertyPtr<PropT,HandleT>::operator[](size_t _idx) {
    assert(h_->ptr_->n_elements() > _idx);
    return (*h_->ptr_)[_idx];
}

template <class PropT, class HandleT>
void PropertyPtr<PropT,HandleT>::resize(unsigned int _size) {
    h_->ptr_->resize(_size);
}

template <class PropT, class HandleT>
const std::string& PropertyPtr<PropT,HandleT>::name() const {
    return h_->ptr_->name();
}

template <class PropT, class HandleT>
PropT* PropertyPtr<PropT,HandleT>::operator->() {
    return h_->ptr_;
}

template <class PropT, class HandleT>
PropT& PropertyPtr<PropT,HandleT>::operator*()  {
    return *h_->ptr_;
}

template <class PropT, class HandleT>
void PropertyPtr<PropT,HandleT>::delete_element(size_t _idx) {
    h_->ptr_->delete_element(_idx);
}

template <class PropT, class HandleT>
void PropertyPtr<PropT,HandleT>::set_handle(const OpenVolumeMeshHandle& _handle) {
    handle_.idx(_handle.idx());
}

template <class PropT, class HandleT>
OpenVolumeMeshHandle PropertyPtr<PropT,HandleT>::handle() const {
    return handle_;
}

} // Namespace OpenVolumeMesh
