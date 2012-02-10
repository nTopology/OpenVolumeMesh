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

#define PROPERTYDEFINEST_CC

#include "PropertyDefines.hh"
#include "PropertyPtr.hh"

namespace OpenVolumeMesh {

/// Property classes for the different entity types
template<class T>
VertexPropertyT<T>::VertexPropertyT(const std::string& _name, ResourceManager& _resMan, VertexPropHandle _handle) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, VertexPropHandle>(new OpenVolumeMeshPropertyT<T>(_name), _resMan, _handle) {

}

template<class T>
std::ostream& VertexPropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "VProp " << std::endl;
    _ostr << PropertyPtr<OpenVolumeMeshPropertyT<T>, VertexPropHandle>::h_->ptr_->serialize(_ostr);
    return _ostr;
}

template<class T>
EdgePropertyT<T>::EdgePropertyT(const std::string& _name, ResourceManager& _resMan, EdgePropHandle _handle) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, EdgePropHandle>(new OpenVolumeMeshPropertyT<T>(_name), _resMan, _handle) {

}

template<class T>
std::ostream& EdgePropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "EProp " << std::endl;
    _ostr << PropertyPtr<OpenVolumeMeshPropertyT<T>, EdgePropHandle>::h_->ptr_->serialize(_ostr);
    return _ostr;
}

template<class T>
HalfEdgePropertyT<T>::HalfEdgePropertyT(const std::string& _name, ResourceManager& _resMan, HalfEdgePropHandle _handle) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfEdgePropHandle>(new OpenVolumeMeshPropertyT<T>(_name), _resMan, _handle) {

}

template<class T>
std::ostream& HalfEdgePropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "HEProp " << std::endl;
    _ostr << PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfEdgePropHandle>::h_->ptr_->serialize(_ostr);
    return _ostr;
}

template<class T>
FacePropertyT<T>::FacePropertyT(const std::string& _name, ResourceManager& _resMan, FacePropHandle _handle) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, FacePropHandle>(new OpenVolumeMeshPropertyT<T>(_name), _resMan, _handle) {

}

template<class T>
std::ostream& FacePropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "FProp " << std::endl;
    _ostr << PropertyPtr<OpenVolumeMeshPropertyT<T>, FacePropHandle>::h_->ptr_->serialize(_ostr);
    return _ostr;
}

template<class T>
HalfFacePropertyT<T>::HalfFacePropertyT(const std::string& _name, ResourceManager& _resMan, HalfFacePropHandle _handle) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfFacePropHandle>(new OpenVolumeMeshPropertyT<T>(_name), _resMan, _handle) {

}

template<class T>
std::ostream& HalfFacePropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "HFProp " << std::endl;
    _ostr << PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfFacePropHandle>::h_->ptr_->serialize(_ostr);
    return _ostr;
}

template<class T>
CellPropertyT<T>::CellPropertyT(const std::string& _name, ResourceManager& _resMan, CellPropHandle _handle) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, CellPropHandle>(new OpenVolumeMeshPropertyT<T>(_name), _resMan, _handle) {

}

template<class T>
std::ostream& CellPropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "CProp " << std::endl;
    _ostr << PropertyPtr<OpenVolumeMeshPropertyT<T>, CellPropHandle>::h_->ptr_->serialize(_ostr);
    return _ostr;
}

template<class T>
MeshPropertyT<T>::MeshPropertyT(const std::string& _name, ResourceManager& _resMan, MeshPropHandle _handle) :
        PropertyPtr<OpenVolumeMeshPropertyT<T>, MeshPropHandle>(new OpenVolumeMeshPropertyT<T>(_name), _resMan, _handle) {

}

template<class T>
std::ostream& MeshPropertyT<T>::serialize(std::ostream& _ostr) const {
    _ostr << "MProp " << std::endl;
    _ostr << PropertyPtr<OpenVolumeMeshPropertyT<T>, MeshPropHandle>::h_->ptr_->serialize(_ostr);
    return _ostr;
}

} // Namespace OpenVolumeMesh
