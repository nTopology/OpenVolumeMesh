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

#ifndef PROPERTYDEFINES_HH_
#define PROPERTYDEFINES_HH_

#include <iostream>
#include <typeinfo>

#include "BaseProperty.hh"
#include "PropertyHandles.hh"

namespace OpenVolumeMesh {

template <class T>
class OpenVolumeMeshPropertyT;
template <class PropT, class HandleT>
class PropertyPtr;

class ResourceManager;

/// Property classes for the different entity types
template<class T>
class VertexPropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, VertexPropHandle> {
public:
    VertexPropertyT(const std::string& _name, ResourceManager& _resMan, VertexPropHandle _handle);
    ~VertexPropertyT() {}
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
};
template<class T>
class EdgePropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, EdgePropHandle> {
public:
    EdgePropertyT(const std::string& _name, ResourceManager& _resMan, EdgePropHandle _handle);
    ~EdgePropertyT() {}
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
};
template<class T>
class HalfEdgePropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfEdgePropHandle> {
public:
    HalfEdgePropertyT(const std::string& _name, ResourceManager& _resMan, HalfEdgePropHandle _handle);
    ~HalfEdgePropertyT() {}
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
};
template<class T>
class FacePropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, FacePropHandle> {
public:
    FacePropertyT(const std::string& _name, ResourceManager& _resMan, FacePropHandle _handle);
    ~FacePropertyT() {}
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
};
template<class T>
class HalfFacePropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfFacePropHandle> {
public:
    HalfFacePropertyT(const std::string& _name, ResourceManager& _resMan, HalfFacePropHandle _handle);
    ~HalfFacePropertyT() {}
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
};
template<class T>
class CellPropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, CellPropHandle> {
public:
    CellPropertyT(const std::string& _name, ResourceManager& _resMan, CellPropHandle _handle);
    ~CellPropertyT() {}
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
};
template<class T>
class MeshPropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, MeshPropHandle> {
public:
    MeshPropertyT(const std::string& _name, ResourceManager& _resMan, MeshPropHandle _handle);
    ~MeshPropertyT() {}
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
};

template <class T>
const std::string typeName() { return typeid(T).name(); }

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(PROPERTYDEFINEST_CC)
#include "PropertyDefinesT.cc"
#endif

#endif /* PROPERTYDEFINES_HH_ */
