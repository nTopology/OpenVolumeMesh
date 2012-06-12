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

#define RESOURCEMANAGERT_CC

#include "ResourceManager.hh"

#include "PropertyDefines.hh"

#include "PropertyPtr.hh"

namespace OpenVolumeMesh {

template<class T>
VertexPropertyT<T> ResourceManager::request_vertex_property(const std::string& _name) {

    return request_property<std::vector<BaseProperty*>,VertexPropertyT<T>,VertexPropHandle>(vertex_props_, _name, n_vertices());
}

template<class T>
EdgePropertyT<T> ResourceManager::request_edge_property(const std::string& _name) {

    return request_property<std::vector<BaseProperty*>,EdgePropertyT<T>,EdgePropHandle>(edge_props_, _name, n_edges());
}

template<class T>
HalfEdgePropertyT<T> ResourceManager::request_halfedge_property(const std::string& _name) {

    return request_property<std::vector<BaseProperty*>,HalfEdgePropertyT<T>,HalfEdgePropHandle>(halfedge_props_, _name, n_edges()*2u);
}

template<class T>
FacePropertyT<T> ResourceManager::request_face_property(const std::string& _name) {

    return request_property<std::vector<BaseProperty*>,FacePropertyT<T>,FacePropHandle>(face_props_, _name, n_faces());
}

template<class T>
HalfFacePropertyT<T> ResourceManager::request_halfface_property(const std::string& _name) {

    return request_property<std::vector<BaseProperty*>,HalfFacePropertyT<T>,HalfFacePropHandle>(halfface_props_, _name, n_faces()*2u);
}

template<class T>
CellPropertyT<T> ResourceManager::request_cell_property(const std::string& _name) {

    return request_property<std::vector<BaseProperty*>,CellPropertyT<T>,CellPropHandle>(cell_props_, _name, n_cells());
}

template<class T>
MeshPropertyT<T> ResourceManager::request_mesh_property(const std::string& _name) {

    return request_property<std::vector<BaseProperty*>,MeshPropertyT<T>,MeshPropHandle>(mesh_props_, _name, 1);
}

template<class StdVecT, class PropT, class HandleT>
PropT ResourceManager::request_property(StdVecT& _vec, const std::string& _name, size_t _size) {

    if(!_name.empty()) {
        for(typename StdVecT::iterator it = _vec.begin();
                it != _vec.end(); ++it) {
            if((*it)->name() == _name) {
                PropT* prop = dynamic_cast<PropT*>(*it);
                if(prop != NULL) return *prop;
                else break;
            }
        }
    }

    HandleT handle(_vec.size());

    PropT* prop = new PropT(_name, *this, handle);
    prop->resize(_size);

    // Store property pointer
    _vec.push_back(prop);

    return *prop;
}

template<class T>
void ResourceManager::set_persistent(VertexPropertyT<T>& _prop, bool _flag) {

    set_persistentT(_prop, _flag);
}

template<class T>
void ResourceManager::set_persistent(EdgePropertyT<T>& _prop, bool _flag) {

    set_persistentT(_prop, _flag);
}

template<class T>
void ResourceManager::set_persistent(HalfEdgePropertyT<T>& _prop, bool _flag) {

    set_persistentT(_prop, _flag);
}

template<class T>
void ResourceManager::set_persistent(FacePropertyT<T>& _prop, bool _flag) {

    set_persistentT(_prop, _flag);
}

template<class T>
void ResourceManager::set_persistent(HalfFacePropertyT<T>& _prop, bool _flag) {

    set_persistentT(_prop, _flag);
}

template<class T>
void ResourceManager::set_persistent(CellPropertyT<T>& _prop, bool _flag) {

    set_persistentT(_prop, _flag);
}

template<class T>
void ResourceManager::set_persistent(MeshPropertyT<T>& _prop, bool _flag) {

    set_persistentT(_prop, _flag);
}

template<class PropT>
void ResourceManager::set_persistentT(PropT& _prop, bool _flag) {

    if(_flag == _prop->persistent()) return;

    _prop->set_persistent(_flag);
}

template<class StdVecT>
void ResourceManager::remove_property(StdVecT& _vec, size_t _idx) {

    (*(_vec.begin() + _idx))->lock();
    delete *(_vec.begin() + _idx);
    _vec.erase(_vec.begin() + _idx);
    size_t n = _vec.size();
    for(size_t i = 0; i < n; ++i) {
        _vec[i]->set_handle(OpenVolumeMeshHandle(i));
    }
}

template<class StdVecT>
void ResourceManager::resize_props(StdVecT& _vec, unsigned int _n) {

    for(typename StdVecT::iterator it = _vec.begin();
            it != _vec.end(); ++it) {
        (*it)->resize(_n);
    }
}

template<class StdVecT>
void ResourceManager::entity_deleted(StdVecT& _vec, const OpenVolumeMeshHandle& _h) {

    for(typename StdVecT::iterator it = _vec.begin();
            it != _vec.end(); ++it) {
        (*it)->delete_element(_h.idx());
    }
}

template<class StdVecT>
void ResourceManager::clearVec(StdVecT& _vec) {

    for(typename StdVecT::iterator it = _vec.begin();
            it != _vec.end(); ++it) {
        if(!(*it)->persistent()) {
            std::cerr << "Could not clear properties since at " <<
                    "least one property is still in use!" << std::endl;
            return;
        }
    }

    for(typename StdVecT::iterator it = _vec.begin();
            it != _vec.end(); ++it) {
        delete *it;
    }
    _vec.clear();
}

} // Namespace OpenVolumeMesh
