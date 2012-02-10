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

    for(std::vector<BaseProperty*>::iterator it = vertex_props_.begin();
            it != vertex_props_.end(); ++it) {
        if((*it)->persistent() && (*it)->name() == _name) {
            VertexPropertyT<T>* prop = dynamic_cast<VertexPropertyT<T>*>(*it);
            if(prop != NULL) {
                return *prop;
            }
        }
    }

    VertexPropHandle handle = vertex_props_.size();

    VertexPropertyT<T>* prop = new VertexPropertyT<T>(_name, *this, handle);
    (*prop)->resize(n_vertices());

    // Store property pointer
    vertex_props_.push_back(prop);

    // Decrease reference counter by one
    // since the property is going to be copied
    // on return
    --(*prop).h_->count_;

    return *prop;
}

template<class T>
EdgePropertyT<T> ResourceManager::request_edge_property(const std::string& _name) {

    for(std::vector<BaseProperty*>::iterator it = edge_props_.begin();
            it != edge_props_.end(); ++it) {
        if((*it)->persistent() && (*it)->name() == _name) {
            EdgePropertyT<T>* prop = dynamic_cast<EdgePropertyT<T>*>(*it);
            if(prop != NULL) {
                return *prop;
            }
        }
    }

    EdgePropHandle handle = edge_props_.size();

    EdgePropertyT<T>* prop = new EdgePropertyT<T>(_name, *this, handle);
    (*prop)->resize(n_edges());

    // Store property pointer
    edge_props_.push_back(prop);

    // Decrease reference counter by one
    // since the property is going to be copied
    // on return
    --(*prop).h_->count_;

    return *prop;
}

template<class T>
HalfEdgePropertyT<T> ResourceManager::request_halfedge_property(const std::string& _name) {

    for(std::vector<BaseProperty*>::iterator it = halfedge_props_.begin();
            it != halfedge_props_.end(); ++it) {
        if((*it)->persistent() && (*it)->name() == _name) {
            HalfEdgePropertyT<T>* prop = dynamic_cast<HalfEdgePropertyT<T>*>(*it);
            if(prop != NULL) {
                return *prop;
            }
        }
    }

    HalfEdgePropHandle handle = halfedge_props_.size();

    HalfEdgePropertyT<T>* prop = new HalfEdgePropertyT<T>(_name, *this, handle);
    (*prop)->resize(n_halfedges());

    // Store property pointer
    halfedge_props_.push_back(prop);

    // Decrease reference counter by one
    // since the property is going to be copied
    // on return
    --(*prop).h_->count_;

    return *prop;
}

template<class T>
FacePropertyT<T> ResourceManager::request_face_property(const std::string& _name) {

    for(std::vector<BaseProperty*>::iterator it = face_props_.begin();
            it != face_props_.end(); ++it) {
        if((*it)->persistent() && (*it)->name() == _name) {
            FacePropertyT<T>* prop = dynamic_cast<FacePropertyT<T>*>(*it);
            if(prop != NULL) {
                return *prop;
            }
        }
    }

    FacePropHandle handle = face_props_.size();

    FacePropertyT<T>* prop = new FacePropertyT<T>(_name, *this, handle);
    (*prop)->resize(n_faces());

    // Store property pointer
    face_props_.push_back(prop);

    // Decrease reference counter by one
    // since the property is going to be copied
    // on return
    --(*prop).h_->count_;

    return *prop;
}

template<class T>
HalfFacePropertyT<T> ResourceManager::request_halfface_property(const std::string& _name) {

    for(std::vector<BaseProperty*>::iterator it = halfface_props_.begin();
            it != halfface_props_.end(); ++it) {
        if((*it)->persistent() && (*it)->name() == _name) {
            HalfFacePropertyT<T>* prop = dynamic_cast<HalfFacePropertyT<T>*>(*it);
            if(prop != NULL) {
                return *prop;
            }
        }
    }

    HalfFacePropHandle handle = halfface_props_.size();

    HalfFacePropertyT<T>* prop = new HalfFacePropertyT<T>(_name, *this, handle);
    (*prop)->resize(n_halffaces());

    // Store property pointer
    halfface_props_.push_back(prop);

    // Decrease reference counter by one
    // since the property is going to be copied
    // on return
    --(*prop).h_->count_;

    return *prop;
}

template<class T>
CellPropertyT<T> ResourceManager::request_cell_property(const std::string& _name) {

    for(std::vector<BaseProperty*>::iterator it = cell_props_.begin();
            it != cell_props_.end(); ++it) {
        if((*it)->persistent() && (*it)->name() == _name) {
            CellPropertyT<T>* prop = dynamic_cast<CellPropertyT<T>*>(*it);
            if(prop != NULL) {
                return *prop;
            }
        }
    }

    CellPropHandle handle = cell_props_.size();

    CellPropertyT<T>* prop = new CellPropertyT<T>(_name, *this, handle);
    (*prop)->resize(n_cells());

    // Store property pointer
    cell_props_.push_back(prop);

    // Decrease reference counter by one
    // since the property is going to be copied
    // on return
    --(*prop).h_->count_;

    return *prop;
}

template<class T>
MeshPropertyT<T> ResourceManager::request_mesh_property(const std::string& _name) {

    for(std::vector<BaseProperty*>::iterator it = mesh_props_.begin();
            it != mesh_props_.end(); ++it) {
        if((*it)->persistent() && (*it)->name() == _name) {
            MeshPropertyT<T>* prop = dynamic_cast<MeshPropertyT<T>*>(*it);
            if(prop != NULL) {
                return *prop;
            }
        }
    }

    MeshPropHandle handle = mesh_props_.size();

    MeshPropertyT<T>* prop = new MeshPropertyT<T>(_name, *this, handle);
    (*prop)->resize(1u);

    // Store property pointer
    mesh_props_.push_back(prop);

    // Decrease reference counter by one
    // since the property is going to be copied
    // on return
    --(*prop).h_->count_;

    return *prop;
}

} // Namespace OpenVolumeMesh
