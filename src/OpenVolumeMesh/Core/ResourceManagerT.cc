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

    if(!_name.empty()) {
        for(std::vector<BaseProperty*>::iterator it = persistent_vprops_.begin();
                it != persistent_vprops_.end(); ++it) {
            if((*it)->name() == _name) {
                VertexPropertyT<T>* prop = dynamic_cast<VertexPropertyT<T>* >(*it);
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

    if(!_name.empty()) {
        for(std::vector<BaseProperty*>::iterator it = persistent_eprops_.begin();
                it != persistent_eprops_.end(); ++it) {
            if((*it)->name() == _name) {
                EdgePropertyT<T>* prop = dynamic_cast<EdgePropertyT<T>* >(*it);
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

    if(!_name.empty()) {
        for(std::vector<BaseProperty*>::iterator it = persistent_heprops_.begin();
                it != persistent_heprops_.end(); ++it) {
            if((*it)->name() == _name) {
                HalfEdgePropertyT<T>* prop = dynamic_cast<HalfEdgePropertyT<T>* >(*it);
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

    if(!_name.empty()) {
        for(std::vector<BaseProperty*>::iterator it = persistent_fprops_.begin();
                it != persistent_fprops_.end(); ++it) {
            if((*it)->name() == _name) {
                FacePropertyT<T>* prop = dynamic_cast<FacePropertyT<T>* >(*it);
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

    if(!_name.empty()) {
        for(std::vector<BaseProperty*>::iterator it = persistent_hfprops_.begin();
                it != persistent_hfprops_.end(); ++it) {
            if((*it)->name() == _name) {
                HalfFacePropertyT<T>* prop = dynamic_cast<HalfFacePropertyT<T>* >(*it);
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

    if(!_name.empty()) {
        for(std::vector<BaseProperty*>::iterator it = persistent_cprops_.begin();
                it != persistent_cprops_.end(); ++it) {
            if((*it)->name() == _name) {
                CellPropertyT<T>* prop = dynamic_cast<CellPropertyT<T>* >(*it);
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

    if(!_name.empty()) {
        for(std::vector<BaseProperty*>::iterator it = persistent_mprops_.begin();
                it != persistent_mprops_.end(); ++it) {
            if((*it)->name() == _name) {
                MeshPropertyT<T>* prop = dynamic_cast<MeshPropertyT<T>* >(*it);
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

template<class T>
void ResourceManager::set_persistent(VertexPropertyT<T>& _prop, bool _flag) {

    if(_flag == _prop.h_->persistent_) return;

    if(!_flag) {
        for(std::vector<BaseProperty*>::iterator it = persistent_vprops_.begin();
                it != persistent_vprops_.end(); ++it) {
            if((*it)->ptr() == _prop.h_->ptr_) {
                _prop.h_->persistent_ = false;
                delete (*it);
                persistent_vprops_.erase(it);
                break;
            }
        }
    } else {
        _prop.h_->persistent_ = true;
        VertexPropertyT<T>* prop = new VertexPropertyT<T>(_prop);
        persistent_vprops_.push_back(prop);
    }
}

template<class T>
void ResourceManager::set_persistent(EdgePropertyT<T>& _prop, bool _flag) {

    if(_flag == _prop.h_->persistent_) return;

    if(!_flag) {
        for(std::vector<BaseProperty*>::iterator it = persistent_eprops_.begin();
                it != persistent_eprops_.end(); ++it) {
            if((*it)->ptr() == _prop.h_->ptr_) {
                _prop.h_->persistent_ = false;
                delete (*it);
                persistent_eprops_.erase(it);
                break;
            }
        }
    } else {
        _prop.h_->persistent_ = true;
        EdgePropertyT<T>* prop = new EdgePropertyT<T>(_prop);
        persistent_eprops_.push_back(prop);
    }
}

template<class T>
void ResourceManager::set_persistent(HalfEdgePropertyT<T>& _prop, bool _flag) {

    if(_flag == _prop.h_->persistent_) return;

    if(!_flag) {
        for(std::vector<BaseProperty*>::iterator it = persistent_heprops_.begin();
                it != persistent_heprops_.end(); ++it) {
            if((*it)->ptr() == _prop.h_->ptr_) {
                _prop.h_->persistent_ = false;
                delete (*it);
                persistent_heprops_.erase(it);
                break;
            }
        }
    } else {
        _prop.h_->persistent_ = true;
        HalfEdgePropertyT<T>* prop = new HalfEdgePropertyT<T>(_prop);
        persistent_heprops_.push_back(prop);
    }
}

template<class T>
void ResourceManager::set_persistent(FacePropertyT<T>& _prop, bool _flag) {

    if(_flag == _prop.h_->persistent_) return;

    if(!_flag) {
        for(std::vector<BaseProperty*>::iterator it = persistent_fprops_.begin();
                it != persistent_fprops_.end(); ++it) {
            if((*it)->ptr() == _prop.h_->ptr_) {
                _prop.h_->persistent_ = false;
                delete (*it);
                persistent_fprops_.erase(it);
                break;
            }
        }
    } else {
        _prop.h_->persistent_ = true;
        FacePropertyT<T>* prop = new FacePropertyT<T>(_prop);
        persistent_fprops_.push_back(prop);
    }
}

template<class T>
void ResourceManager::set_persistent(HalfFacePropertyT<T>& _prop, bool _flag) {

    if(_flag == _prop.h_->persistent_) return;

    if(!_flag) {
        for(std::vector<BaseProperty*>::iterator it = persistent_hfprops_.begin();
                it != persistent_hfprops_.end(); ++it) {
            if((*it)->ptr() == _prop.h_->ptr_) {
                _prop.h_->persistent_ = false;
                delete (*it);
                persistent_hfprops_.erase(it);
                break;
            }
        }
    } else {
        _prop.h_->persistent_ = true;
        HalfFacePropertyT<T>* prop = new HalfFacePropertyT<T>(_prop);
        persistent_hfprops_.push_back(prop);
    }
}

template<class T>
void ResourceManager::set_persistent(CellPropertyT<T>& _prop, bool _flag) {

    if(_flag == _prop.h_->persistent_) return;

    if(!_flag) {
        for(std::vector<BaseProperty*>::iterator it = persistent_cprops_.begin();
                it != persistent_cprops_.end(); ++it) {
            if((*it)->ptr() == _prop.h_->ptr_) {
                _prop.h_->persistent_ = false;
                delete (*it);
                persistent_cprops_.erase(it);
                break;
            }
        }
    } else {
        _prop.h_->persistent_ = true;
        CellPropertyT<T>* prop = new CellPropertyT<T>(_prop);
        persistent_cprops_.push_back(prop);
    }
}

template<class T>
void ResourceManager::set_persistent(MeshPropertyT<T>& _prop, bool _flag) {

    if(_flag == _prop.h_->persistent_) return;

    if(!_flag) {
        for(std::vector<BaseProperty*>::iterator it = persistent_mprops_.begin();
                it != persistent_mprops_.end(); ++it) {
            if((*it)->ptr() == _prop.h_->ptr_) {
                _prop.h_->persistent_ = false;
                delete (*it);
                persistent_mprops_.erase(it);
                break;
            }
        }
    } else {
        _prop.h_->persistent_ = true;
        MeshPropertyT<T>* prop = new MeshPropertyT<T>(_prop);
        persistent_mprops_.push_back(prop);
    }
}

} // Namespace OpenVolumeMesh
