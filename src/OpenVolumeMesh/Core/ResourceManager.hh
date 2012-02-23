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

#ifndef RESOURCEMANAGER_HH_
#define RESOURCEMANAGER_HH_

#include <string>
#include <vector>

#include "OpenVolumeMeshProperty.hh"
#include "PropertyHandles.hh"

namespace OpenVolumeMesh {

// Forward declarations
class BaseProperty;
template <class T>
class VertexPropertyT;
template <class T>
class EdgePropertyT;
template <class T>
class HalfEdgePropertyT;
template <class T>
class FacePropertyT;
template <class T>
class HalfFacePropertyT;
template <class T>
class CellPropertyT;
template <class T>
class MeshPropertyT;

class ResourceManager {
public:
    ResourceManager();
    virtual ~ResourceManager();

    /// Change size of stored vertex properties
    void resize_vprops(unsigned int _nv);

    /// Change size of stored edge properties
    void resize_eprops(unsigned int _ne);

    /// Change size of stored face properties
    void resize_fprops(unsigned int _nf);

    /// Change size of stored cell properties
    void resize_cprops(unsigned int _nc);

    void vertex_deleted(const VertexHandle& _h);

    void edge_deleted(const EdgeHandle& _h);

    void face_deleted(const FaceHandle& _h);

    void cell_deleted(const CellHandle& _h);

    /// Get number of vertices in mesh
    virtual unsigned int n_vertices() const = 0;
    /// Get number of edges in mesh
    virtual unsigned int n_edges() const = 0;
    /// Get number of halfedges in mesh
    virtual unsigned int n_halfedges() const = 0;
    /// Get number of faces in mesh
    virtual unsigned int n_faces() const = 0;
    /// Get number of halffaces in mesh
    virtual unsigned int n_halffaces() const = 0;
    /// Get number of cells in mesh
    virtual unsigned int n_cells() const = 0;

    template<class T> VertexPropertyT<T> request_vertex_property(const std::string& _name);

    template<class T> EdgePropertyT<T> request_edge_property(const std::string& _name);

    template<class T> HalfEdgePropertyT<T> request_halfedge_property(const std::string& _name);

    template<class T> FacePropertyT<T> request_face_property(const std::string& _name);

    template<class T> HalfFacePropertyT<T> request_halfface_property(const std::string& _name);

    template<class T> CellPropertyT<T> request_cell_property(const std::string& _name);

    template<class T> MeshPropertyT<T> request_mesh_property(const std::string& _name);

    void released_property(VertexPropHandle _handle);

    void released_property(EdgePropHandle _handle);

    void released_property(HalfEdgePropHandle _handle);

    void released_property(FacePropHandle _handle);

    void released_property(HalfFacePropHandle _handle);

    void released_property(CellPropHandle _handle);

    void released_property(MeshPropHandle _handle);

    unsigned int n_vertex_props() const { return vertex_props_.size(); }

    unsigned int n_edge_props() const { return edge_props_.size(); }

    unsigned int n_halfedge_props() const { return halfedge_props_.size(); }

    unsigned int n_face_props() const { return face_props_.size(); }

    unsigned int n_halfface_props() const { return halfface_props_.size(); }

    unsigned int n_cell_props() const { return cell_props_.size(); }

    unsigned int n_mesh_props() const { return mesh_props_.size(); }

    template<class T> void set_persistent(VertexPropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(EdgePropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(HalfEdgePropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(FacePropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(HalfFacePropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(CellPropertyT<T>& _prop, bool _flag = true);

    template<class T> void set_persistent(MeshPropertyT<T>& _prop, bool _flag = true);

private:

    std::vector<BaseProperty*> vertex_props_;

    std::vector<BaseProperty*> edge_props_;

    std::vector<BaseProperty*> halfedge_props_;

    std::vector<BaseProperty*> face_props_;

    std::vector<BaseProperty*> halfface_props_;

    std::vector<BaseProperty*> cell_props_;

    std::vector<BaseProperty*> mesh_props_;

    std::vector<BaseProperty*> persistent_vprops_;

    std::vector<BaseProperty*> persistent_eprops_;

    std::vector<BaseProperty*> persistent_heprops_;

    std::vector<BaseProperty*> persistent_fprops_;

    std::vector<BaseProperty*> persistent_hfprops_;

    std::vector<BaseProperty*> persistent_cprops_;

    std::vector<BaseProperty*> persistent_mprops_;
};

}

#if defined(INCLUDE_TEMPLATES) && !defined(RESOURCEMANAGERT_CC)
#include "ResourceManagerT.cc"
#endif

#endif /* RESOURCEMANAGER_HH_ */
