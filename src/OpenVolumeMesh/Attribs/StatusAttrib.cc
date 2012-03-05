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

#include "StatusAttrib.hh"

#include "../Core/TopologyKernel.hh"
#include "../Core/PropertyDefines.hh"

namespace OpenVolumeMesh {

StatusAttrib::StatusAttrib(TopologyKernel& _kernel) :
kernel_(_kernel),
v_status_(_kernel.request_vertex_property<OpenVolumeMeshStatus>("vertex_status")),
e_status_(_kernel.request_edge_property<OpenVolumeMeshStatus>("edge_status")),
he_status_(_kernel.request_halfedge_property<OpenVolumeMeshStatus>("halfedge_status")),
f_status_(_kernel.request_face_property<OpenVolumeMeshStatus>("face_status")),
hf_status_(_kernel.request_halfface_property<OpenVolumeMeshStatus>("halfface_status")),
c_status_(_kernel.request_cell_property<OpenVolumeMeshStatus>("cell_status")),
m_status_(_kernel.request_mesh_property<OpenVolumeMeshStatus>("mesh_status")) {

}

//========================================================================================

StatusAttrib::~StatusAttrib() {

}

void StatusAttrib::garbage_collection(bool _preserveManifoldness) {

    /*
     * This is not a real garbage collection in its conventional
     * sense. What happens in this routine are the following steps:
     *
     * 1. Delete all entities marked as deleted from bottom to top.
     * 2. Preserve manifoldness (optionally) by deleting all
     *    isolated entities in a top-down fashion afterwards.
     */

    for(VertexIter v_it = kernel_.vertices_begin(); v_it != kernel_.vertices_end();) {

        if(!v_status_[v_it->idx()].deleted()) {
            ++v_it;
        } else {
            v_it = kernel_.delete_vertex(*v_it);
        }
    }

    for(EdgeIter e_it = kernel_.edges_begin(); e_it != kernel_.edges_end();) {

        if(!e_status_[e_it->idx()].deleted()) {
            ++e_it;
        } else {
            e_it = kernel_.delete_edge(*e_it);
        }
    }

    for(FaceIter f_it = kernel_.faces_begin(); f_it != kernel_.faces_end();) {

        if(!f_status_[f_it->idx()].deleted()) {
            ++f_it;
        } else {
            f_it = kernel_.delete_face(*f_it);
        }
    }

    for(CellIter c_it = kernel_.cells_begin(); c_it != kernel_.cells_end();) {

        if(!c_status_[c_it->idx()].deleted()) {
            ++c_it;
        } else {
            c_it = kernel_.delete_cell(*c_it);
        }
    }

    // Step 5
    if(kernel_.has_bottom_up_adjacencies()) {
        kernel_.update_adjacencies();
    }

    // Step 6
    if(_preserveManifoldness) {
        if(kernel_.has_bottom_up_adjacencies()) {

            // Go over all faces and find those
            // that are not incident to any cell
            for(FaceIter f_it = kernel_.faces_begin(); f_it != kernel_.faces_end();) {

                // Get half-faces
                HalfFaceHandle hf0 = kernel_.halfface_handle(*f_it, 0);
                HalfFaceHandle hf1 = kernel_.halfface_handle(*f_it, 1);

                // If neither of the half-faces is incident to a cell, delete face
                if(kernel_.incident_cell(hf0) == TopologyKernel::InvalidCellHandle &&
                        kernel_.incident_cell(hf1) == TopologyKernel::InvalidCellHandle) {

                    f_it = kernel_.delete_face(*f_it);

                    kernel_.update_face_adjacencies();

                } else {
                    ++f_it;
                }
            }

            kernel_.update_edge_adjacencies();

            // Go over all edges and find those
            // whose half-edges are not incident to any half-face
            for(EdgeIter e_it = kernel_.edges_begin(); e_it != kernel_.edges_end();) {

                // Get half-edges
                HalfEdgeHandle he0 = kernel_.halfedge_handle(*e_it, 0);
                HalfEdgeHandle he1 = kernel_.halfedge_handle(*e_it, 1);

                // If neither of the half-edges is incident to a half-face, delete edge
                HalfEdgeHalfFaceIter he0hf_it = kernel_.hehf_iter(he0);
                HalfEdgeHalfFaceIter he1hf_it = kernel_.hehf_iter(he1);

                if(!he0hf_it.valid() && !he1hf_it.valid()) {

                    e_it = kernel_.delete_edge(*e_it);

                    kernel_.update_edge_adjacencies();

                } else {
                     ++e_it;
                }
            }

            // Vertex caches have to be re-computed because the face/half-face
            // indices have changed since the last deletions
            kernel_.update_vertex_adjacencies();

            // Go over all vertices and find those
            // that are not incident to any edge
            for(VertexIter v_it = kernel_.vertices_begin(); v_it != kernel_.vertices_end();) {

                // If neither of the half-edges is incident to a half-face, delete edge
                VertexOHalfEdgeIter voh_it = kernel_.voh_iter(*v_it);

                if(!voh_it.valid()) {

                    v_it = kernel_.delete_vertex(*v_it);

                    kernel_.update_vertex_adjacencies();

                } else {
                     ++v_it;
                }
            }

        } else {
            std::cerr << "Preservation of three-manifoldness in garbage_collection() "
                    << "requires bottom-up adjacencies!" << std::endl;
            return;
        }
    }
}

} // Namespace OpenVolumeMesh
