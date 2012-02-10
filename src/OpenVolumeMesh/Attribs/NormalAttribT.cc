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

#define NORMALATTRIBT_CC

#include <set>

#include "NormalAttrib.hh"

#include "../Core/GeometryKernel.hh"

namespace OpenVolumeMesh {

template <class VecT>
NormalAttrib<VecT>::NormalAttrib(GeometryKernel<VecT>& _kernel) :
kernel_(_kernel),
v_normals_(_kernel.template request_vertex_property<VecT>("vertex_normals")),
f_normals_(_kernel.template request_face_property<VecT>("face_normals")) {

}

template <class VecT>
NormalAttrib<VecT>::~NormalAttrib() {

}

template <class VecT>
void NormalAttrib<VecT>::update_vertex_normals() {

    if(!kernel_.has_bottom_up_adjacencies()) {
        std::cerr << "Error: update_vertex_normals() needs bottom-up adjacencies!" << std::endl;
        return;
    }

    // Compute face normals
    update_face_normals();

    for(VertexIter v_it = kernel_.v_iter(); v_it.valid(); ++v_it) {
        compute_vertex_normal(*v_it);
    }
}

template <class VecT>
void NormalAttrib<VecT>::update_face_normals() {

    for(FaceIter f_it = kernel_.f_iter(); f_it.valid(); ++f_it) {
        // Assume the face is planar, so just take the
        // first two edges
        compute_face_normal(*f_it);
    }
}

template <typename VecT>
void NormalAttrib<VecT>::compute_vertex_normal(const VertexHandle& _vh) {

    std::set<FaceHandle> faces;
    for(VertexOHalfEdgeIter voh_it = kernel_.voh_iter(_vh);
            voh_it.valid(); ++voh_it) {

        for(HalfEdgeHalfFaceIter hehf_it = kernel_.hehf_iter(*voh_it);
                hehf_it.valid(); ++hehf_it) {
            if(kernel_.is_boundary(*hehf_it)) {
                faces.insert(kernel_.face_handle(*hehf_it));
            }
        }
    }
    VecT normal;
    for(std::set<FaceHandle>::const_iterator f_it = faces.begin();
            f_it != faces.end(); ++f_it) {
        normal += f_normals_[f_it->idx()];
    }

    normal.normalize();

    v_normals_[_vh.idx()] = normal;
}

template <typename VecT>
void NormalAttrib<VecT>::compute_face_normal(const FaceHandle& _fh) {

    if(kernel_.face(_fh).halfedges().size() < 3) {
        std::cerr << "Warning: Degenerate face detected!" << std::endl;
        return;
    }

    // Always try to compute the outside normals
    HalfFaceHandle hfh = kernel_.is_boundary(kernel_.halfface_handle(_fh, 0)) ?
            kernel_.halfface_handle(_fh, 0) : kernel_.halfface_handle(_fh, 1);

    std::vector<HalfEdgeHandle> halfedges = kernel_.halfface(hfh).halfedges();
    std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();

    VecT p1 = kernel_.vertex(kernel_.halfedge(*he_it).from_vertex());
    VecT p2 = kernel_.vertex(kernel_.halfedge(*he_it).to_vertex());
    ++he_it;
    VecT p3 = kernel_.vertex(kernel_.halfedge(*he_it).to_vertex());

    VecT n = (p3 - p2) % (p1 - p2);
    n.normalize();

    f_normals_[_fh.idx()] = n;
}

} // Namespace OpenVolumeMesh
