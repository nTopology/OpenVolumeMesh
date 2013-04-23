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
 *   $Revision: 36 $                                                         *
 *   $Date: 2012-01-10 18:00:06 +0100 (Di, 10 Jan 2012) $                    *
 *   $LastChangedBy: kremer $                                                *
 *                                                                           *
\*===========================================================================*/

#define COLORATTRIBT_CC

#include "ColorAttrib.hh"

namespace OpenVolumeMesh {

template <class ColT>
ColorAttrib<ColT>::ColorAttrib(TopologyKernel& _kernel, const ColT _def) :
        vcolor_prop_(_kernel.request_vertex_property<ColT>("vertex_color", _def)),
        ecolor_prop_(_kernel.request_edge_property<ColT>("edge_color", _def)),
        hecolor_prop_(_kernel.request_halfedge_property<ColT>("halfedge_color", _def)),
        fcolor_prop_(_kernel.request_face_property<ColT>("face_color", _def)),
        hfcolor_prop_(_kernel.request_halfface_property<ColT>("halfface_color", _def)),
        ccolor_prop_(_kernel.request_cell_property<ColT>("cell_color", _def)),
        kernel_(_kernel),
        vertex_colors_available_(false),
        halfedge_colors_available_(false),
        edge_colors_available_(false),
        halfface_colors_available_(false),
        face_colors_available_(false),
        cell_colors_available_(false),
        default_color_ (_def)
{

}

template <class ColT>
ColorAttrib<ColT>::~ColorAttrib() {

}

template <class ColT>
void ColorAttrib<ColT>::clear_vertex_colors()
{
    for (VertexIter v_it = kernel_.vertices_begin(); v_it != kernel_.vertices_end(); ++v_it)
        vcolor_prop_[v_it->idx()] = default_color_;
    vertex_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_halfedge_colors()
{
    for (HalfEdgeIter he_it = kernel_.halfedges_begin(); he_it != kernel_.halfedges_end(); ++he_it)
        hecolor_prop_[he_it->idx()] = default_color_;
    halfedge_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_edge_colors()
{
    for (EdgeIter e_it = kernel_.edges_begin(); e_it != kernel_.edges_end(); ++e_it)
        ecolor_prop_[e_it->idx()] = default_color_;
    edge_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_halfface_colors()
{
    for (HalfFaceIter hf_it = kernel_.halffaces_begin(); hf_it != kernel_.halffaces_end(); ++hf_it)
        hfcolor_prop_[hf_it->idx()] = default_color_;
    halfface_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_face_colors()
{
    for (FaceIter f_it = kernel_.faces_begin(); f_it != kernel_.faces_end(); ++f_it)
        fcolor_prop_[f_it->idx()] = default_color_;
    face_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_cell_colors()
{
    for (CellIter c_it = kernel_.cells_begin(); c_it != kernel_.cells_end(); ++c_it)
        ccolor_prop_[c_it->idx()] = default_color_;
    cell_colors_available_   = false;
}

} // Namespace OpenVolumeMesh
