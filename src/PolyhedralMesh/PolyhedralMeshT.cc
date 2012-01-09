/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                           www.openmesh.org                                *
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
 *   $Revision: 1 $                                                          *
 *   $Date: 2011-01-09 12:46:45 +0100 (Mo, 09. Jan 2011) $                   *
 *   $LastChangedBy: kremer $                                                *
 *                                                                           *
\*===========================================================================*/

#define POLYHEDRALMESHT_CC

#include <set>

#include "PolyhedralMesh.hh"

namespace OpenVolumeMesh {

//========================================================================================

template <typename VecT>
PolyhedralMesh<VecT>::PolyhedralMesh() :
has_bottom_up_adjacencies_(false),
has_face_normals_(false),
has_vertex_status_(false),
has_edge_status_(false),
has_halfedge_status_(false),
has_face_status_(false),
has_halfface_status_(false),
has_cell_status_(false) {

}

//========================================================================================

/// Build adjacency list in top-down-order
template <typename VecT>
void PolyhedralMesh<VecT>::update_adjacencies() {

    // Clear adjacencies
    outgoing_hes_per_vertex_.clear();
    outgoing_hes_per_vertex_.resize(vertices_.size());

    // Store outgoing halfedges per vertex
    for(unsigned int i = 0; i < edges_.size(); ++i) {

        VertexHandle from = edges_[i].from_vertex();
        if((unsigned int)from >= outgoing_hes_per_vertex_.size()) {
            std::cerr << "update_adjacencies(): Vertex handle is out of bounds!" << std::endl;
            return;
        }
        outgoing_hes_per_vertex_[from].push_back(halfedge_handle(EdgeHandle(i), 0));

        VertexHandle to = edges_[i].to_vertex();
        if((unsigned int)to >= outgoing_hes_per_vertex_.size()) {
            std::cerr << "update_adjacencies(): Vertex handle is out of bounds!" << std::endl;
            return;
        }
        // Store opposite halfedge handle
        outgoing_hes_per_vertex_[to].push_back(halfedge_handle(EdgeHandle(i), 1));
    }

    // Clear
    incident_hfs_per_he_.clear();
    incident_hfs_per_he_.resize(edges_.size() * 2u);

    // Store incident halffaces per halfedge
    for(unsigned int i = 0; i < faces_.size(); ++i) {

        std::vector<HalfEdgeHandle> halfedges = faces_[i].halfedges();

        // Go over all halfedges
        for(typename std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();
                he_it != halfedges.end(); ++he_it) {

            incident_hfs_per_he_[*he_it].push_back(halfface_handle(FaceHandle(i), 0));
            incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].push_back(
                    halfface_handle(FaceHandle(i), 1));
        }
    }

    // Clear
    incident_cell_per_hf_.clear();
    incident_cell_per_hf_.resize(faces_.size() * 2u, InvalidCellHandle);

    for(unsigned int i = 0; i < cells_.size(); ++i) {

        std::vector<HalfFaceHandle> halffaces = cells_[i].halffaces();

        // Go over all halffaces
        for(typename std::vector<HalfFaceHandle>::const_iterator hf_it = halffaces.begin();
                hf_it != halffaces.end(); ++hf_it) {

            if(incident_cell_per_hf_[*hf_it] == InvalidCellHandle) {

                incident_cell_per_hf_[*hf_it] = CellHandle(i);

            } else {
                std::cerr << "Detected non-three-manifold configuration!" << std::endl;
                std::cerr << "Connectivity probably won't work." << std::endl;
                continue;
            }
        }
    }

    /* Put halffaces in clockwise order via the
     * same cell property which now exists.
     * Note, this only works for manifold configurations though.
     * Proceed as follows: Pick one starting halfface. Assuming
     * that all halfface normals point into the incident cell,
     * we find the adjacent halfface within the incident cell
     * along the considered halfedge. We set the found halfface
     * to be the one to be processed next. If we reach an outside
     * region, we try to go back from the starting halfface in reverse
     * order. If the complex is properly connected (the pairwise
     * intersection of two adjacent 3-dimensional cells is always
     * a 2-dimensional entity, namely a facet), such an ordering
     * always exists and will be found. If not, a correct order
     * can not be given and the related iterators will address
     * the related entities arbitrarily.
     */

    for(unsigned int i = 0; i < edges_.size(); ++i) {

    	for(unsigned char s = 0; s <= 1; s++) {

    		HalfEdgeHandle cur_he = halfedge_handle(i, s);
			std::vector<HalfFaceHandle> new_halffaces;
			HalfFaceHandle start_hf = InvalidHalfFaceHandle;
			HalfFaceHandle cur_hf = InvalidHalfFaceHandle;

			// Start with one incident halfface and go
			// into the first direction
			if(incident_hfs_per_he_[cur_he].size() != 0) {

				// Get start halfface
				cur_hf = *incident_hfs_per_he_[cur_he].begin();
				start_hf = cur_hf;

				while(cur_hf != InvalidHalfFaceHandle) {

					// Add halfface
					new_halffaces.push_back(cur_hf);

					// Go to next halfface
					cur_hf = adjacent_halfface_in_cell(cur_hf, cur_he);

					if(cur_hf != InvalidHalfFaceHandle)
						cur_hf = opposite_halfface_handle(cur_hf);

					// End when we're through
					if(cur_hf == start_hf) break;
				}

				// First direction has terminated
				// If new_halffaces has the same size as old (unordered)
				// vector of incident halffaces, we are done here
				// If not, try the other way round
				if(new_halffaces.size() != incident_hfs_per_he_[cur_he].size()) {

					// Get opposite of start halfface
					cur_hf = start_hf;

					 while(cur_hf != InvalidHalfFaceHandle) {

						 cur_hf = opposite_halfface_handle(cur_hf);
						 cur_hf = adjacent_halfface_in_cell(cur_hf, cur_he);

						 if(cur_hf == start_hf) break;

						 if(cur_hf != InvalidHalfFaceHandle)
							 new_halffaces.insert(new_halffaces.begin(), cur_hf);
						 else break;
					}
				}

				// Everything worked just fine, set the new ordered vector
				if(new_halffaces.size() == incident_hfs_per_he_[cur_he].size()) {
					incident_hfs_per_he_[cur_he] = new_halffaces;
				}
			}
    	}
    }


    // Clear
    boundary_faces_.clear();

    // Get boundary faces
    for(unsigned int i = 0; i < faces_.size(); ++i) {

        if(incident_cell_per_hf_[halfface_handle(FaceHandle(i), 0)] == InvalidCellHandle ||
           incident_cell_per_hf_[halfface_handle(FaceHandle(i), 1)] == InvalidCellHandle) {

            // If at least one of two halffaces does not have an
            // incident cell it is a boundary face
            boundary_faces_.push_back(FaceHandle(i));
        }
    }

    has_bottom_up_adjacencies_ = true;
}

//========================================================================================

/// Add vertex
template <typename VecT>
typename PolyhedralMesh<VecT>::VertexHandle
PolyhedralMesh<VecT>::add_vertex(const VecT& _p) {

    // Create vertex object
    OpenVolumeMeshVertex<VecT> v(_p);

    // Store vertex in list
    vertices_.push_back(v);

    // Resize props
    vprops_resize(n_vertices());

    // Resize status
    if(has_vertex_status_) {
        vertex_status_.resize(vertices_.size(), OpenVolumeMeshStatus());
    }

    // Get handle of recently created vertex
    return vertices_.size()-1;
}

//========================================================================================

/// Add edge
template <typename VecT>
typename PolyhedralMesh<VecT>::EdgeHandle
PolyhedralMesh<VecT>::add_edge(const VertexHandle& _fromVertex,
                                 const VertexHandle& _toVertex) {

    if((unsigned int)_fromVertex >= vertices_.size() || (unsigned int)_toVertex >= vertices_.size()) {
        std::cerr << "Vertex handle is out of bounds!" << std::endl;
        return InvalidEdgeHandle;
    }

    // Test if edge does not exist, yet
    for(unsigned int i = 0; i < edges_.size(); ++i) {
        if(edge(EdgeHandle(i)).from_vertex() == _fromVertex && edge(EdgeHandle(i)).to_vertex() == _toVertex) {
            return EdgeHandle(i);
        } else if(edge(EdgeHandle(i)).from_vertex() == _toVertex && edge(EdgeHandle(i)).to_vertex() == _fromVertex) {
            return EdgeHandle(i);
        }
    }

    // Create edge object
    OpenVolumeMeshEdge<VecT> e(_fromVertex, _toVertex);

    // Store edge locally
    edges_.push_back(e);

    // Resize props
    eprops_resize(n_edges());
    heprops_resize(n_halfedges());

    // Resize status
    if(has_edge_status_) {
        edge_status_.resize(edges_.size(), OpenVolumeMeshStatus());
    }
    if(has_halfedge_status_) {
        halfedge_status_.resize(edges_.size() * 2u, OpenVolumeMeshStatus());
    }

    // Get handle of recently created edge
    return edges_.size()-1;
}

//========================================================================================

/// Add face via incident edges
template <typename VecT>
typename PolyhedralMesh<VecT>::FaceHandle
PolyhedralMesh<VecT>::add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck) {

    // Test if all edges are valid
    for(typename std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin();
            it != _halfedges.end(); ++it) {
        if((unsigned int)*it >= edges_.size() * 2u) {
            std::cerr << "Halfedge handle out of bounds!" << std::endl;
            return InvalidFaceHandle;
        }
    }

    // Perform topology check
    if(_topologyCheck) {

        /*
         * Test if halfedges are connected
         *
         * The test works as follows:
         * For every edge in the parameter vector
         * put all incident vertices into a
         * set of either "from"-vertices or "to"-vertices,
         * respectively.
         * If and only if all edges are connected,
         * then both sets are identical.
         */

        std::set<VertexHandle> fromVertices;
        std::set<VertexHandle> toVertices;

        for(typename std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin();
                it != _halfedges.end(); ++it) {

            fromVertices.insert(halfedge(*it).from_vertex());
            toVertices.insert(halfedge(*it).to_vertex());
        }

        // Debug
        assert(fromVertices.size() == toVertices.size());

        for(typename std::set<VertexHandle>::const_iterator v_it = fromVertices.begin();
                v_it != fromVertices.end(); ++v_it) {
            if(toVertices.count(*v_it) != 1) {
                std::cerr << "The specified halfedges are not connected!" << std::endl;
                return InvalidFaceHandle;
            }
        }

        // The halfedges are now guaranteed to be connected
    }

    // Create face
    OpenVolumeMeshFace<VecT> face(_halfedges);

    faces_.push_back(face);

    // Get added face's handle
    FaceHandle fh(faces_.size() - 1);

    // Resize props
    fprops_resize(n_faces());
    hfprops_resize(n_halffaces());

    if(has_face_status_) {
        face_status_.resize(faces_.size(), OpenVolumeMeshStatus());
    }
    if(has_halfface_status_) {
        halfface_status_.resize(faces_.size() * 2u, OpenVolumeMeshStatus());
    }

    // Compute normal if necessary
    if(has_face_normals_) {
        face_normals_.resize(faces_.size(), VecT());
        compute_face_normal(fh);
    }

    // Return handle of recently created face
    return fh;
}

//========================================================================================

/// Add face via incident vertices
/// Define the _vertices in counter-clockwise order (from the "outside")
template <typename VecT>
typename PolyhedralMesh<VecT>::FaceHandle
PolyhedralMesh<VecT>::add_face(const std::vector<VertexHandle>& _vertices) {

    // Test if all vertices exist
    for(typename std::vector<VertexHandle>::const_iterator it = _vertices.begin();
            it != _vertices.end(); ++it) {
        if((unsigned int)*it >= vertices_.size()) {
            std::cerr << "Vertex handle out of bounds!" << std::endl;
            return InvalidFaceHandle;
        }
    }

    // Add edge for each pair of vertices
    std::vector<HalfEdgeHandle> halfedges;
    typename std::vector<VertexHandle>::const_iterator it = _vertices.begin();
    for(; (it+1) != _vertices.end(); ++it) {
        EdgeHandle e_idx = add_edge(*it, *(it+1));

        // Swap halfedge if edge already existed and
        // has been initially defined in reverse orientation
        int swap = 0;
        if(edge(e_idx).to_vertex() == *it) swap = 1;

        halfedges.push_back(halfedge_handle(e_idx, swap));
    }
    EdgeHandle e_idx = add_edge(*it, *_vertices.begin());
    int swap = 0;
    if(edge(e_idx).to_vertex() == *it) swap = 1;
    halfedges.push_back(halfedge_handle(e_idx, swap));

    // Add face
#ifndef NDEBUG
    return add_face(halfedges, true);
#else
    return add_face(halfedges, false);
#endif
}

//========================================================================================

/// Add cell via incident halffaces
template <typename VecT>
typename PolyhedralMesh<VecT>::CellHandle
PolyhedralMesh<VecT>::add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck) {

    // Test if _halffaces contains at least four halffaces
    if(_halffaces.size() < 4u) {
        std::cerr << "A cell must consist of at least four faces!" << std::endl;
        return InvalidCellHandle;
    }

    // Test if halffaces have valid indices
    for(typename std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
            it != _halffaces.end(); ++it) {
        if((unsigned int)*it >= faces_.size() * 2u) {
            std::cerr << "HalfFace handle is out of bounds!" << std::endl;
            return InvalidCellHandle;
        }
    }

    // Perform topology check
    if(_topologyCheck) {

        /*
         * Test if all halffaces are connected => Cell is closed
         *
         * The test works as follows: Put every incident halfedge of
         * all faces into a set. Only if all faces are pairwise
         * connected via one common edge (two opposing halfedges)
         * every pair of halfedges that belong to each commonly shared
         * edge (between each pair of faces) will appear in the set.
         * This means, the number of edges, belonging to each of
         * the halfedges is exactly the number of halfedges divided
         * by 2. This follows from the Euler-Poincare theorem.
         */

        std::set<HalfEdgeHandle> incidentHalfedges;

        for(typename std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
                it != _halffaces.end(); ++it) {

            OpenVolumeMeshFace<VecT> hface = halfface(*it);
            for(typename std::vector<HalfEdgeHandle>::const_iterator he_it = hface.halfedges().begin();
                    he_it != hface.halfedges().end(); ++he_it) {
                incidentHalfedges.insert(*he_it);
            }
        }

        // Now test if all pairs of halfedges were added for all incident edges
        std::set<EdgeHandle> commonEdges;
        for(typename std::set<HalfEdgeHandle>::const_iterator he_it = incidentHalfedges.begin();
                he_it != incidentHalfedges.end(); ++he_it) {
            commonEdges.insert(edge_handle(*he_it));
        }

        if((commonEdges.size() * 2u != incidentHalfedges.size()) ||
                (commonEdges.size() != _halffaces.size() * 2u)) {
            std::cerr << "The specified halffaces are not connected!" << std::endl;
            return InvalidCellHandle;
        }
    }

    // Create new cell
    OpenVolumeMeshCell<VecT> cell(_halffaces);

    cells_.push_back(cell);

    // Resize props
    cprops_resize(n_cells());

    // Resize status
    if(has_cell_status_) {
        cell_status_.resize(cells_.size(), OpenVolumeMeshStatus());
    }

    return cells_.size()-1;
}

//========================================================================================

/// Get vertex with handle _vertexHandle
template <typename VecT>
const OpenVolumeMeshVertex<VecT>&
PolyhedralMesh<VecT>::vertex(const VertexHandle& _vertexHandle) const {

	// Test if vertex is valid
    assert((unsigned int)_vertexHandle < vertices_.size());
    assert(_vertexHandle >= 0);

    return vertices_[_vertexHandle];
}

//========================================================================================

/// Get edge with handle _edgeHandle
template <typename VecT>
const OpenVolumeMeshEdge<VecT>&
PolyhedralMesh<VecT>::edge(const EdgeHandle& _edgeHandle) const {

    // Test if edge is valid
    assert((unsigned int)_edgeHandle < edges_.size());
    assert(_edgeHandle >= 0);

    return edges_[_edgeHandle];
}

//========================================================================================

/// Get face with handle _faceHandle
template <typename VecT>
const OpenVolumeMeshFace<VecT>&
PolyhedralMesh<VecT>::face(const FaceHandle& _faceHandle) const {

    // Test if face is valid
    assert((unsigned int)_faceHandle < faces_.size());
    assert(_faceHandle >= 0);

    return faces_[_faceHandle];
}

//========================================================================================

/// Get cell with handle _cellHandle
template <typename VecT>
const OpenVolumeMeshCell<VecT>&
PolyhedralMesh<VecT>::cell(const CellHandle& _cellHandle) const {

    // Test if cell is valid
    assert((unsigned int)_cellHandle < cells_.size());
    assert(_cellHandle >= 0);

    return cells_[_cellHandle];
}

//========================================================================================

/// Get vertex with handle _vertexHandle
template <typename VecT>
OpenVolumeMeshVertex<VecT>&
PolyhedralMesh<VecT>::vertex(const VertexHandle& _vertexHandle) {

    // Test if vertex is valid
    assert((unsigned int)_vertexHandle < vertices_.size());
    assert(_vertexHandle >= 0);

    return vertices_[_vertexHandle];
}

//========================================================================================

/// Get edge with handle _edgeHandle
template <typename VecT>
OpenVolumeMeshEdge<VecT>&
PolyhedralMesh<VecT>::edge(const EdgeHandle& _edgeHandle) {

    // Test if edge is valid
    assert((unsigned int)_edgeHandle < edges_.size());
    assert(_edgeHandle >= 0);

    return edges_[_edgeHandle];
}

//========================================================================================

/// Get face with handle _faceHandle
template <typename VecT>
OpenVolumeMeshFace<VecT>&
PolyhedralMesh<VecT>::face(const FaceHandle& _faceHandle) {

    // Test if face is valid
    assert((unsigned int)_faceHandle < faces_.size());
    assert(_faceHandle >= 0);

    return faces_[_faceHandle];
}

//========================================================================================

/// Get cell with handle _cellHandle
template <typename VecT>
OpenVolumeMeshCell<VecT>&
PolyhedralMesh<VecT>::cell(const CellHandle& _cellHandle) {

    // Test if cell is valid
    assert((unsigned int)_cellHandle < cells_.size());
    assert(_cellHandle >= 0);

    return cells_[_cellHandle];
}

//========================================================================================

/// Get edge that corresponds to halfedge with handle _halfEdgeHandle
template <typename VecT>
const OpenVolumeMeshEdge<VecT>
PolyhedralMesh<VecT>::halfedge(const HalfEdgeHandle& _halfEdgeHandle) const {

    // Is handle in range?
	assert((unsigned int)_halfEdgeHandle < (edges_.size() * 2));
    assert(_halfEdgeHandle >= 0);

    // In case the handle is even, just return the corresponding edge
    /// Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle % 2 == 0)
        return edges_[(int)(_halfEdgeHandle / 2)];
    else
    	return edges_[(int)(_halfEdgeHandle / 2)].opposite();
}

//========================================================================================

/// Get face that corresponds to halfface with handle _halfFaceHandle
template <typename VecT>
const OpenVolumeMeshFace<VecT>
PolyhedralMesh<VecT>::halfface(const HalfFaceHandle& _halfFaceHandle) const {

    // Is handle in range?
	assert((unsigned int)_halfFaceHandle < (faces_.size() * 2));
    assert(_halfFaceHandle >= 0);

    // In case the handle is not even, just return the corresponding face
    // Otherwise return the opposite halfface via opposite()
    if(_halfFaceHandle % 2 == 0)
        return faces_[(int)(_halfFaceHandle / 2)];
    else
    	return faces_[(int)(_halfFaceHandle / 2)].opposite();
}

//========================================================================================

/// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
template <typename VecT>
const OpenVolumeMeshEdge<VecT>
PolyhedralMesh<VecT>::opposite_halfedge(const HalfEdgeHandle& _halfEdgeHandle) const {

    // Is handle in range?
    assert(_halfEdgeHandle >= 0);
    assert((unsigned int)_halfEdgeHandle < (edges_.size() * 2));

    // In case the handle is not even, just return the corresponding edge
    // Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle % 2 != 0)
        return edges_[(int)(_halfEdgeHandle / 2)];
    else
    	return edges_[(int)(_halfEdgeHandle / 2)].opposite();
}

//========================================================================================

/// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
template <typename VecT>
const OpenVolumeMeshFace<VecT>
PolyhedralMesh<VecT>::opposite_halfface(const HalfFaceHandle& _halfFaceHandle) const {

    // Is handle in range?
    assert(_halfFaceHandle >= 0);
    assert((unsigned int)_halfFaceHandle < (faces_.size() * 2));

    // In case the handle is not even, just return the corresponding face
    // Otherwise return the opposite via the first face's opposite() function
    if(_halfFaceHandle % 2 != 0)
        return faces_[(int)(_halfFaceHandle / 2)];
    else
    	return faces_[(int)(_halfFaceHandle / 2)].opposite();
}

//========================================================================================

template <typename VecT>
const typename PolyhedralMesh<VecT>::HalfEdgeHandle
PolyhedralMesh<VecT>::halfedge(const VertexHandle& _vh1, const VertexHandle& _vh2) const {

    assert(_vh1.idx() < vertices_.size());
    assert(_vh2.idx() < vertices_.size());

    for(VertexOHalfedgeIter voh_it = voh_iter(_vh1); voh_it.valid(); ++voh_it) {
        if(halfedge(*voh_it).to_vertex() == _vh2) {
            return *voh_it;
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

template <typename VecT>
const typename PolyhedralMesh<VecT>::HalfEdgeHandle
PolyhedralMesh<VecT>::next_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const {

    assert((unsigned int)_hfh < faces_.size() * 2u);
    assert((unsigned int)_heh < edges_.size() * 2u);

    std::vector<HalfEdgeHandle> hes = halfface(_hfh).halfedges();

    for(typename std::vector<HalfEdgeHandle>::const_iterator it = hes.begin();
            it != hes.end(); ++it) {
        if(*it == _heh) {
            if((it + 1) != hes.end()) return *(it + 1);
            else return *hes.begin();
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

template <typename VecT>
const typename PolyhedralMesh<VecT>::HalfEdgeHandle
PolyhedralMesh<VecT>::prev_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const {

    assert((unsigned int)_hfh < faces_.size() * 2u);
    assert((unsigned int)_heh < edges_.size() * 2u);

    std::vector<HalfEdgeHandle> hes = halfface(_hfh).halfedges();

    for(typename std::vector<HalfEdgeHandle>::const_iterator it = hes.begin();
            it != hes.end(); ++it) {
        if(*it == _heh) {
            if(it != hes.begin()) return *(it - 1);
            else return *(hes.end() - 1);
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

template <typename VecT>
bool PolyhedralMesh<VecT>::set_vertex(const VertexHandle& _vertexHandle, const VecT& _position) {

    // Check if handle is in range
    if(_vertexHandle < 0 || _vertexHandle >= vertices_.size()) return false;

    // Set position
    vertices_[_vertexHandle].set_position(_position);
    return true;
}

//========================================================================================

template <typename VecT>
bool PolyhedralMesh<VecT>::set_edge(const EdgeHandle& _edgeHandle,
        const VertexHandle& _fromVertex, const VertexHandle& _toVertex) {

    // Check if handle is in range
    if(_edgeHandle < 0 || _edgeHandle >= edges_.size()) return false;

    // Set edge
    edges_[_edgeHandle].set_edge(_fromVertex, _toVertex);
    return true;
}

//========================================================================================

template <typename VecT>
bool PolyhedralMesh<VecT>::set_face(const FaceHandle& _faceHandle, const std::vector<HalfEdgeHandle>& _halfedges) {

    // Check if handle is in range
    if(_faceHandle < 0 || _faceHandle >= faces_.size()) return false;

    // Set face
    faces_[_faceHandle].set_halfedges(_halfedges);
    return true;
}

//========================================================================================

template <typename VecT>
bool PolyhedralMesh<VecT>::set_cell(const CellHandle& _cellHandle, const std::vector<HalfFaceHandle>& _halffaces) {

    // Check if handle is in range
    if(_cellHandle < 0 || _cellHandle >= cells_.size()) return false;

    // Set cell
    cells_[_cellHandle].set_halffaces(_halffaces);
    return true;
}

//========================================================================================

template <typename VecT>
typename PolyhedralMesh<VecT>::HalfFaceHandle
PolyhedralMesh<VecT>::adjacent_halfface_in_cell(const HalfFaceHandle& _halfFaceHandle, const HalfEdgeHandle& _halfEdgeHandle) const {

	if((unsigned int)_halfFaceHandle >= incident_cell_per_hf_.size() || _halfFaceHandle < 0) {
		return InvalidHalfFaceHandle;
	}
	if(incident_cell_per_hf_[_halfFaceHandle] == InvalidCellHandle) {
		// Specified halfface is on the outside of the complex
		return InvalidHalfFaceHandle;
	}

	OpenVolumeMeshCell<VecT> c = cell(incident_cell_per_hf_[_halfFaceHandle]);

	// Make sure that _halfFaceHandle is incident to _halfEdgeHandle
	bool skipped = false;
	bool found = false;
	HalfFaceHandle idx = InvalidHalfFaceHandle;
	for(typename std::vector<HalfFaceHandle>::const_iterator hf_it = c.halffaces().begin();
			hf_it != c.halffaces().end(); ++hf_it) {

		if(*hf_it == _halfFaceHandle) {
			skipped = true;
			continue;
		}

		OpenVolumeMeshFace<VecT> hf = halfface(*hf_it);
		for(typename std::vector<HalfEdgeHandle>::const_iterator he_it = hf.halfedges().begin();
			he_it != hf.halfedges().end(); ++he_it) {

			if(edge_handle(*he_it) == edge_handle(_halfEdgeHandle)) {
				found = true;
				idx = *hf_it;
			}
			if(skipped && found) break;
		}
		if(skipped && found) break;
	}
	return ((skipped && found) ? idx : InvalidHalfFaceHandle);
}

//========================================================================================

template <typename VecT>
typename PolyhedralMesh<VecT>::CellHandle
PolyhedralMesh<VecT>::incident_cell(const HalfFaceHandle& _halfFaceHandle) const {

	if((unsigned int)_halfFaceHandle >= incident_cell_per_hf_.size() || _halfFaceHandle < 0) {
		return InvalidCellHandle;
	}

	return incident_cell_per_hf_[_halfFaceHandle];
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::garbage_collection() {

    if(!has_vertex_status() || !has_edge_status() || !has_face_status() || !has_cell_status()) {
        std::cerr << "garbage_collection() requires all " <<
                "entity types to have status properties!" << std::endl;
        return;
    }

    /*
     * Perform these steps from vertex to cell:
     * ===========================================
     * 1. Delete entity with handle h from vector
     * 2. Delete all higher dimensional entities
     *    containing a handle to the delete entity (set delete flag)
     * 3. Replace all handle indices h_i > h with (h_i - 1) in all
     *    higher dimensional entities (steps 2 and 3 can be combined)
     * 4. Delete property and status data of h
     * ===========================================
     * 5. At last, call update_adjacencies()
     */

    VertexHandle vh(0);
    for(typename Vertices::iterator v_it = vertices_.begin(); v_it != vertices_.end();) {

        if(!status(vh).deleted()) {
            ++v_it;
            vh.idx(vh.idx() + 1);
            continue;
        }

        // Step 1
        v_it = vertices_.erase(v_it);

        // Step 2 + 3
        for(EdgeIter e_it = e_iter(); e_it.valid(); ++e_it) {

            if(status(*e_it).deleted()) continue;

            // Delete edges that are incident to vh
            Edge& e = edge(*e_it);
            if(e.from_vertex() == vh || e.to_vertex() == vh) {
                status(*e_it).set_deleted(true);
                continue;
            }

            // Replace handles with higher index
            if(e.from_vertex().idx() > vh.idx()) {
                e.set_from_vertex(VertexHandle(e.from_vertex().idx() - 1));
            }

            if(e.to_vertex().idx() > vh.idx()) {
                e.set_to_vertex(VertexHandle(e.to_vertex().idx() - 1));
            }
        }

        // Step 4
        vertex_status_.erase(vertex_status_.begin() + vh.idx());
        remove_vprop_element((size_t)vh.idx());
    }

    EdgeHandle eh(0);
    for(typename Edges::iterator e_it = edges_.begin(); e_it != edges_.end();) {

        if(!status(eh).deleted()) {
            ++e_it;
            eh.idx(eh.idx() + 1);
            continue;
        }

        // Step 1
        e_it = edges_.erase(e_it);

        // Step 2 + 3
        for(FaceIter f_it = f_iter(); f_it.valid(); ++f_it) {

            if(status(*f_it).deleted()) continue;

            // Delete faces that are incident to eh
            Face& f = face(*f_it);
            std::vector<HalfEdgeHandle> hes = face(*f_it).halfedges();
            for(typename std::vector<HalfEdgeHandle>::iterator he_it = hes.begin();
                    he_it != hes.end(); ++he_it) {

                EdgeHandle t_eh = edge_handle(*he_it);

                if(t_eh == eh) {
                    status(*f_it).set_deleted(true);
                    continue;
                }

                // Replace handles with higher index
                if(t_eh > eh) {
                    bool orig = (*he_it == halfedge_handle(t_eh, 0));
                    *he_it = halfedge_handle(t_eh - 1, (orig ? 0 : 1));
                }
            }
            f.set_halfedges(hes);
            f.compute_opposite_halfedges();
        }

        // Step 4
        edge_status_.erase(edge_status_.begin() + eh.idx());
        remove_eprop_element((size_t)eh.idx());
        remove_heprop_element((size_t)halfedge_handle(eh, 0));
        remove_heprop_element((size_t)halfedge_handle(eh, 1));
    }

    FaceHandle fh(0);
    for(typename Faces::iterator f_it = faces_.begin(); f_it != faces_.end();) {

        if(!status(fh).deleted()) {
            ++f_it;
            fh.idx(fh.idx() + 1);
            continue;
        }

        // Step 1
        f_it = faces_.erase(f_it);

        // Step 2 + 3
        for(CellIter c_it = c_iter(); c_it.valid(); ++c_it) {

            if(status(*c_it).deleted()) continue;

            // Delete cells that are incident to fh
            Cell& c = cell(*c_it);
            std::vector<HalfFaceHandle> hfs = cell(*c_it).halffaces();
            for(typename std::vector<HalfFaceHandle>::iterator hf_it = hfs.begin();
                    hf_it != hfs.end(); ++hf_it) {

                FaceHandle t_fh = face_handle(*hf_it);

                if(t_fh == fh) {
                    status(*c_it).set_deleted(true);
                    continue;
                }

                // Replace handles with higher index
                if(t_fh > fh) {
                    bool orig = (*hf_it == halfface_handle(t_fh, 0));
                    *hf_it = halfface_handle(t_fh - 1, (orig ? 0 : 1));
                }
            }
            c.set_halffaces(hfs);
        }

        // Step 4
        face_status_.erase(face_status_.begin() + fh.idx());
        remove_fprop_element((size_t)fh.idx());
        remove_hfprop_element((size_t)halfface_handle(fh, 0));
        remove_hfprop_element((size_t)halfface_handle(fh, 1));
    }

    CellHandle ch(0);
    for(typename Cells::iterator c_it = cells_.begin(); c_it != cells_.end();) {

        if(!status(ch).deleted()) {
            ++c_it;
            ch.idx(ch.idx() + 1);
            continue;
        }

        // Step 1
        c_it = cells_.erase(c_it);

        // No step 2 and 3 needed

        // Step 4
        cell_status_.erase(cell_status_.begin() + ch.idx());
        remove_cprop_element((size_t)ch.idx());
    }

    // Step 5
    update_adjacencies();
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::request_face_normals() {

    // Resize face normals vector
    face_normals_.resize(faces_.size(), VecT());

    // Set flag
    has_face_normals_ = true;

    // Fill normals
    update_face_normals();
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::update_face_normals() {

    assert(face_normals_.size() == faces_.size());

    for(FaceIter f_it = f_iter(); f_it.valid(); ++f_it) {
        // Assume the face is planar, so just take the
        // first two edges
        compute_face_normal(*f_it);
    }
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::compute_face_normal(const FaceHandle& _fh) {

    assert(face_normals_.size() == faces_.size());

    if(face(_fh).halfedges().size() < 2) {
        std::cerr << "Warning: Degenerate face detected!" << std::endl;
        return;
    }

    // Always try to compute the outside normals
    HalfFaceHandle hfh = is_boundary(halfface_handle(_fh, 0)) ?
            halfface_handle(_fh, 0) : halfface_handle(_fh, 1);

    std::vector<HalfEdgeHandle> halfedges = halfface(hfh).halfedges();
    typename std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();

    VecT p1 = vertex(halfedge(*he_it).from_vertex()).position();
    VecT p2 = vertex(halfedge(*he_it).to_vertex()).position();
    ++he_it;
    VecT p3 = vertex(halfedge(*he_it).to_vertex()).position();

    VecT n = (p3 - p2) % (p1 - p2);
    n.normalize();

    face_normals_[_fh.idx()] = n;
}

//========================================================================================

template <typename VecT>
const VecT& PolyhedralMesh<VecT>::face_normal(const FaceHandle& _fh) const {

    assert((unsigned int)_fh < face_normals_.size());

    return face_normals_[_fh];
}

//========================================================================================

template <typename VecT>
const VecT PolyhedralMesh<VecT>::halfface_normal(const HalfFaceHandle& _hfh) const {

    assert((unsigned int)_hfh < face_normals_.size() * 2u);

    if((unsigned int)_hfh % 2u == 1)
        return face_normals_[face_handle(_hfh)] * (-1);

    return face_normals_[face_handle(_hfh)];
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::request_vertex_status() {

    // Resize vector
    vertex_status_.resize(vertices_.size(), OpenVolumeMeshStatus());

    has_vertex_status_ = true;
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::request_edge_status() {

    // Resize vector
    edge_status_.resize(edges_.size(), OpenVolumeMeshStatus());

    has_edge_status_ = true;
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::request_halfedge_status() {

    // Resize vector
    halfedge_status_.resize(edges_.size() * 2u, OpenVolumeMeshStatus());

    has_halfedge_status_ = true;
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::request_face_status() {

    // Resize vector
    face_status_.resize(faces_.size(), OpenVolumeMeshStatus());

    has_face_status_ = true;
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::request_halfface_status() {

    // Resize vector
    halfface_status_.resize(faces_.size() * 2u, OpenVolumeMeshStatus());

    has_halfface_status_ = true;
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::request_cell_status() {

    // Resize vector
    cell_status_.resize(cells_.size(), OpenVolumeMeshStatus());

    has_cell_status_ = true;
}

//========================================================================================

template <typename VecT>
void PolyhedralMesh<VecT>::request_status() {

    // Request all status types
    request_vertex_status();
    request_edge_status();
    request_halfedge_status();
    request_face_status();
    request_halfface_status();
    request_cell_status();
}

//========================================================================================

template <typename VecT>
OpenVolumeMeshStatus& PolyhedralMesh<VecT>::status(const VertexHandle& _vh) {

    assert((unsigned int)_vh < vertex_status_.size());

    return vertex_status_[_vh];
}

//========================================================================================

template <typename VecT>
OpenVolumeMeshStatus& PolyhedralMesh<VecT>::status(const EdgeHandle& _eh) {

    assert((unsigned int)_eh < edge_status_.size());

    return edge_status_[_eh];
}

//========================================================================================

template <typename VecT>
OpenVolumeMeshStatus& PolyhedralMesh<VecT>::status(const HalfEdgeHandle& _heh) {

    assert((unsigned int)_heh < halfedge_status_.size());

    return halfedge_status_[_heh];
}

//========================================================================================

template <typename VecT>
OpenVolumeMeshStatus& PolyhedralMesh<VecT>::status(const FaceHandle& _fh) {

    assert((unsigned int)_fh < face_status_.size());

    return face_status_[_fh];
}

//========================================================================================

template <typename VecT>
OpenVolumeMeshStatus& PolyhedralMesh<VecT>::status(const HalfFaceHandle& _hfh) {

    assert((unsigned int)_hfh < halfface_status_.size());

    return halfface_status_[_hfh];
}

//========================================================================================

template <typename VecT>
OpenVolumeMeshStatus& PolyhedralMesh<VecT>::status(const CellHandle& _ch) {

    assert((unsigned int)_ch < cell_status_.size());

    return cell_status_[_ch];
}

//========================================================================================

} // Namespace OpenVolumeMesh
