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

#include "TopologyKernel.hh"

namespace OpenVolumeMesh {

// Initialize constants
const VertexHandle      TopologyKernel::InvalidVertexHandle   = VertexHandle(-1);
const EdgeHandle        TopologyKernel::InvalidEdgeHandle     = EdgeHandle(-1);
const HalfEdgeHandle    TopologyKernel::InvalidHalfEdgeHandle = HalfEdgeHandle(-1);
const FaceHandle        TopologyKernel::InvalidFaceHandle     = FaceHandle(-1);
const HalfFaceHandle    TopologyKernel::InvalidHalfFaceHandle = HalfFaceHandle(-1);
const CellHandle        TopologyKernel::InvalidCellHandle     = CellHandle(-1);

TopologyKernel::TopologyKernel() :
    n_vertices_(0u),
    has_vertex_adjacencies_(false),
    has_edge_adjacencies_(false),
    has_face_adjacencies_(false) {

}

TopologyKernel::~TopologyKernel() {
}

//========================================================================================

/// Add edge
EdgeHandle TopologyKernel::add_edge(const VertexHandle& _fromVertex,
                                    const VertexHandle& _toVertex) {

    if((unsigned int)_fromVertex >= n_vertices() || (unsigned int)_toVertex >= n_vertices()) {
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
    OpenVolumeMeshEdge e(_fromVertex, _toVertex);

    // Store edge locally
    edges_.push_back(e);

    // Resize props
    resize_eprops(n_edges());

    // Get handle of recently created edge
    return EdgeHandle((int)edges_.size()-1);
}

//========================================================================================

/// Add face via incident edges
FaceHandle TopologyKernel::add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck) {

    // Test if all edges are valid
    for(std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin();
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

        for(std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin();
                it != _halfedges.end(); ++it) {

            fromVertices.insert(halfedge(*it).from_vertex());
            toVertices.insert(halfedge(*it).to_vertex());
        }

        for(std::set<VertexHandle>::const_iterator v_it = fromVertices.begin();
                v_it != fromVertices.end(); ++v_it) {
            if(toVertices.count(*v_it) != 1) {
                std::cerr << "The specified halfedges are not connected!" << std::endl;
                return InvalidFaceHandle;
            }
        }

        // The halfedges are now guaranteed to be connected
    }

    // Create face
    OpenVolumeMeshFace face(_halfedges);

    faces_.push_back(face);

    // Get added face's handle
    FaceHandle fh(faces_.size() - 1);

    // Resize props
    resize_fprops(n_faces());

    // Return handle of recently created face
    return fh;
}

//========================================================================================

/// Add face via incident vertices
/// Define the _vertices in counter-clockwise order (from the "outside")
FaceHandle TopologyKernel::add_face(const std::vector<VertexHandle>& _vertices) {

    // Test if all vertices exist
    for(std::vector<VertexHandle>::const_iterator it = _vertices.begin();
            it != _vertices.end(); ++it) {
        if((unsigned int)*it >= n_vertices()) {
            std::cerr << "Vertex handle out of bounds!" << std::endl;
            return InvalidFaceHandle;
        }
    }

    // Add edge for each pair of vertices
    std::vector<HalfEdgeHandle> halfedges;
    std::vector<VertexHandle>::const_iterator it = _vertices.begin();
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
CellHandle TopologyKernel::add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck) {

    // Test if _halffaces contains at least four halffaces
    if(_halffaces.size() < 4u) {
        std::cerr << "A cell must consist of at least four faces!" << std::endl;
        return InvalidCellHandle;
    }

    // Test if halffaces have valid indices
    for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
            it != _halffaces.end(); ++it) {
        if((unsigned int)*it >= faces_.size() * 2u) {
            std::cerr << "HalfFace handle is out of bounds!" << std::endl;
            return InvalidCellHandle;
        }
    }

    // Perform topology check
    if(_topologyCheck) {

        /*
         * Test if all halffaces are connected and form a two-manifold
         * => Cell is closed
         *
         * This test is simple: The number of involved half-edges has to be
         * exactly twice the number of involved edges.
         */

        std::set<HalfEdgeHandle> incidentHalfedges;
        std::set<EdgeHandle>     incidentEdges;

        for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
                it != _halffaces.end(); ++it) {

            OpenVolumeMeshFace hface = halfface(*it);
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = hface.halfedges().begin();
                    he_it != hface.halfedges().end(); ++he_it) {
                incidentHalfedges.insert(*he_it);
                incidentEdges.insert(edge_handle(*he_it));
            }
        }

        if(incidentHalfedges.size() != (incidentEdges.size() * 2u)) {
            std::cerr << "The specified halffaces are not connected!" << std::endl;
            return InvalidCellHandle;
        }

        // The halffaces are now guaranteed to form a two-manifold
    }

    // Create new cell
    OpenVolumeMeshCell cell(_halffaces);

    cells_.push_back(cell);

    // Resize props
    resize_cprops(n_cells());

    return CellHandle((int)cells_.size()-1);
}

//========================================================================================

/// Get edge with handle _edgeHandle
const OpenVolumeMeshEdge& TopologyKernel::edge(const EdgeHandle& _edgeHandle) const {

    // Test if edge is valid
    assert((unsigned int)_edgeHandle < edges_.size());
    assert(_edgeHandle >= 0);

    return edges_[_edgeHandle];
}

//========================================================================================

/// Get face with handle _faceHandle
const OpenVolumeMeshFace& TopologyKernel::face(const FaceHandle& _faceHandle) const {

    // Test if face is valid
    assert((unsigned int)_faceHandle < faces_.size());
    assert(_faceHandle >= 0);

    return faces_[_faceHandle];
}

//========================================================================================

/// Get cell with handle _cellHandle
const OpenVolumeMeshCell& TopologyKernel::cell(const CellHandle& _cellHandle) const {

    // Test if cell is valid
    assert((unsigned int)_cellHandle < cells_.size());
    assert(_cellHandle >= 0);

    return cells_[_cellHandle];
}

//========================================================================================

/// Get edge with handle _edgeHandle
OpenVolumeMeshEdge& TopologyKernel::edge(const EdgeHandle& _edgeHandle) {

    // Test if edge is valid
    assert((unsigned int)_edgeHandle < edges_.size());
    assert(_edgeHandle >= 0);

    return edges_[_edgeHandle];
}

//========================================================================================

/// Get face with handle _faceHandle
OpenVolumeMeshFace& TopologyKernel::face(const FaceHandle& _faceHandle) {

    // Test if face is valid
    assert((unsigned int)_faceHandle < faces_.size());
    assert(_faceHandle >= 0);

    return faces_[_faceHandle];
}

//========================================================================================

/// Get cell with handle _cellHandle
OpenVolumeMeshCell& TopologyKernel::cell(const CellHandle& _cellHandle) {

    // Test if cell is valid
    assert((unsigned int)_cellHandle < cells_.size());
    assert(_cellHandle >= 0);

    return cells_[_cellHandle];
}

//========================================================================================

/// Get edge that corresponds to halfedge with handle _halfEdgeHandle
const OpenVolumeMeshEdge TopologyKernel::halfedge(const HalfEdgeHandle& _halfEdgeHandle) const {

    // Is handle in range?
    assert((unsigned int)_halfEdgeHandle < (edges_.size() * 2));
    assert(_halfEdgeHandle >= 0);

    // In case the handle is even, just return the corresponding edge
    /// Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle % 2 == 0)
        return edges_[(int)(_halfEdgeHandle / 2)];
    else
        return opposite_halfedge(edges_[(int)(_halfEdgeHandle / 2)]);
}

//========================================================================================

/// Get face that corresponds to halfface with handle _halfFaceHandle
const OpenVolumeMeshFace TopologyKernel::halfface(const HalfFaceHandle& _halfFaceHandle) const {

    // Is handle in range?
    assert((unsigned int)_halfFaceHandle < (faces_.size() * 2));
    assert(_halfFaceHandle >= 0);

    // In case the handle is not even, just return the corresponding face
    // Otherwise return the opposite halfface via opposite()
    if(_halfFaceHandle % 2 == 0)
        return faces_[(int)(_halfFaceHandle / 2)];
    else
        return opposite_halfface(faces_[(int)(_halfFaceHandle / 2)]);
}

//========================================================================================

/// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
const OpenVolumeMeshEdge TopologyKernel::opposite_halfedge(const HalfEdgeHandle& _halfEdgeHandle) const {

    // Is handle in range?
    assert(_halfEdgeHandle >= 0);
    assert((unsigned int)_halfEdgeHandle < (edges_.size() * 2));

    // In case the handle is not even, just return the corresponding edge
    // Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle % 2 != 0)
        return edges_[(int)(_halfEdgeHandle / 2)];
    else
        return opposite_halfedge(edges_[(int)(_halfEdgeHandle / 2)]);
}

//========================================================================================

/// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
const OpenVolumeMeshFace TopologyKernel::opposite_halfface(const HalfFaceHandle& _halfFaceHandle) const {

    // Is handle in range?
    assert(_halfFaceHandle >= 0);
    assert((unsigned int)_halfFaceHandle < (faces_.size() * 2));

    // In case the handle is not even, just return the corresponding face
    // Otherwise return the opposite via the first face's opposite() function
    if(_halfFaceHandle % 2 != 0)
        return faces_[(int)(_halfFaceHandle / 2)];
    else
        return opposite_halfface(faces_[(int)(_halfFaceHandle / 2)]);
}

//========================================================================================

const HalfEdgeHandle TopologyKernel::halfedge(const VertexHandle& _vh1, const VertexHandle& _vh2) const {

    assert(_vh1.idx() < (int)n_vertices());
    assert(_vh2.idx() < (int)n_vertices());

    for(VertexOHalfEdgeIter voh_it = voh_iter(_vh1); voh_it.valid(); ++voh_it) {
        if(halfedge(*voh_it).to_vertex() == _vh2) {
            return *voh_it;
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

const HalfEdgeHandle TopologyKernel::next_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const {

    assert((unsigned int)_hfh < faces_.size() * 2u);
    assert((unsigned int)_heh < edges_.size() * 2u);

    std::vector<HalfEdgeHandle> hes = halfface(_hfh).halfedges();

    for(std::vector<HalfEdgeHandle>::const_iterator it = hes.begin();
            it != hes.end(); ++it) {
        if(*it == _heh) {
            if((it + 1) != hes.end()) return *(it + 1);
            else return *hes.begin();
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

const HalfEdgeHandle TopologyKernel::prev_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const {

    assert((unsigned int)_hfh < faces_.size() * 2u);
    assert((unsigned int)_heh < edges_.size() * 2u);

    std::vector<HalfEdgeHandle> hes = halfface(_hfh).halfedges();

    for(std::vector<HalfEdgeHandle>::const_iterator it = hes.begin();
            it != hes.end(); ++it) {
        if(*it == _heh) {
            if(it != hes.begin()) return *(it - 1);
            else return *(hes.end() - 1);
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

void TopologyKernel::update_adjacencies() {

    update_vertex_adjacencies();

    update_edge_adjacencies();

    update_face_adjacencies();
}

//========================================================================================

void TopologyKernel::update_vertex_adjacencies() {

    // Clear adjacencies
    outgoing_hes_per_vertex_.clear();
    outgoing_hes_per_vertex_.resize(n_vertices());

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

    has_vertex_adjacencies_ = true;
}

//========================================================================================

void TopologyKernel::update_edge_adjacencies() {

    // Clear
    incident_hfs_per_he_.clear();
    incident_hfs_per_he_.resize(edges_.size() * 2u);

    // Store incident halffaces per halfedge
    for(unsigned int i = 0; i < faces_.size(); ++i) {

        std::vector<HalfEdgeHandle> halfedges = faces_[i].halfedges();

        // Go over all halfedges
        for(std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();
                he_it != halfedges.end(); ++he_it) {

            incident_hfs_per_he_[*he_it].push_back(halfface_handle(FaceHandle(i), 0));
            incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].push_back(
                    halfface_handle(FaceHandle(i), 1));
        }
    }

    has_edge_adjacencies_ = true;
}

//========================================================================================

void TopologyKernel::update_face_adjacencies() {

    // Clear
    incident_cell_per_hf_.clear();
    incident_cell_per_hf_.resize(faces_.size() * 2u, InvalidCellHandle);

    for(unsigned int i = 0; i < cells_.size(); ++i) {

        std::vector<HalfFaceHandle> halffaces = cells_[i].halffaces();

        // Go over all halffaces
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = halffaces.begin();
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
     * can not be given and, as a result, the related iterators
     * will address the related entities in an arbitrary fashion.
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

    // Compute boundary faces

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

    has_face_adjacencies_ = true;
}

//========================================================================================

HalfFaceHandle
TopologyKernel::adjacent_halfface_in_cell(const HalfFaceHandle& _halfFaceHandle, const HalfEdgeHandle& _halfEdgeHandle) const {

    if((unsigned int)_halfFaceHandle >= incident_cell_per_hf_.size() || _halfFaceHandle < 0) {
        return InvalidHalfFaceHandle;
    }
    if(incident_cell_per_hf_[_halfFaceHandle] == InvalidCellHandle) {
        // Specified halfface is on the outside of the complex
        return InvalidHalfFaceHandle;
    }

    OpenVolumeMeshCell c = cell(incident_cell_per_hf_[_halfFaceHandle]);

    // Make sure that _halfFaceHandle is incident to _halfEdgeHandle
    bool skipped = false;
    bool found = false;
    HalfFaceHandle idx = InvalidHalfFaceHandle;
    for(std::vector<HalfFaceHandle>::const_iterator hf_it = c.halffaces().begin();
            hf_it != c.halffaces().end(); ++hf_it) {

        if(*hf_it == _halfFaceHandle) {
            skipped = true;
            continue;
        }

        OpenVolumeMeshFace hf = halfface(*hf_it);
        for(std::vector<HalfEdgeHandle>::const_iterator he_it = hf.halfedges().begin();
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

CellHandle TopologyKernel::incident_cell(const HalfFaceHandle& _halfFaceHandle) const {

    if((unsigned int)_halfFaceHandle >= incident_cell_per_hf_.size() || _halfFaceHandle < 0) {
        return InvalidCellHandle;
    }

    return incident_cell_per_hf_[_halfFaceHandle];
}

} // Namespace OpenVolumeMesh
