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

#include <OpenVolumeMesh/System/FunctionalInclude.hh>
#include <queue>

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
    v_bottom_up_(true),
    e_bottom_up_(true),
    f_bottom_up_(true) {

}

TopologyKernel::~TopologyKernel() {
}

//========================================================================================

VertexHandle TopologyKernel::add_vertex() {

    ++n_vertices_;

    // Create item for vertex bottom-up incidences
    if(v_bottom_up_) {
        outgoing_hes_per_vertex_.resize(n_vertices_);
    }

    // Resize vertex props
    resize_vprops(n_vertices_);

    // Return 0-indexed handle
    return VertexHandle((int)(n_vertices_ - 1));
}

//========================================================================================

/// Add edge
EdgeHandle TopologyKernel::add_edge(const VertexHandle& _fromVertex,
                                    const VertexHandle& _toVertex,
                                    bool _allowDuplicates) {

#ifndef NDEBUG
    if((unsigned int)_fromVertex.idx() >= n_vertices() || (unsigned int)_toVertex.idx() >= n_vertices()) {
        std::cerr << "Vertex handle is out of bounds!" << std::endl;
        return InvalidEdgeHandle;
    }
#endif

    // Test if edge does not exist, yet
    if(!_allowDuplicates) {
        if(v_bottom_up_) {

            assert(outgoing_hes_per_vertex_.size() > (unsigned int)_fromVertex.idx());
            std::vector<HalfEdgeHandle>& ohes = outgoing_hes_per_vertex_[_fromVertex.idx()];
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = ohes.begin(),
                    he_end = ohes.end(); he_it != he_end; ++he_it) {
                if(halfedge(*he_it).to_vertex() == _toVertex) {
                    return edge_handle(*he_it);
                }
            }
        } else {
            for(unsigned int i = 0; i < edges_.size(); ++i) {
                if(edge(EdgeHandle(i)).from_vertex() == _fromVertex && edge(EdgeHandle(i)).to_vertex() == _toVertex) {
                    return EdgeHandle(i);
                } else if(edge(EdgeHandle(i)).from_vertex() == _toVertex && edge(EdgeHandle(i)).to_vertex() == _fromVertex) {
                    return EdgeHandle(i);
                }
            }
        }
    }

    // Create edge object
    OpenVolumeMeshEdge e(_fromVertex, _toVertex);

    // Store edge locally
    edges_.push_back(e);

    // Resize props
    resize_eprops(n_edges());

    EdgeHandle eh((int)edges_.size()-1);

    // Update vertex bottom-up incidences
    if(v_bottom_up_) {
        assert(outgoing_hes_per_vertex_.size() > (unsigned int)_fromVertex.idx());
        assert(outgoing_hes_per_vertex_.size() > (unsigned int)_toVertex.idx());
        outgoing_hes_per_vertex_[_fromVertex.idx()].push_back(halfedge_handle(eh, 0));
        outgoing_hes_per_vertex_[_toVertex.idx()].push_back(halfedge_handle(eh, 1));
    }

    // Create item for edge bottom-up incidences
    if(e_bottom_up_) {
        incident_hfs_per_he_.resize(n_halfedges());
    }

    // Get handle of recently created edge
    return eh;
}

//========================================================================================

/// Add face via incident edges
FaceHandle TopologyKernel::add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck) {

#ifndef NDEBUG
    // Test if all edges are valid
    for(std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin(),
            end = _halfedges.end(); it != end; ++it) {
        if((unsigned int)it->idx() >= edges_.size() * 2u) {
            std::cerr << "Halfedge handle out of bounds!" << std::endl;
            return InvalidFaceHandle;
        }
    }
#endif

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

        for(std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin(),
            end = _halfedges.end(); it != end; ++it) {

            fromVertices.insert(halfedge(*it).from_vertex());
            toVertices.insert(halfedge(*it).to_vertex());
        }

        for(std::set<VertexHandle>::const_iterator v_it = fromVertices.begin(),
                v_end = fromVertices.end(); v_it != v_end; ++v_it) {
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

    // Update edge bottom-up incidences
    if(e_bottom_up_) {

        for(std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin(),
            end = _halfedges.end(); it != end; ++it) {
            assert(incident_hfs_per_he_.size() > (unsigned int)it->idx());
            assert(incident_hfs_per_he_.size() > (unsigned int)opposite_halfedge_handle(*it).idx());
            incident_hfs_per_he_[it->idx()].push_back(halfface_handle(fh, 0));
            incident_hfs_per_he_[opposite_halfedge_handle(*it).idx()].push_back(halfface_handle(fh, 1));
        }
    }

    // Create item for face bottom-up incidences
    if(f_bottom_up_) {
        incident_cell_per_hf_.resize(n_halffaces(), InvalidCellHandle);
    }

    // Return handle of recently created face
    return fh;
}

//========================================================================================

/// Add face via incident vertices
/// Define the _vertices in counter-clockwise order (from the "outside")
FaceHandle TopologyKernel::add_face(const std::vector<VertexHandle>& _vertices) {

#ifndef NDEBUG
    // Test if all vertices exist
    for(std::vector<VertexHandle>::const_iterator it = _vertices.begin(),
            end = _vertices.end(); it != end; ++it) {
        if((unsigned int)it->idx() >= n_vertices()) {
            std::cerr << "Vertex handle out of bounds!" << std::endl;
            return InvalidFaceHandle;
        }
    }
#endif

    // Add edge for each pair of vertices
    std::vector<HalfEdgeHandle> halfedges;
    std::vector<VertexHandle>::const_iterator it = _vertices.begin();
    std::vector<VertexHandle>::const_iterator end = _vertices.end();
    for(; (it+1) != end; ++it) {
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

void TopologyKernel::reorder_incident_halffaces(const EdgeHandle& _eh) {

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

    for(unsigned char s = 0; s <= 1; s++) {

        HalfEdgeHandle cur_he = halfedge_handle(_eh, s);
        std::vector<HalfFaceHandle> new_halffaces;
        HalfFaceHandle start_hf = InvalidHalfFaceHandle;
        HalfFaceHandle cur_hf = InvalidHalfFaceHandle;

        // Start with one incident halfface and go
        // into the first direction
        assert(incident_hfs_per_he_.size() > (unsigned int)cur_he.idx());

        if(incident_hfs_per_he_[cur_he.idx()].size() != 0) {

            // Get start halfface
            cur_hf = *incident_hfs_per_he_[cur_he.idx()].begin();
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
            if(new_halffaces.size() != incident_hfs_per_he_[cur_he.idx()].size()) {

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
            if(new_halffaces.size() == incident_hfs_per_he_[cur_he.idx()].size()) {
                incident_hfs_per_he_[cur_he.idx()] = new_halffaces;
            }
        }
    }
}

//========================================================================================

/// Add cell via incident halffaces
CellHandle TopologyKernel::add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck) {

#ifndef NDEBUG
    // Test if halffaces have valid indices
    for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin(),
            end = _halffaces.end(); it != end; ++it) {
        if((unsigned int)it->idx() >= faces_.size() * 2u) {
            std::cerr << "HalfFace handle is out of bounds!" << std::endl;
            return InvalidCellHandle;
        }
    }
#endif

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

        for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin(),
                end = _halffaces.end(); it != end; ++it) {

            OpenVolumeMeshFace hface = halfface(*it);
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = hface.halfedges().begin(),
                    he_end = hface.halfedges().end(); he_it != he_end; ++he_it) {
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

    CellHandle ch((int)cells_.size()-1);

    // Update face bottom-up incidences
    if(f_bottom_up_) {

        std::set<EdgeHandle> cell_edges;
        for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin(),
                end = _halffaces.end(); it != end; ++it) {
            assert(incident_cell_per_hf_.size() > (unsigned int)it->idx());

            if(_topologyCheck) {
                if(incident_cell_per_hf_[it->idx()] != InvalidCellHandle) {
                    std::cerr << "Warning: One of the specified half-faces is already incident to another cell!" << std::endl;
                }
            }

            // Overwrite incident cell for current half-face
            incident_cell_per_hf_[it->idx()] = ch;

            // Collect all edges of cell
            const std::vector<HalfEdgeHandle> hes = halfface(*it).halfedges();
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                    he_end = hes.end(); he_it != he_end; ++he_it) {
                cell_edges.insert(edge_handle(*he_it));
            }
        }

        if(e_bottom_up_) {

            // Try to reorder all half-faces w.r.t.
            // their incident half-edges such that all
            // half-faces are in cyclic order around
            // a half-edge
            for(std::set<EdgeHandle>::const_iterator e_it = cell_edges.begin(),
                    e_end = cell_edges.end(); e_it != e_end; ++e_it) {
                reorder_incident_halffaces(*e_it);
            }
        }
    }

    return ch;
}

//========================================================================================

/// Set the vertices of an edge
void TopologyKernel::set_edge(const EdgeHandle& _eh, const VertexHandle& _fromVertex, const VertexHandle& _toVertex) {

    Edge& e = edge(_eh);

    // Update bottom-up entries
    if(has_vertex_bottom_up_incidences()) {

        const VertexHandle& fv = e.from_vertex();
        const VertexHandle& tv = e.to_vertex();

        const HalfEdgeHandle heh0 = halfedge_handle(_eh, 0);
        const HalfEdgeHandle heh1 = halfedge_handle(_eh, 1);

        std::remove(outgoing_hes_per_vertex_[fv.idx()].begin(), outgoing_hes_per_vertex_[fv.idx()].end(), heh0);
        std::remove(outgoing_hes_per_vertex_[tv.idx()].begin(), outgoing_hes_per_vertex_[tv.idx()].end(), heh1);

        outgoing_hes_per_vertex_[_fromVertex.idx()].push_back(heh0);
        outgoing_hes_per_vertex_[_toVertex.idx()].push_back(heh1);
    }

    e.set_from_vertex(_fromVertex);
    e.set_to_vertex(_toVertex);
}

//========================================================================================

/// Set the half-edges of a face
void TopologyKernel::set_face(const FaceHandle& _fh, const std::vector<HalfEdgeHandle>& _hes) {

    Face& f = face(_fh);

    if(has_edge_bottom_up_incidences()) {

        const HalfFaceHandle hf0 = halfface_handle(_fh, 0);
        const HalfFaceHandle hf1 = halfface_handle(_fh, 1);

        const std::vector<HalfEdgeHandle>& hes = f.halfedges();

        for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                he_end = hes.end(); he_it != he_end; ++he_it) {

            std::remove(incident_hfs_per_he_[he_it->idx()].begin(),
                        incident_hfs_per_he_[he_it->idx()].end(), hf0);
            std::remove(incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].begin(),
                        incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].end(), hf1);
        }

        for(std::vector<HalfEdgeHandle>::const_iterator he_it = _hes.begin(),
                he_end = _hes.end(); he_it != he_end; ++he_it) {

            incident_hfs_per_he_[he_it->idx()].push_back(hf0);
            incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].push_back(hf1);
        }

        // TODO: Reorder incident half-faces
    }

    f.set_halfedges(_hes);
}

//========================================================================================

/// Set the half-faces of a cell
void TopologyKernel::set_cell(const CellHandle& _ch, const std::vector<HalfFaceHandle>& _hfs) {

    Cell& c = cell(_ch);

    if(has_face_bottom_up_incidences()) {

        const std::vector<HalfFaceHandle>& hfs = c.halffaces();
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {

            incident_cell_per_hf_[*hf_it] = InvalidCellHandle;
        }

        for(std::vector<HalfFaceHandle>::const_iterator hf_it = _hfs.begin(),
                hf_end = _hfs.end(); hf_it != hf_end; ++hf_it) {

            incident_cell_per_hf_[*hf_it] = _ch;
        }
    }

    c.set_halffaces(_hfs);
}

//========================================================================================

/**
 * \brief Delete vertex from mesh
 *
 * Get all incident higher-dimensional entities and delete the complete
 * subtree of the mesh incident to vertex _h.
 * In this function all incident entities are gathered
 * and deleted using the delete_*_core functions
 * that do the actual deletion including the update
 * of the bottom-up incidences, etc.
 *
 * @param _h The handle to the vertex to be deleted
 */
VertexIter TopologyKernel::delete_vertex(const VertexHandle& _h) {

    std::vector<VertexHandle> vs;
    vs.push_back(_h);

    std::set<EdgeHandle> incidentEdges_s;
    get_incident_edges(vs, incidentEdges_s);

    std::set<FaceHandle> incidentFaces_s;
    get_incident_faces(incidentEdges_s, incidentFaces_s);

    std::set<CellHandle> incidentCells_s;
    get_incident_cells(incidentFaces_s, incidentCells_s);

    // Delete cells
    for(std::set<CellHandle>::const_reverse_iterator c_it = incidentCells_s.rbegin(),
            c_end = incidentCells_s.rend(); c_it != c_end; ++c_it) {
        delete_cell_core(*c_it);
    }

    // Delete faces
    for(std::set<FaceHandle>::const_reverse_iterator f_it = incidentFaces_s.rbegin(),
            f_end = incidentFaces_s.rend(); f_it != f_end; ++f_it) {
        delete_face_core(*f_it);
    }

    // Delete edges
    for(std::set<EdgeHandle>::const_reverse_iterator e_it = incidentEdges_s.rbegin(),
            e_end = incidentEdges_s.rend(); e_it != e_end; ++e_it) {
        delete_edge_core(*e_it);
    }

    // Delete vertex
    return delete_vertex_core(_h);
}

//========================================================================================

/**
 * \brief Delete edge from mesh
 *
 * Get all incident higher-dimensional entities and delete the complete
 * subtree of the mesh incident to edge _h.
 * In this function all incident entities are gathered
 * and deleted using the delete_*_core functions
 * that do the actual deletion including the update
 * of the bottom-up incidences, etc.
 *
 * @param _h The handle to the edge to be deleted
 */
EdgeIter TopologyKernel::delete_edge(const EdgeHandle& _h) {

    std::vector<EdgeHandle> es;
    es.push_back(_h);

    std::set<FaceHandle> incidentFaces_s;
    get_incident_faces(es, incidentFaces_s);

    std::set<CellHandle> incidentCells_s;
    get_incident_cells(incidentFaces_s, incidentCells_s);

    // Delete cells
    for(std::set<CellHandle>::const_reverse_iterator c_it = incidentCells_s.rbegin(),
            c_end = incidentCells_s.rend(); c_it != c_end; ++c_it) {
        delete_cell_core(*c_it);
    }

    // Delete faces
    for(std::set<FaceHandle>::const_reverse_iterator f_it = incidentFaces_s.rbegin(),
            f_end = incidentFaces_s.rend(); f_it != f_end; ++f_it) {
        delete_face_core(*f_it);
    }

    // Delete edge
    return delete_edge_core(_h);
}

//========================================================================================

/**
 * \brief Delete face from mesh
 *
 * Get all incident higher-dimensional entities and delete the complete
 * subtree of the mesh incident to face _h.
 * In this function all incident entities are gathered
 * and deleted using the delete_*_core functions
 * that do the actual deletion including the update
 * of the bottom-up incidences, etc.
 *
 * @param _h The handle to the face to be deleted
 */
FaceIter TopologyKernel::delete_face(const FaceHandle& _h) {

    std::vector<FaceHandle> fs;
    fs.push_back(_h);

    std::set<CellHandle> incidentCells_s;
    get_incident_cells(fs, incidentCells_s);

    // Delete cells
    for(std::set<CellHandle>::const_reverse_iterator c_it = incidentCells_s.rbegin(),
            c_end = incidentCells_s.rend(); c_it != c_end; ++c_it) {
        delete_cell_core(*c_it);
    }

    // Delete face
    return delete_face_core(_h);
}

//========================================================================================

/**
 * \brief Delete cell from mesh
 *
 * Since there's no higher dimensional incident
 * entity to a cell, we can safely delete it from the
 * mesh.
 *
 * @param _h The handle to the cell to be deleted
 */
CellIter TopologyKernel::delete_cell(const CellHandle& _h) {

    return delete_cell_core(_h);
}

//========================================================================================

template <class ContainerT>
void TopologyKernel::get_incident_edges(const ContainerT& _vs,
                                        std::set<EdgeHandle>& _es) const {

    _es.clear();

    if(v_bottom_up_) {

        for(typename ContainerT::const_iterator v_it = _vs.begin(),
                v_end = _vs.end(); v_it != v_end; ++v_it) {

            const std::vector<HalfEdgeHandle>& inc_hes = outgoing_hes_per_vertex_[v_it->idx()];

            for(std::vector<HalfEdgeHandle>::const_iterator he_it = inc_hes.begin(),
                    he_end = inc_hes.end(); he_it != he_end; ++he_it) {

                _es.insert(edge_handle(*he_it));
            }
        }
    } else {

        for(typename ContainerT::const_iterator v_it = _vs.begin(),
                v_end = _vs.end(); v_it != v_end; ++v_it) {

            for(EdgeIter e_it = edges_begin(), e_end = edges_end(); e_it != e_end; ++e_it) {

                const Edge& e = edge(*e_it);

                if(e.from_vertex() == *v_it || e.to_vertex() == *v_it) {
                    _es.insert(*e_it);
                }
            }
        }
    }
}

//========================================================================================

template <class ContainerT>
void TopologyKernel::get_incident_faces(const ContainerT& _es,
                                        std::set<FaceHandle>& _fs) const {

    _fs.clear();

    if(e_bottom_up_) {

        for(typename ContainerT::const_iterator e_it = _es.begin(),
                e_end = _es.end(); e_it != e_end; ++e_it) {

            for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(halfedge_handle(*e_it, 0));
                    hehf_it.valid(); ++hehf_it) {

                const FaceHandle fh = face_handle(*hehf_it);

                if(_fs.count(fh) == 0) {
                    _fs.insert(fh);
                }
            }
        }
    } else {

        for(typename ContainerT::const_iterator e_it = _es.begin(),
                e_end = _es.end(); e_it != e_end; ++e_it) {

            for(FaceIter f_it = faces_begin(),
                    f_end = faces_end(); f_it != f_end; ++f_it) {

                const std::vector<HalfEdgeHandle>& hes = face(*f_it).halfedges();

                for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                        he_end = hes.end(); he_it != he_end; ++he_it) {

                    if(edge_handle(*he_it) == *e_it) {
                        _fs.insert(*f_it);
                        break;
                    }
                }
            }
        }
    }
}

//========================================================================================

template <class ContainerT>
void TopologyKernel::get_incident_cells(const ContainerT& _fs,
                                        std::set<CellHandle>& _cs) const {

    _cs.clear();

    if(f_bottom_up_) {

        for(typename ContainerT::const_iterator f_it = _fs.begin(),
            f_end = _fs.end(); f_it != f_end; ++f_it) {

            const HalfFaceHandle hfh0 = halfface_handle(*f_it, 0);
            const HalfFaceHandle hfh1 = halfface_handle(*f_it, 1);

            const CellHandle c0 = incident_cell(hfh0);
            const CellHandle c1 = incident_cell(hfh1);

            if(c0.is_valid()) _cs.insert(c0);
            if(c1.is_valid()) _cs.insert(c1);
        }
    } else {

        for(typename ContainerT::const_iterator f_it = _fs.begin(),
            f_end = _fs.end(); f_it != f_end; ++f_it) {

            for(CellIter c_it = cells_begin(), c_end = cells_end();
                c_it != c_end; ++c_it) {

                const std::vector<HalfFaceHandle>& hfs = cell(*c_it).halffaces();

                for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                        hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {

                    if(face_handle(*hf_it) == *f_it) {
                        _cs.insert(*c_it);
                        break;
                    }
                }
            }
        }
    }
}

//========================================================================================

/**
 * \brief Delete vertex from mesh
 *
 * After performing this operation, all vertices
 * following vertex _h in the array will be accessible
 * through their old handle decreased by one.
 * This function directly fixes the vertex links
 * in all edges. These steps are performed:
 *
 * 1) Decrease all vertex handles > _h in incident edges
 * 2) Delete entry in bottom-up list: V -> HE
 * 3) Delete vertex itself (not necessary here since
 *    a vertex is only represented by a number)
 * 4) Delete property entry
 *
 * @param _h A vertex's handle
 */
VertexIter TopologyKernel::delete_vertex_core(const VertexHandle& _h) {

    assert(_h.idx() < (int)n_vertices());

    // 1)
    if(v_bottom_up_) {

        // Decrease all vertex handles >= _h in all edge definitions
        for(int i = _h.idx(), end = n_vertices(); i < end; ++i) {
            const std::vector<HalfEdgeHandle>& hes = outgoing_hes_per_vertex_[i];
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                    he_end = hes.end(); he_it != he_end; ++he_it) {

                Edge& e = edge(edge_handle(*he_it));
                if(e.from_vertex().idx() == i) {
                    e.set_from_vertex(VertexHandle(i-1));
                }
                if(e.to_vertex().idx() == i) {
                    e.set_to_vertex(VertexHandle(i-1));
                }
            }
        }

    } else {

        // Iterate over all edges
        for(EdgeIter e_it = edges_begin(), e_end = edges_end();
                e_it != e_end; ++e_it) {

            // Decrease all vertex handles in edge definitions that are greater than _h
            if(edge(*e_it).from_vertex() > _h) {
                edge(*e_it).set_from_vertex(VertexHandle(edge(*e_it).from_vertex().idx() - 1));
            }
            if(edge(*e_it).to_vertex() > _h) {
                edge(*e_it).set_to_vertex(VertexHandle(edge(*e_it).to_vertex().idx() - 1));
            }
        }
    }

    // 2)
    if(v_bottom_up_) {
        assert(outgoing_hes_per_vertex_.size() > (unsigned int)_h.idx());
        outgoing_hes_per_vertex_.erase(outgoing_hes_per_vertex_.begin() + _h.idx());
    }

    // 3)
    --n_vertices_;

    // 4)
    vertex_deleted(_h);

    // Iterator to next element in vertex list
    return (vertices_begin() + _h.idx());
}

//========================================================================================

/**
 * \brief Delete edge from mesh
 *
 * After performing this operation, all edges
 * following edge _h in the array will be accessible
 * through their old handle decreased by one.
 * This function directly fixes the edge links
 * in all faces. These steps are performed:
 *
 * 1) Delete bottom-up links from incident vertices
 * 2) Decrease all half-edge handles > _h in incident faces
 * 3) Delete entry in bottom-up list: HE -> HF
 * 4) Decrease all half-edge handles > 2*_h.idx() in
 *    vertex bottom-up list
 * 5) Delete edge itself
 * 6) Delete property entry
 *
 * @param _h An edge's handle
 */
EdgeIter TopologyKernel::delete_edge_core(const EdgeHandle& _h) {

    assert(_h.idx() < (int)edges_.size());

    // 1)
    if(v_bottom_up_) {

        VertexHandle v0 = edge(_h).from_vertex();
        VertexHandle v1 = edge(_h).to_vertex();
        assert(outgoing_hes_per_vertex_.size() > (unsigned int)std::max(v0.idx(), v1.idx()));

        outgoing_hes_per_vertex_[v0.idx()].erase(
                std::remove(outgoing_hes_per_vertex_[v0.idx()].begin(),
                            outgoing_hes_per_vertex_[v0.idx()].end(),
                            halfedge_handle(_h, 0)),
                            outgoing_hes_per_vertex_[v0.idx()].end());

        outgoing_hes_per_vertex_[v1.idx()].erase(
                std::remove(outgoing_hes_per_vertex_[v1.idx()].begin(),
                            outgoing_hes_per_vertex_[v1.idx()].end(),
                            halfedge_handle(_h, 1)),
                            outgoing_hes_per_vertex_[v1.idx()].end());
    }

    // 2)
    if(e_bottom_up_) {

        assert(incident_hfs_per_he_.size() > (unsigned int)halfedge_handle(_h, 0).idx());

        // Decrease all half-edge handles > he and
        // delete all half-edge handles == he in face definitions
        // Get all faces that need updates
        std::set<FaceHandle> update_faces;
        for(std::vector<std::vector<HalfFaceHandle> >::const_iterator iit =
                (incident_hfs_per_he_.begin() + halfedge_handle(_h, 0).idx()),
                iit_end = incident_hfs_per_he_.end(); iit != iit_end; ++iit) {
            for(std::vector<HalfFaceHandle>::const_iterator it = iit->begin(),
                    end = iit->end(); it != end; ++it) {
                update_faces.insert(face_handle(*it));
            }
        }
        // Update respective handles
        HEHandleCorrection cor(halfedge_handle(_h, 1));
        for(std::set<FaceHandle>::iterator f_it = update_faces.begin(),
                f_end = update_faces.end(); f_it != f_end; ++f_it) {

            std::vector<HalfEdgeHandle> hes = face(*f_it).halfedges();

            // Delete current half-edge from face's half-edge list
            hes.erase(std::remove(hes.begin(), hes.end(), halfedge_handle(_h, 0)), hes.end());
            hes.erase(std::remove(hes.begin(), hes.end(), halfedge_handle(_h, 1)), hes.end());

            std::for_each(hes.begin(), hes.end(),
                          fun::bind(&HEHandleCorrection::correctValue, &cor, fun::placeholders::_1));
            face(*f_it).set_halfedges(hes);
        }
    } else {

        // Iterate over all faces
        for(FaceIter f_it = faces_begin(), f_end = faces_end();
                f_it != f_end; ++f_it) {

            // Get face's half-edges
            std::vector<HalfEdgeHandle> hes = face(*f_it).halfedges();

            // Delete current half-edge from face's half-edge list
            hes.erase(std::remove(hes.begin(), hes.end(), halfedge_handle(_h, 0)), hes.end());
            hes.erase(std::remove(hes.begin(), hes.end(), halfedge_handle(_h, 1)), hes.end());

            // Decrease all half-edge handles greater than _h in face
            HEHandleCorrection cor(halfedge_handle(_h, 1));
            std::for_each(hes.begin(), hes.end(),
                          fun::bind(&HEHandleCorrection::correctValue, &cor, fun::placeholders::_1));
            face(*f_it).set_halfedges(hes);
        }
    }

    // 3)
    if(e_bottom_up_) {
        assert(incident_hfs_per_he_.size() > (unsigned int)halfedge_handle(_h, 1).idx());
        incident_hfs_per_he_.erase(incident_hfs_per_he_.begin() + halfedge_handle(_h, 1).idx());
        incident_hfs_per_he_.erase(incident_hfs_per_he_.begin() + halfedge_handle(_h, 0).idx());
    }

    // 4)
    if(v_bottom_up_) {
        HEHandleCorrection cor(halfedge_handle(_h, 1));
        std::for_each(outgoing_hes_per_vertex_.begin(),
                      outgoing_hes_per_vertex_.end(),
                      fun::bind(&HEHandleCorrection::correctVecValue, &cor, fun::placeholders::_1));
    }

    // 5)
    edges_.erase(edges_.begin() + _h.idx());

    // 6)
    edge_deleted(_h);

    // Return iterator to next element in list
    return (edges_begin() + _h.idx());
}

//========================================================================================

/**
 * \brief Delete face from mesh
 *
 * After performing this operation, all faces
 * following face _h in the array will be accessible
 * through their old handle decreased by one.
 * This function directly fixes the face links
 * in all cells. These steps are performed:
 *
 * 1) Delete bottom-up links from incident edges
 * 2) Decrease all half-face handles > _h in incident cells
 * 3) Delete entry in bottom-up list: HF -> C
 * 4) Decrease all half-face handles > 2*_h.idx() in
 *    half-edge bottom-up list
 * 5) Delete face itself
 * 6) Delete property entry
 *
 * @param _h An face's handle
 */
FaceIter TopologyKernel::delete_face_core(const FaceHandle& _h) {

    assert(_h.idx() < (int)faces_.size());

    // 1)
    if(e_bottom_up_) {

        const std::vector<HalfEdgeHandle>& hes = face(_h).halfedges();
        for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                he_end = hes.end(); he_it != he_end; ++he_it) {

            assert(incident_hfs_per_he_.size() > (unsigned int)std::max(he_it->idx(), opposite_halfedge_handle(*he_it).idx()));

            incident_hfs_per_he_[he_it->idx()].erase(
                    std::remove(incident_hfs_per_he_[he_it->idx()].begin(),
                                incident_hfs_per_he_[he_it->idx()].end(),
                                halfface_handle(_h, 0)), incident_hfs_per_he_[he_it->idx()].end());


            incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].erase(
                    std::remove(incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].begin(),
                                incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].end(),
                                halfface_handle(_h, 1)), incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].end());
        }
    }

    // 2)
    if(f_bottom_up_) {

        // Decrease all half-face handles > _h in all cells
        // and delete all half-face handles == _h
        std::set<CellHandle> update_cells;
        for(std::vector<CellHandle>::const_iterator c_it = (incident_cell_per_hf_.begin() + halfface_handle(_h, 0).idx()),
                c_end = incident_cell_per_hf_.end(); c_it != c_end; ++c_it) {
            if(!c_it->is_valid()) continue;
            update_cells.insert(*c_it);
        }
        for(std::set<CellHandle>::const_iterator c_it = update_cells.begin(),
                c_end = update_cells.end(); c_it != c_end; ++c_it) {

            std::vector<HalfFaceHandle> hfs = cell(*c_it).halffaces();

            // Delete current half-faces from cell's half-face list
            hfs.erase(std::remove(hfs.begin(), hfs.end(), halfface_handle(_h, 0)), hfs.end());
            hfs.erase(std::remove(hfs.begin(), hfs.end(), halfface_handle(_h, 1)), hfs.end());

            HFHandleCorrection cor(halfface_handle(_h, 1));
            std::for_each(hfs.begin(), hfs.end(),
                          fun::bind(&HFHandleCorrection::correctValue, &cor, fun::placeholders::_1));
            cell(*c_it).set_halffaces(hfs);
        }

    } else {

        // Iterate over all cells
        for(CellIter c_it = cells_begin(), c_end = cells_end(); c_it != c_end; ++c_it) {

            std::vector<HalfFaceHandle> hfs = cell(*c_it).halffaces();

            // Delete current half-faces from cell's half-face list
            hfs.erase(std::remove(hfs.begin(), hfs.end(), halfface_handle(_h, 0)), hfs.end());
            hfs.erase(std::remove(hfs.begin(), hfs.end(), halfface_handle(_h, 1)), hfs.end());

            HFHandleCorrection cor(halfface_handle(_h, 1));
            std::for_each(hfs.begin(), hfs.end(),
                          fun::bind(&HFHandleCorrection::correctValue, &cor, fun::placeholders::_1));
            cell(*c_it).set_halffaces(hfs);
        }
    }

    // 3)
    if(f_bottom_up_) {
        assert(incident_cell_per_hf_.size() > (unsigned int)halfface_handle(_h, 1).idx());
        incident_cell_per_hf_.erase(incident_cell_per_hf_.begin() + halfface_handle(_h, 1).idx());
        incident_cell_per_hf_.erase(incident_cell_per_hf_.begin() + halfface_handle(_h, 0).idx());
    }

    // 4)
    if(e_bottom_up_) {
        HFHandleCorrection cor(halfface_handle(_h, 1));
        std::for_each(incident_hfs_per_he_.begin(),
                      incident_hfs_per_he_.end(),
                      fun::bind(&HFHandleCorrection::correctVecValue, &cor, fun::placeholders::_1));
    }

    // 5)
    faces_.erase(faces_.begin() + _h.idx());

    // 6)
    face_deleted(_h);

    // Return iterator to next element in list
    return (faces_begin() + _h.idx());
}

//========================================================================================

/**
 * \brief Delete cell from mesh
 *
 * After performing this operation, all cells
 * following cell _h in the array will be accessible
 * through their old handle decreased by one.
 * These steps are performed:
 *
 * 1) Delete links in BU: HF -> C
 * 2) Decrease all entries > c in BU: HF -> C
 * 3) Delete cell from storage array
 * 4) Delete property item
 *
 * @param _h A cell handle
 */
CellIter TopologyKernel::delete_cell_core(const CellHandle& _h) {

    assert(_h.idx() < (int)cells_.size());

    // 1)
    if(f_bottom_up_) {
        const std::vector<HalfFaceHandle>& hfs = cell(_h).halffaces();
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {
            assert(incident_cell_per_hf_.size() > (unsigned int)hf_it->idx());

            incident_cell_per_hf_[hf_it->idx()] = InvalidCellHandle;
        }
    }

    // 2)
    if(f_bottom_up_) {
        CHandleCorrection cor(_h);
        std::for_each(incident_cell_per_hf_.begin(),
                      incident_cell_per_hf_.end(),
                      fun::bind(&CHandleCorrection::correctValue, &cor, fun::placeholders::_1));
    }

    // 3)
    cells_.erase(cells_.begin() + _h.idx());

    // 4)
    cell_deleted(_h);

    return (cells_begin() + _h.idx());
}

//========================================================================================

void TopologyKernel::delete_multiple_vertices(const std::vector<bool>& _tag) {

    assert(_tag.size() == n_vertices_);

    std::vector<int> newIndices(n_vertices(), -1);
    int curIdx = 0;

    std::vector<int>::iterator idx_it = newIndices.begin();
    for(std::vector<bool>::const_iterator t_it = _tag.begin(),
            t_end = _tag.end(); t_it != t_end; ++t_it, ++idx_it) {

        if(!(*t_it)) {
            // Not marked as deleted
            *idx_it = curIdx;
            ++curIdx;
        } else {
            --n_vertices_;
        }
    }

    // Delete properties accordingly
    delete_multiple_vertex_props(_tag);

    EdgeCorrector corrector(newIndices);
    std::for_each(edges_.begin(), edges_.end(), corrector);
}

//========================================================================================

void TopologyKernel::delete_multiple_edges(const std::vector<bool>& _tag) {

    assert(_tag.size() == n_edges());

    std::vector<int> newIndices(n_edges(), -1);
    int curIdx = 0;

    std::vector<Edge> newEdges;

    std::vector<int>::iterator idx_it = newIndices.begin();
    std::vector<Edge>::const_iterator e_it = edges_.begin();

    for(std::vector<bool>::const_iterator t_it = _tag.begin(),
            t_end = _tag.end(); t_it != t_end; ++t_it, ++idx_it, ++e_it) {

        if(!(*t_it)) {
            // Not marked as deleted

            newEdges.push_back(*e_it);

            *idx_it = curIdx;
            ++curIdx;
        }
    }

    // Swap edges
    edges_.swap(newEdges);

    // Delete properties accordingly
    delete_multiple_edge_props(_tag);

    FaceCorrector corrector(newIndices);
    std::for_each(faces_.begin(), faces_.end(), corrector);
}

//========================================================================================

void TopologyKernel::delete_multiple_faces(const std::vector<bool>& _tag) {

    assert(_tag.size() == n_faces());

    std::vector<int> newIndices(n_faces(), -1);
    int curIdx = 0;

    std::vector<Face> newFaces;

    std::vector<int>::iterator idx_it = newIndices.begin();
    std::vector<Face>::const_iterator f_it = faces_.begin();

    for(std::vector<bool>::const_iterator t_it = _tag.begin(),
            t_end = _tag.end(); t_it != t_end; ++t_it, ++idx_it, ++f_it) {

        if(!(*t_it)) {
            // Not marked as deleted

            newFaces.push_back(*f_it);

            *idx_it = curIdx;
            ++curIdx;
        }
    }

    // Swap faces
    faces_.swap(newFaces);

    // Delete properties accordingly
    delete_multiple_face_props(_tag);

    CellCorrector corrector(newIndices);
    std::for_each(cells_.begin(), cells_.end(), corrector);
}

//========================================================================================

void TopologyKernel::delete_multiple_cells(const std::vector<bool>& _tag) {

    assert(_tag.size() == n_cells());

    std::vector<Cell> newCells;

    std::vector<Cell>::const_iterator c_it = cells_.begin();

    for(std::vector<bool>::const_iterator t_it = _tag.begin(),
            t_end = _tag.end(); t_it != t_end; ++t_it, ++c_it) {

        if(!(*t_it)) {
            // Not marked as deleted

            newCells.push_back(*c_it);
        }
    }

    // Swap cells
    cells_.swap(newCells);

    // Delete properties accordingly
    delete_multiple_cell_props(_tag);
}

//========================================================================================

CellIter TopologyKernel::delete_cell_range(const CellIter& _first, const CellIter& _last) {

    assert(_first >= cells_begin());
    assert(_last < cells_end());

    std::vector<Cell>::iterator it = cells_.erase(cells_.begin() + _first->idx(), cells_.begin() + _last->idx());

    // Re-compute face bottom-up incidences if necessary
    if(f_bottom_up_) {
        f_bottom_up_ = false;
        enable_face_bottom_up_incidences(true);
    }

    return CellIter(this, CellHandle(it - cells_.begin()));
}

//========================================================================================

/// Get edge with handle _edgeHandle
const OpenVolumeMeshEdge& TopologyKernel::edge(const EdgeHandle& _edgeHandle) const {

    // Test if edge is valid
    assert((unsigned int)_edgeHandle.idx() < edges_.size());
    assert(_edgeHandle.idx() >= 0);

    return edges_[_edgeHandle.idx()];
}

//========================================================================================

/// Get face with handle _faceHandle
const OpenVolumeMeshFace& TopologyKernel::face(const FaceHandle& _faceHandle) const {

    // Test if face is valid
    assert((unsigned int)_faceHandle.idx() < faces_.size());
    assert(_faceHandle.idx() >= 0);

    return faces_[_faceHandle.idx()];
}

//========================================================================================

/// Get cell with handle _cellHandle
const OpenVolumeMeshCell& TopologyKernel::cell(const CellHandle& _cellHandle) const {

    // Test if cell is valid
    assert((unsigned int)_cellHandle.idx() < cells_.size());
    assert(_cellHandle.idx() >= 0);

    return cells_[_cellHandle.idx()];
}

//========================================================================================

/// Get edge with handle _edgeHandle
OpenVolumeMeshEdge& TopologyKernel::edge(const EdgeHandle& _edgeHandle) {

    // Test if edge is valid
    assert((unsigned int)_edgeHandle.idx() < edges_.size());
    assert(_edgeHandle.idx() >= 0);

    return edges_[_edgeHandle.idx()];
}

//========================================================================================

/// Get face with handle _faceHandle
OpenVolumeMeshFace& TopologyKernel::face(const FaceHandle& _faceHandle) {

    // Test if face is valid
    assert((unsigned int)_faceHandle.idx() < faces_.size());
    assert(_faceHandle.idx() >= 0);

    return faces_[_faceHandle.idx()];
}

//========================================================================================

/// Get cell with handle _cellHandle
OpenVolumeMeshCell& TopologyKernel::cell(const CellHandle& _cellHandle) {

    // Test if cell is valid
    assert((unsigned int)_cellHandle.idx() < cells_.size());
    assert(_cellHandle.idx() >= 0);

    return cells_[_cellHandle.idx()];
}

//========================================================================================

/// Get edge that corresponds to halfedge with handle _halfEdgeHandle
OpenVolumeMeshEdge TopologyKernel::halfedge(const HalfEdgeHandle& _halfEdgeHandle) const {

    // Is handle in range?
    assert((unsigned int)_halfEdgeHandle.idx() < (edges_.size() * 2));
    assert(_halfEdgeHandle.idx() >= 0);

    // In case the handle is even, just return the corresponding edge
    /// Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle.idx() % 2 == 0)
        return edges_[(int)(_halfEdgeHandle.idx() / 2)];
    else
        return opposite_halfedge(edges_[(int)(_halfEdgeHandle.idx() / 2)]);
}

//========================================================================================

/// Get face that corresponds to halfface with handle _halfFaceHandle
OpenVolumeMeshFace TopologyKernel::halfface(const HalfFaceHandle& _halfFaceHandle) const {

    // Is handle in range?
    assert((unsigned int)_halfFaceHandle.idx() < (faces_.size() * 2));
    assert(_halfFaceHandle.idx() >= 0);

    // In case the handle is not even, just return the corresponding face
    // Otherwise return the opposite halfface via opposite()
    if(_halfFaceHandle.idx() % 2 == 0)
        return faces_[(int)(_halfFaceHandle.idx() / 2)];
    else
        return opposite_halfface(faces_[(int)(_halfFaceHandle.idx() / 2)]);
}

//========================================================================================

/// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
OpenVolumeMeshEdge TopologyKernel::opposite_halfedge(const HalfEdgeHandle& _halfEdgeHandle) const {

    // Is handle in range?
    assert(_halfEdgeHandle.idx() >= 0);
    assert((unsigned int)_halfEdgeHandle.idx() < (edges_.size() * 2));

    // In case the handle is not even, just return the corresponding edge
    // Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle.idx() % 2 != 0)
        return edges_[(int)(_halfEdgeHandle.idx() / 2)];
    else
        return opposite_halfedge(edges_[(int)(_halfEdgeHandle.idx() / 2)]);
}

//========================================================================================

/// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
OpenVolumeMeshFace TopologyKernel::opposite_halfface(const HalfFaceHandle& _halfFaceHandle) const {

    // Is handle in range?
    assert(_halfFaceHandle.idx() >= 0);
    assert((unsigned int)_halfFaceHandle.idx() < (faces_.size() * 2));

    // In case the handle is not even, just return the corresponding face
    // Otherwise return the opposite via the first face's opposite() function
    if(_halfFaceHandle.idx() % 2 != 0)
        return faces_[(int)(_halfFaceHandle.idx() / 2)];
    else
        return opposite_halfface(faces_[(int)(_halfFaceHandle.idx() / 2)]);
}

//========================================================================================

HalfEdgeHandle TopologyKernel::halfedge(const VertexHandle& _vh1, const VertexHandle& _vh2) const {

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

HalfFaceHandle TopologyKernel::halfface(const std::vector<VertexHandle>& _vs) const {

    assert(_vs.size() > 2);

    VertexHandle v0 = _vs[0], v1 = _vs[1], v2 = _vs[2];

    assert(v0.is_valid() && v1.is_valid() && v2.is_valid());

    HalfEdgeHandle he0 = halfedge(v0, v1);
    if(!he0.is_valid()) return InvalidHalfFaceHandle;
    HalfEdgeHandle he1 = halfedge(v1, v2);
    if(!he1.is_valid()) return InvalidHalfFaceHandle;

    std::vector<HalfEdgeHandle> hes;
    hes.push_back(he0);
    hes.push_back(he1);

    return halfface(hes);
}

//========================================================================================

HalfFaceHandle TopologyKernel::halfface(const std::vector<HalfEdgeHandle>& _hes) const {

    assert(_hes.size() >= 2);

    HalfEdgeHandle he0 = _hes[0], he1 = _hes[1];

    assert(he0.is_valid() && he1.is_valid());

    for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(he0); hehf_it.valid(); ++hehf_it) {

        std::vector<HalfEdgeHandle> hes = halfface(*hehf_it).halfedges();
        if(std::find(hes.begin(), hes.end(), he1) != hes.end()) {
            return *hehf_it;
        }
    }

    return InvalidHalfFaceHandle;
}

//========================================================================================

HalfEdgeHandle TopologyKernel::next_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const {

    assert((unsigned int)_hfh.idx() < faces_.size() * 2u);
    assert((unsigned int)_heh.idx() < edges_.size() * 2u);

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

HalfEdgeHandle TopologyKernel::prev_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const {

    assert((unsigned int)_hfh.idx() < faces_.size() * 2u);
    assert((unsigned int)_heh.idx() < edges_.size() * 2u);

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

HalfFaceHandle
TopologyKernel::adjacent_halfface_in_cell(const HalfFaceHandle& _halfFaceHandle, const HalfEdgeHandle& _halfEdgeHandle) const {

#ifndef NDEBUG
    if((unsigned int)_halfFaceHandle.idx() >= incident_cell_per_hf_.size() || _halfFaceHandle.idx() < 0) {
        return InvalidHalfFaceHandle;
    }
#endif
    if(!has_face_bottom_up_incidences()) {
        std::cerr << "Error: Function adjacent_halfface_in_cell() needs face bottom-up incidences!" << std::endl;
        return InvalidHalfFaceHandle;
    }
    if(incident_cell_per_hf_[_halfFaceHandle.idx()] == InvalidCellHandle) {
        // Specified halfface is on the outside of the complex
        return InvalidHalfFaceHandle;
    }

    OpenVolumeMeshCell c = cell(incident_cell_per_hf_[_halfFaceHandle.idx()]);

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

    if(!has_face_bottom_up_incidences()) {
        return InvalidCellHandle;
    }
    if((unsigned int)_halfFaceHandle.idx() >= incident_cell_per_hf_.size() || _halfFaceHandle.idx() < 0) {
        return InvalidCellHandle;
    }

    return incident_cell_per_hf_[_halfFaceHandle.idx()];
}

//========================================================================================

void TopologyKernel::compute_vertex_bottom_up_incidences() {

    // Clear incidences
    outgoing_hes_per_vertex_.clear();
    outgoing_hes_per_vertex_.resize(n_vertices());

    // Store outgoing halfedges per vertex
    unsigned int n_edges = edges_.size();
    for(unsigned int i = 0; i < n_edges; ++i) {

        VertexHandle from = edges_[i].from_vertex();
        if((unsigned int)from.idx() >= outgoing_hes_per_vertex_.size()) {
            std::cerr << "update_incidences(): Vertex handle is out of bounds!" << std::endl;
            return;
        }
        outgoing_hes_per_vertex_[from.idx()].push_back(halfedge_handle(EdgeHandle(i), 0));

        VertexHandle to = edges_[i].to_vertex();
        if((unsigned int)to.idx() >= outgoing_hes_per_vertex_.size()) {
            std::cerr << "update_incidences(): Vertex handle is out of bounds!" << std::endl;
            return;
        }
        // Store opposite halfedge handle
        outgoing_hes_per_vertex_[to.idx()].push_back(halfedge_handle(EdgeHandle(i), 1));
    }
}

//========================================================================================

void TopologyKernel::compute_edge_bottom_up_incidences() {

    // Clear
    incident_hfs_per_he_.clear();
    incident_hfs_per_he_.resize(edges_.size() * 2u);

    // Store incident halffaces per halfedge
    unsigned int n_faces = faces_.size();
    for(unsigned int i = 0; i < n_faces; ++i) {

        std::vector<HalfEdgeHandle> halfedges = faces_[i].halfedges();

        // Go over all halfedges
        for(std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();
                he_it != halfedges.end(); ++he_it) {

            incident_hfs_per_he_[he_it->idx()].push_back(halfface_handle(FaceHandle(i), 0));
            incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].push_back(
                    halfface_handle(FaceHandle(i), 1));
        }
    }
}

//========================================================================================

void TopologyKernel::compute_face_bottom_up_incidences() {

    // Clear
    incident_cell_per_hf_.clear();
    incident_cell_per_hf_.resize(faces_.size() * 2u, InvalidCellHandle);

    unsigned int n_cells = cells_.size();
    for(unsigned int i = 0; i < n_cells; ++i) {

        std::vector<HalfFaceHandle> halffaces = cells_[i].halffaces();

        // Go over all halffaces
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = halffaces.begin();
                hf_it != halffaces.end(); ++hf_it) {

            if(incident_cell_per_hf_[hf_it->idx()] == InvalidCellHandle) {

                incident_cell_per_hf_[hf_it->idx()] = CellHandle(i);

            } else {

                std::cerr << "Detected non-three-manifold configuration!" << std::endl;
                std::cerr << "Connectivity probably won't work." << std::endl;
                continue;
            }
        }
    }
}

} // Namespace OpenVolumeMesh
