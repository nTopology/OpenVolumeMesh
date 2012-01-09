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

#ifndef POLYHEDRALMESH_HH
#define POLYHEDRALMESH_HH

//== INCLUDES =================================================================

#include <vector>
#include <set>

#include "../Core/OpenVolumeMeshHandle.hh"
#include "../Core/OpenVolumeMeshBaseKernel.hh"
#include "../Core/OpenVolumeMeshStatus.hh"
#include "Utils/Iterators.hh"

namespace OpenVolumeMesh {

//== FORWARD DECLARATIONS ======================================================

template <typename VecT>
class OpenVolumeMeshVertex;
template <typename VecT>
class OpenVolumeMeshEdge;
template <typename VecT>
class OpenVolumeMeshFace;
template <typename VecT>
class OpenVolumeMeshCell;

//== CLASS DEFINITION =========================================================

/** \class PolyhedralMesh PolyhedralMesh.hh

    This is a general index based data structure for easy storage and processing
    of volume meshes. Meshes are built defining entities of higher dimensions
    via their incident entities of lower dimension. The following properties
    are being stored intrinsically:

    Vertex (0-dimensional):
        - Position

    Edge (1-dimensional):
        - From vertex id
        - To vertex id

    Face (2-dimensional):
        - Sorted list of incident halfedges (in counter-clockwise order)

    Cell (3-dimensional):
        - List of incident half-faces

    Half-edges as well as half-faces are uniquely determined as follows:
    Each edge is implicitly given a direction via the from and to vertex
    indices. By swapping them the opposite half-edge can be fetched.
    Half-faces are uniquely determined by their incident half-edges.
    The opposite half-face consequently contains the opposite half-edges.

    Once the mesh is built, it is possible to compute the adjacency list
    by calling the update_adjacencies() function which computes all
    necessary adjacency relations in top-down order (with respect to
    the dimensionality of the entities). These adjacency relations are
    fundamental for the use of circulators which allow for comfortable
    addressing of incident/adjacent entities. The update_adjacencies()
    function will generate the following properties in addition to
    the ones mentioned above:

    Vertex (0-dimensional):
        - List of all outgoing half-edges

    Edge (1-dimensional):
        - List of all adjacent half faces (per half-edge)

    Face (2-dimensional):
        - Both incident cells (restriction to two incident cells in order
            (to preserve three manifoldness)

    Half-edge as well as half-face indices are computed as follows:
    The half-edge index of the first half-edge of edge k is 2*k,
    the half-edge index of the second half-edge of edge k is 2*k+1.
    This works analogously for half-faces.
*/

/*
 * TODO:
 *
 * - Implement garbage collection
 * - Make vertex points, edges, faces, etc.
 * 	 be default properties
 * - Make status be a dynamic property for each entity
 * - Put all OpenVolumeMesh stuff into its own namespace
 */

//***************************************************************************

template <typename VecT>
class PolyhedralMesh: public OpenVolumeMeshBaseKernel {
public:
    //=====================================================================
    // Defines
    //=====================================================================

	typedef class VertexHandle		VertexHandle;
	typedef class EdgeHandle		EdgeHandle;
	typedef class HalfEdgeHandle	HalfEdgeHandle;
	typedef class FaceHandle		FaceHandle;
	typedef class HalfFaceHandle	HalfFaceHandle;
	typedef class CellHandle		CellHandle;

    static const VertexHandle   InvalidVertexHandle;
    static const EdgeHandle     InvalidEdgeHandle;
    static const FaceHandle     InvalidFaceHandle;
    static const CellHandle     InvalidCellHandle;
    static const HalfEdgeHandle InvalidHalfEdgeHandle;
    static const HalfFaceHandle InvalidHalfFaceHandle;

    typedef OpenVolumeMeshVertex<VecT>  Vertex;
    typedef OpenVolumeMeshEdge<VecT>    Edge;
    typedef OpenVolumeMeshFace<VecT>    Face;
    typedef OpenVolumeMeshCell<VecT>    Cell;

    typedef std::vector<Vertex>  Vertices;
    typedef std::vector<Edge>    Edges;
    typedef std::vector<Face>    Faces;
    typedef std::vector<Cell>    Cells;

    //=====================================================================
    // Constructors/Destructor
    //=====================================================================

    PolyhedralMesh();
    virtual ~PolyhedralMesh() {};

    //=====================================================================
    // Iterators
    //=====================================================================

    friend class VertexOHalfedgeIter<VecT>;
    friend class HalfEdgeHalfFaceIter<VecT>;
    friend class VertexCellIter<VecT>;
    friend class HalfedgeCellIter<VecT>;
    friend class CellVertexIter<VecT>;
    friend class CellCellIter<VecT>;
    friend class BoundaryFaceIter<VecT>;
    friend class VertexIter<VecT>;
    friend class EdgeIter<VecT>;
    friend class HalfEdgeIter<VecT>;
    friend class FaceIter<VecT>;
    friend class HalfFaceIter<VecT>;
    friend class CellIter<VecT>;

    typedef class VertexOHalfedgeIter<VecT> 	VertexOHalfedgeIter;
    typedef class HalfEdgeHalfFaceIter<VecT> 	HalfEdgeHalfFaceIter;
    typedef class VertexCellIter<VecT> 		    VertexCellIter;
    typedef class HalfedgeCellIter<VecT> 		  HalfedgeCellIter;
    typedef class CellVertexIter<VecT> 		    CellVertexIter;
    typedef class CellCellIter<VecT> 			    CellCellIter;
    typedef class BoundaryFaceIter<VecT> 		  BoundaryFaceIter;
    typedef class VertexIter<VecT> 			      VertexIter;
    typedef class EdgeIter<VecT> 				      EdgeIter;
    typedef class HalfEdgeIter<VecT> 			    HalfEdgeIter;
    typedef class FaceIter<VecT> 				      FaceIter;
    typedef class HalfFaceIter<VecT> 			    HalfFaceIter;
    typedef class CellIter<VecT> 				      CellIter;

    VertexOHalfedgeIter voh_iter(const VertexHandle& _idx) const {
    	return VertexOHalfedgeIter(_idx, this);
    }

    HalfEdgeHalfFaceIter hehf_iter(const HalfEdgeHandle& _idx) const {
        return HalfEdgeHalfFaceIter(_idx, this);
    }

    VertexCellIter vc_iter(const VertexHandle& _idx) const {
        return VertexCellIter(_idx, this);
    }

    HalfedgeCellIter hec_iter(const HalfEdgeHandle& _idx) const {
        return HalfedgeCellIter(_idx, this);
    }

    CellVertexIter cv_iter(const CellHandle& _idx) const {
        return CellVertexIter(_idx, this);
    }

    CellCellIter cc_iter(const CellHandle& _idx) const {
        return CellCellIter(_idx, this);
    }

    BoundaryFaceIter bf_iter() const {
        return BoundaryFaceIter(this);
    }

    VertexIter v_iter() const {
        return VertexIter(this);
    }

    VertexIter vertices_begin() const {
        return VertexIter(this, VertexHandle(0));
    }

    VertexIter vertices_end() const {
        return VertexIter(this, VertexHandle(vertices_.size()));
    }

    EdgeIter e_iter() const {
        return EdgeIter(this);
    }

    EdgeIter edges_begin() const {
        return EdgeIter(this, EdgeHandle(0));
    }

    EdgeIter edges_end() const {
        return EdgeIter(this, EdgeHandle(edges_.size()));
    }

    HalfEdgeIter he_iter() const {
        return HalfEdgeIter(this);
    }

    HalfEdgeIter halfedges_begin() const {
        return HalfEdgeIter(this, HalfEdgeHandle(0));
    }

    HalfEdgeIter halfedges_end() const {
        return HalfEdgeIter(this, HalfEdgeHandle(edges_.size() * 2));
    }

    FaceIter f_iter() const {
        return FaceIter(this);
    }

    FaceIter faces_begin() const {
        return FaceIter(this, FaceHandle(0));
    }

    FaceIter faces_end() const {
        return FaceIter(this, FaceHandle(faces_.size()));
    }

    HalfFaceIter hf_iter() const {
        return HalfFaceIter(this);
    }

    HalfFaceIter halffaces_begin() const {
        return HalfFaceIter(this, HalfFaceHandle(0));
    }

    HalfFaceIter halffaces_end() const {
        return HalfFaceIter(this, HalfFaceHandle(faces_.size() * 2));
    }

    CellIter c_iter() const {
        return CellIter(this);
    }

    CellIter cells_begin() const {
        return CellIter(this, CellHandle(0));
    }

    CellIter cells_end() const {
        return CellIter(this, CellHandle(cells_.size()));
    }

    //=====================================================================
    // Public functions
    //=====================================================================

    /// Build adjacency list in top-down-order
    virtual void update_adjacencies();

    /// Add vertex
    virtual VertexHandle add_vertex(const VecT& _p);

    /// Add edge
    virtual EdgeHandle add_edge(const VertexHandle& _fromVertex, const VertexHandle& _toHandle);

    /// Add face via incident edges
    virtual FaceHandle add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck = true);

    /// Add face via incident vertices
    virtual FaceHandle add_face(const std::vector<VertexHandle>& _vertices);

    /// Add cell via incident halffaces
    virtual CellHandle add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck = true);

    /// Get vertex with handle _vertexHandle
    inline const OpenVolumeMeshVertex<VecT>& vertex(const VertexHandle& _vertexHandle) const;

    /// Get edge with handle _edgeHandle
    inline const OpenVolumeMeshEdge<VecT>& edge(const EdgeHandle& _edgeHandle) const;

    /// Get face with handle _faceHandle
    inline const OpenVolumeMeshFace<VecT>& face(const FaceHandle& _faceHandle) const;

    /// Get cell with handle _cellHandle
    inline const OpenVolumeMeshCell<VecT>& cell(const CellHandle& _cellHandle) const;

    /// Get vertex with handle _vertexHandle
    inline OpenVolumeMeshVertex<VecT>& vertex(const VertexHandle& _vertexHandle);

    /// Get edge with handle _edgeHandle
    inline OpenVolumeMeshEdge<VecT>& edge(const EdgeHandle& _edgeHandle);

    /// Get face with handle _faceHandle
    inline OpenVolumeMeshFace<VecT>& face(const FaceHandle& _faceHandle);

    /// Get cell with handle _cellHandle
    inline OpenVolumeMeshCell<VecT>& cell(const CellHandle& _cellHandle);

    /// Get edge that corresponds to halfedge with handle _halfEdgeHandle
    const OpenVolumeMeshEdge<VecT> halfedge(const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get face that corresponds to halfface with handle _halfFaceHandle
    const OpenVolumeMeshFace<VecT> halfface(const HalfFaceHandle& _halfFaceHandle) const;

    /// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
    const OpenVolumeMeshEdge<VecT> opposite_halfedge(const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
    const OpenVolumeMeshFace<VecT> opposite_halfface(const HalfFaceHandle& _halfFaceHandle) const;

    // Get halfedge from vertex _vh1 to _vh2
    const HalfEdgeHandle halfedge(const VertexHandle& _vh1, const VertexHandle& _vh2) const;

    /// Get next halfedge within a halfface
    const HalfEdgeHandle next_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const;

    /// Get previous halfedge within a halfface
    const HalfEdgeHandle prev_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const;

    /// Get valence of vertex
    inline unsigned int valence(const VertexHandle& _vh) const {
        if(!has_bottom_up_adjacencies_) {
            std::cerr << "Could not get vertex valence: No bottom-up adjacencies available!" << std::endl;
            return 0u;
        }
        assert((unsigned int)_vh < outgoing_hes_per_vertex_.size());
        return outgoing_hes_per_vertex_[_vh].size();
    }

    /// Get valence of edge
    inline unsigned int valence(const EdgeHandle& _eh) const {
        if(!has_bottom_up_adjacencies_) {
            std::cerr << "Could not get vertex valence: No bottom-up adjacencies available!" << std::endl;
            return 0u;
        }
        assert((unsigned int)halfedge_handle(_eh, 0) < incident_hfs_per_he_.size());
        return incident_hfs_per_he_[halfedge_handle(_eh, 0)].size();
    }

    /// Get valence of face
    inline unsigned int valence(const FaceHandle& _fh) const {

        assert((unsigned int)_fh < faces_.size());
        return face(_fh).halfedges().size();
    }

    /// Get valence of cell
    inline unsigned int valence(const CellHandle& _ch) const {

        assert((unsigned int)_ch < cells_.size());
        return cell(_ch).halffaces().size();
    }

    //=====================================================================
    // Change entities
    //=====================================================================

    bool set_vertex(const VertexHandle& _vertexHandle, const VecT& _position);

    bool set_edge(const EdgeHandle& _edgeHandle, const VertexHandle& _fromVertex, const VertexHandle& _toVertex);

    bool set_face(const FaceHandle& _faceHandle, const std::vector<HalfEdgeHandle>& _halfedges);

    bool set_cell(const CellHandle& _cellHandle, const std::vector<HalfFaceHandle>& _halffaces);

    /// Clear whole volume mesh
    void clear() {

        vertices_.clear();
        edges_.clear();
        faces_.clear();
        cells_.clear();
        outgoing_hes_per_vertex_.clear();
        incident_hfs_per_he_.clear();
        incident_cell_per_hf_.clear();
        boundary_faces_.clear();
        face_normals_.clear();

        // Clear props
        vprops_clear();
        eprops_clear();
        heprops_clear();
        fprops_clear();
        hfprops_clear();
        cprops_clear();
        mprops_clear();
    }

    //=====================================================================
    // Connectivity
    //=====================================================================

    /// Get halfface that is adjacent (w.r.t. a common halfedge) within the same cell
    HalfFaceHandle adjacent_halfface_in_cell(const HalfFaceHandle& _halfFaceHandle, const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get cell that is incident to the given halfface
    CellHandle incident_cell(const HalfFaceHandle& _halfFaceHandle) const;

    //=====================================================================
    // Handle conversion
    //=====================================================================

    /// Conversion function
    static inline HalfEdgeHandle halfedge_handle(const EdgeHandle& _idx, const unsigned char _subIdx) {
        // Is handle in range?
        if(_idx < 0 || _subIdx > 1u) return InvalidHalfEdgeHandle;
        return (2u * _idx) + (_subIdx ? 1u : 0u);
    }

    /// Conversion function
    static inline HalfFaceHandle halfface_handle(const FaceHandle& _idx, const unsigned char _subIdx) {
        // Is handle in range?
        if(_idx < 0 || _subIdx > 1u) return InvalidHalfFaceHandle;
        return (2u * _idx) + (_subIdx ? 1u : 0u);
    }

    /// Handle conversion
    static inline EdgeHandle edge_handle(const HalfEdgeHandle& _idx) {
        // Is handle in range?
        if(_idx < 0) return InvalidEdgeHandle;
        return (int)(_idx / 2u);
    }

    static inline FaceHandle face_handle(const HalfFaceHandle& _idx) {
        // Is handle in range?
        if(_idx < 0) return InvalidFaceHandle;
        return (int)(_idx / 2u);
    }

    static inline HalfEdgeHandle opposite_halfedge_handle(const HalfEdgeHandle& _idx) {
        // Is handle in range?
        if(_idx < 0) return InvalidHalfEdgeHandle;

        // Is handle even?
        if(_idx % 2u == 0u) {
            return _idx + 1u;
        }
        return _idx - 1u;
    }

    static inline HalfFaceHandle opposite_halfface_handle(const HalfFaceHandle& _idx) {
        // Is handle in range?
        if(_idx < 0) return InvalidHalfFaceHandle;

        // Is handle even?
        if(_idx % 2u == 0u) {
            return _idx + 1u;
        }
        return _idx - 1u;
    }

    //=====================================================================
    // Stats
    //=====================================================================

    unsigned int n_vertices() 	const { return vertices_.size(); }
    unsigned int n_edges()    	const { return edges_.size(); }
    unsigned int n_halfedges()	const { return edges_.size() * 2u; }
    unsigned int n_faces()    	const { return faces_.size(); }
    unsigned int n_halffaces()  const { return faces_.size() * 2u; }
    unsigned int n_cells()    	const { return cells_.size(); }

    unsigned int n_boundary_faces() const {
    	if(!has_bottom_up_adjacencies_) {
    		std::cerr << "Warning: This function needs bottom-up adjacencies!" << std::endl;
    		return 0;
    	}
    	return boundary_faces_.size();
    }

    bool is_boundary(const HalfFaceHandle& _halfFaceHandle) const {
    	return _halfFaceHandle >= 0 && (unsigned int)_halfFaceHandle < incident_cell_per_hf_.size() &&
    			incident_cell_per_hf_[_halfFaceHandle] == InvalidCellHandle;
    }

    bool is_boundary(const FaceHandle& _faceHandle) const {
    	return  is_boundary(halfface_handle(_faceHandle, 0)) ||
    			is_boundary(halfface_handle(_faceHandle, 1));
    }

    bool is_boundary(const EdgeHandle& _edgeHandle) const {
        if(!has_bottom_up_adjacencies()) {
            std::cerr << "Error: is_boundary() needs bottom-up adjacencies!" << std::endl;
            return false;
        }
        for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(halfedge_handle(_edgeHandle, 0));
                hehf_it.valid(); ++hehf_it) {
            if(is_boundary(face_handle(*hehf_it))) {
                return true;
            }
        }
        return false;
    }

    bool is_boundary(const HalfEdgeHandle& _halfedgeHandle) const {
        if(!has_bottom_up_adjacencies()) {
            std::cerr << "Error: is_boundary() needs bottom-up adjacencies!" << std::endl;
            return false;
        }
        for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(_halfedgeHandle);
                hehf_it.valid(); ++hehf_it) {
            if(is_boundary(face_handle(*hehf_it))) {
                return true;
            }
        }
        return false;
    }

    bool is_boundary(const VertexHandle& _vertexHandle) const {
        if(!has_bottom_up_adjacencies()) {
            std::cerr << "Error: is_boundary() needs bottom-up adjacencies!" << std::endl;
            return false;
        }
        for(VertexOHalfedgeIter voh_it = voh_iter(_vertexHandle); voh_it.valid(); ++voh_it) {
            if(is_boundary(*voh_it)) return true;
        }
        return false;
    }

    unsigned int n_vertices_in_cell(const CellHandle& _ch) const {

        std::set<VertexHandle> vertices;
        std::vector<HalfFaceHandle> hfs = cell(_ch).halffaces();
        for(typename std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin();
                hf_it != hfs.end(); ++hf_it) {
            std::vector<HalfEdgeHandle> hes = halfface(*hf_it).halfedges();
            for(typename std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin();
                he_it != hes.end(); ++he_it) {
                vertices.insert(halfedge(*he_it).to_vertex());
            }
        }
        return vertices.size();
    }

    // Delete entities
    void delete_vertex(const VertexHandle& _vh) {
        if(!has_vertex_status_) request_vertex_status();
        status(_vh).set_deleted(true);
    }

    void delete_edge(const EdgeHandle& _eh) {
        if(!has_edge_status_) request_edge_status();
        status(_eh).set_deleted(true);
    }

    void delete_face(const FaceHandle& _fh) {
        if(!has_face_status_) request_face_status();
        status(_fh).set_deleted(true);
    }

    void delete_cell(const CellHandle& _ch) {
        if(!has_cell_status_) request_cell_status();
        status(_ch).set_deleted(true);
    }

    // Really delete entities and reset handles in the
    // respective higher dimensional entities
    void garbage_collection();

    //=====================================================================
    // Additional properties
    //=====================================================================

    /// Request and compute face normals
    void request_face_normals();

    /// Just update all face normals
    void update_face_normals();

    /// Compute normal for one face
    void compute_face_normal(const FaceHandle& _fh);

    /// Get face normal
    const VecT& face_normal(const FaceHandle& _fh) const;

    /// Get halfface normal
    const VecT halfface_normal(const HalfFaceHandle& _hfh) const;

    /// Get status flag
    bool has_face_normals() const { return has_face_normals_; }

    /// Request vertex status property
    void request_vertex_status();

    /// Request edge status property
    void request_edge_status();

    /// Request halfedge status property
    void request_halfedge_status();

    /// Request face status property
    void request_face_status();

    /// Request halfface status property
    void request_halfface_status();

    /// Request cell status property
    void request_cell_status();

    /// Request status property for all entities
    void request_status();

    /// Get vertex status object
    OpenVolumeMeshStatus& status(const VertexHandle& _vh);

    /// Get edge status object
    OpenVolumeMeshStatus& status(const EdgeHandle& _eh);

    /// Get halfedge status object
    OpenVolumeMeshStatus& status(const HalfEdgeHandle& _heh);

    /// Get face status object
    OpenVolumeMeshStatus& status(const FaceHandle& _fh);

    /// Get halfface status object
    OpenVolumeMeshStatus& status(const HalfFaceHandle& _hfh);

    /// Get edge status object
    OpenVolumeMeshStatus& status(const CellHandle& _ch);

    // Get status flag
    bool has_vertex_status() const { return has_vertex_status_; }

    // Get status flag
    bool has_edge_status() const { return has_edge_status_; }

    // Get status flag
    bool has_halfedge_status() const { return has_halfedge_status_; }

    // Get status flag
    bool has_face_status() const { return has_face_status_; }

    // Get status flag
    bool has_halfface_status() const { return has_halfface_status_; }

    // Get status flag
    bool has_cell_status() const { return has_cell_status_; }

    // Get status flag
    bool has_bottom_up_adjacencies() const { return has_bottom_up_adjacencies_; }

private:

    //=====================================================================
    // Private member functions
    //=====================================================================

protected:

    //=====================================================================
    // Member variables
    //=====================================================================

    // List of vertices
    Vertices vertices_;

    // List of edges
    Edges edges_;

    // List of faces
    Faces faces_;

    // List of cells
    Cells cells_;

    //=======================
    // Bottom-up-adjacencies
    //=======================

    // Outgoing halfedges per vertex
    std::vector<std::vector<HalfEdgeHandle> > outgoing_hes_per_vertex_;

    // Incident halffaces per (directed) halfedge
    std::vector<std::vector<HalfFaceHandle> > incident_hfs_per_he_;

    // Incident cell (at most one) per halfface
    std::vector<CellHandle> incident_cell_per_hf_;

    // Store boundary faces
    std::vector<FaceHandle> boundary_faces_;

    // Store indicator whether bottom-up adjacencies
    // have been computed
    bool has_bottom_up_adjacencies_;

    //=======================
    // Additional properties
    //=======================

    // Face normals
    std::vector<VecT> face_normals_;

    // Flag that indicates whether face's have normals
    bool has_face_normals_;

    // Entity status container
    std::vector<OpenVolumeMeshStatus> vertex_status_;

    // Entity status container
    std::vector<OpenVolumeMeshStatus> edge_status_;

    // Entity status container
    std::vector<OpenVolumeMeshStatus> halfedge_status_;

    // Entity status container
    std::vector<OpenVolumeMeshStatus> face_status_;

    // Entity status container
    std::vector<OpenVolumeMeshStatus> halfface_status_;

    // Entity status container
    std::vector<OpenVolumeMeshStatus> cell_status_;

    // Flag that indicates whether vertex status is available
    bool has_vertex_status_;

    // Flag that indicates whether edge status is available
    bool has_edge_status_;

    // Flag that indicates whether halfedge status is available
    bool has_halfedge_status_;

    // Flag that indicates whether face status is available
    bool has_face_status_;

    // Flag that indicates whether halfface status is available
    bool has_halfface_status_;

    // Flag that indicates whether cell status is available
    bool has_cell_status_;
};

// Initialize constants
template <typename VecT>
const typename PolyhedralMesh<VecT>::VertexHandle   PolyhedralMesh<VecT>::InvalidVertexHandle(-1);
template <typename VecT>
const typename PolyhedralMesh<VecT>::EdgeHandle     PolyhedralMesh<VecT>::InvalidEdgeHandle(-1);
template <typename VecT>
const typename PolyhedralMesh<VecT>::FaceHandle     PolyhedralMesh<VecT>::InvalidFaceHandle(-1);
template <typename VecT>
const typename PolyhedralMesh<VecT>::CellHandle     PolyhedralMesh<VecT>::InvalidCellHandle(-1);
template <typename VecT>
const typename PolyhedralMesh<VecT>::HalfEdgeHandle PolyhedralMesh<VecT>::InvalidHalfEdgeHandle(-1);
template <typename VecT>
const typename PolyhedralMesh<VecT>::HalfFaceHandle PolyhedralMesh<VecT>::InvalidHalfFaceHandle(-1);

//***************************************************************************

template <typename VecT>
class OpenVolumeMeshVertex {
public:
    OpenVolumeMeshVertex(const VecT& _position) : position_(_position) {}

    virtual ~OpenVolumeMeshVertex() {}

    const VecT& position() const { return position_; }

    void set_position(const VecT& _position) { position_ = _position; }

private:
    VecT position_;
};

// Stream operator for vertices
template <typename VecT>
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshVertex<VecT>& _vertex) {
    return _os << "(" << _vertex.position() << ")";
}

//***************************************************************************

template <typename VecT>
class OpenVolumeMeshEdge {
public:
    OpenVolumeMeshEdge(const typename OpenVolumeMesh<VecT>::VertexHandle& _fromVertex,
                       const typename OpenVolumeMesh<VecT>::VertexHandle& _toVertex) :
        fromVertex_(_fromVertex),
        toVertex_(_toVertex) {}

    virtual ~OpenVolumeMeshEdge() {}

    typename OpenVolumeMesh<VecT>::VertexHandle from_vertex() const { return fromVertex_; }
    typename OpenVolumeMesh<VecT>::VertexHandle to_vertex()   const { return toVertex_;   }

    void set_from_vertex(const typename OpenVolumeMesh<VecT>::VertexHandle& _vertex) { fromVertex_ = _vertex; }
    void set_to_vertex(const typename OpenVolumeMesh<VecT>::VertexHandle& _vertex)   { toVertex_ = _vertex; }

    const OpenVolumeMeshEdge<VecT> opposite() const {
    	return OpenVolumeMeshEdge<VecT>(toVertex_, fromVertex_);
    }
private:
    typename OpenVolumeMesh<VecT>::VertexHandle fromVertex_;
    typename OpenVolumeMesh<VecT>::VertexHandle toVertex_;
};

// Stream operator for edges
template <typename VecT>
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshEdge<VecT>& _edge) {
    return _os << "(" << _edge.from_vertex() << ", " << _edge.to_vertex() << ")";
}

//***************************************************************************

template <typename VecT>
class OpenVolumeMeshFace {
public:
    OpenVolumeMeshFace(const std::vector<typename OpenVolumeMesh<VecT>::HalfEdgeHandle>& _halfedges) :
    halfedges_(_halfedges) {
        compute_opposite_halfedges();
    }

    virtual ~OpenVolumeMeshFace() {}

    const std::vector<typename OpenVolumeMesh<VecT>::HalfEdgeHandle>& halfedges() const { return halfedges_; }

    void set_halfedges(const std::vector<typename OpenVolumeMesh<VecT>::HalfEdgeHandle>& _halfedges) {
        halfedges_ = _halfedges;
        opp_halfedges_.clear();
        for(typename std::vector<typename OpenVolumeMesh<VecT>::HalfEdgeHandle>::const_iterator it = halfedges_.begin();
                it != halfedges_.end(); ++it) {
            opp_halfedges_.insert(opp_halfedges_.begin(), OpenVolumeMesh<VecT>::opposite_halfedge_handle(*it));
        }
    }

    const OpenVolumeMeshFace<VecT> opposite() const {
    	return OpenVolumeMeshFace<VecT>(opp_halfedges_);
    }

    void compute_opposite_halfedges() {
        opp_halfedges_.clear();
        for(typename std::vector<typename OpenVolumeMesh<VecT>::HalfEdgeHandle>::const_iterator it = halfedges_.begin();
                it != halfedges_.end(); ++it) {
            opp_halfedges_.insert(opp_halfedges_.begin(), OpenVolumeMesh<VecT>::opposite_halfedge_handle(*it));
        }
    }

private:
    std::vector<typename OpenVolumeMesh<VecT>::HalfEdgeHandle> halfedges_;
    std::vector<typename OpenVolumeMesh<VecT>::HalfEdgeHandle> opp_halfedges_;
};

// Stream operator for faces
template <typename VecT>
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshFace<VecT>& _face) {
    _os << "(";
    for(typename std::vector<typename OpenVolumeMesh<VecT>::HalfEdgeHandle>::const_iterator it = _face.halfedges().begin();
            it < _face.halfedges().end(); ++it) {
        _os << *it;
        if(it+1 < _face.halfedges().end())
            _os << ", ";
    }
    _os << ")";
    return _os;
}

//***************************************************************************

template <typename VecT>
class OpenVolumeMeshCell {
public:
    OpenVolumeMeshCell(const std::vector<typename OpenVolumeMesh<VecT>::HalfFaceHandle>& _halffaces) :
        halffaces_(_halffaces) {}

    virtual ~OpenVolumeMeshCell() {}

    const std::vector<typename OpenVolumeMesh<VecT>::HalfFaceHandle>& halffaces() const { return halffaces_; }

    void set_halffaces(const std::vector<typename OpenVolumeMesh<VecT>::HalfFaceHandle>& _halffaces) { halffaces_ = _halffaces; }

private:
    std::vector<typename OpenVolumeMesh<VecT>::HalfFaceHandle> halffaces_;
};

// Stream operator for cells
template <typename VecT>
std::ostream& operator<<(std::ostream& _os, const OpenVolumeMeshCell<VecT>& _cell) {
    _os << "(";
    for(typename std::vector<typename OpenVolumeMesh<VecT>::HalfFaceHandle>::const_iterator it = _cell.halffaces().begin();
            it < _cell.halffaces().end(); ++it) {
        _os << *it;
        if(it+1 < _cell.halffaces().end())
            _os << ", ";
    }
    _os << ")";
    return _os;
}

} // Namespace OpenVolumeMesh

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(POLYHEDRALMESHT_CC)
#include "PolyhedralMeshT.cc"
#endif
//=============================================================================

#endif // POLYHEDRALMESH_HH
