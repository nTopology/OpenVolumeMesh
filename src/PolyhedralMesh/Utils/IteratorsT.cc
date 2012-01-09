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

#define ITERATORST_CC

#include <iostream>
#include <set>

#include "Iterators.hh"

namespace OpenVolumeMesh {

//================================================================================================
// VertexOHalfedgeIter
//================================================================================================

template <class VecT>
VertexOHalfedgeIter<VecT>::VertexOHalfedgeIter(const VertexHandle& _ref_h,
		const PolyhedralMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h),
cur_index_(0) {

	if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

	if((unsigned int)_ref_h >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
		BaseIter::valid(false);
	}

	if(BaseIter::valid()) {
		if((unsigned int)cur_index_ >= BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h].size()) {
			BaseIter::valid(false);
		}
	}

	if(BaseIter::valid()) {
		BaseIter::cur_handle((
				BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h])[cur_index_]);
	}
}

template <class VecT>
VertexOHalfedgeIter<VecT>& VertexOHalfedgeIter<VecT>::operator--() {

	--cur_index_;

	if(cur_index_ < 0) {
		BaseIter::valid(false);
	} else {
		BaseIter::cur_handle((
				BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle()])[cur_index_]);
	}

	return *this;
}

template <class VecT>
VertexOHalfedgeIter<VecT>& VertexOHalfedgeIter<VecT>::operator++() {

	++cur_index_;

	if((unsigned int)cur_index_ >= BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle()].size()) {
		BaseIter::valid(false);
	} else {
		BaseIter::cur_handle((
				BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle()])[cur_index_]);
	}

	return *this;
}

////================================================================================================
//// HalfEdgeHalfFaceIter
////================================================================================================

template <class VecT>
HalfEdgeHalfFaceIter<VecT>::HalfEdgeHalfFaceIter(const HalfEdgeHandle& _ref_h,
        const PolyhedralMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h),
cur_index_(0) {

	if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

	if((unsigned int)_ref_h >= BaseIter::mesh()->incident_hfs_per_he_.size()) {
		BaseIter::valid(false);
	}

	if(BaseIter::valid()) {
		if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[_ref_h].size()) {
			BaseIter::valid(false);
		}
	}

	if(BaseIter::valid()) {
		BaseIter::cur_handle((
				BaseIter::mesh()->incident_hfs_per_he_[_ref_h])[cur_index_]);
	}
}

template <class VecT>
HalfEdgeHalfFaceIter<VecT>& HalfEdgeHalfFaceIter<VecT>::operator--() {

	--cur_index_;

	if(cur_index_ < 0) {
		BaseIter::valid(false);
	} else {
		BaseIter::cur_handle((
				BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()])[cur_index_]);
	}

	return *this;
}

template <class VecT>
HalfEdgeHalfFaceIter<VecT>& HalfEdgeHalfFaceIter<VecT>::operator++() {

	++cur_index_;

	if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()].size()) {
		BaseIter::valid(false);
	} else {
		BaseIter::cur_handle((
				BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()])[cur_index_]);
	}

	return *this;
}

////================================================================================================
//// VertexCellIter
////================================================================================================

template <class VecT>
VertexCellIter<VecT>::VertexCellIter(const VertexHandle& _ref_h,
        const PolyhedralMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h) {

	if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
        BaseIter::valid(false);
        return;
    }

    // Build up cell list
    std::vector<HalfEdgeHandle> incidentHalfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h];
    for(typename std::vector<HalfEdgeHandle>::const_iterator it = incidentHalfedges.begin(); it != incidentHalfedges.end(); ++it) {

    	if(*it < 0 || (unsigned int)*it >= BaseIter::mesh()->incident_hfs_per_he_.size()) continue;
    	std::vector<HalfFaceHandle> incidentHalfFaces = BaseIter::mesh()->incident_hfs_per_he_[*it];

    	for(typename std::vector<HalfFaceHandle>::const_iterator hf_it = incidentHalfFaces.begin();
    			hf_it != incidentHalfFaces.end(); ++hf_it) {
    		if((unsigned int)*hf_it < BaseIter::mesh()->incident_cell_per_hf_.size()) {
    			CellHandle c_idx = BaseIter::mesh()->incident_cell_per_hf_[*hf_it];
    			if(c_idx != PolyhedralMesh<VecT>::InvalidCellHandle)
    				cells_.insert(c_idx);
    		}
    	}
    }
    cell_iter_ = cells_.begin();
    BaseIter::valid(cell_iter_ != cells_.end());
    if(BaseIter::valid()) {
    	BaseIter::cur_handle(*cell_iter_);
    }
}

template <class VecT>
VertexCellIter<VecT>& VertexCellIter<VecT>::operator--() {

	--cell_iter_;
	if(cell_iter_ >= cells_.begin()) {
		BaseIter::cur_handle(*cell_iter_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

template <class VecT>
VertexCellIter<VecT>& VertexCellIter<VecT>::operator++() {

	++cell_iter_;
	if(cell_iter_ != cells_.end()) {
		BaseIter::cur_handle(*cell_iter_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

////================================================================================================
//// HalfedgeCellIter
////================================================================================================

template <class VecT>
HalfedgeCellIter<VecT>::HalfedgeCellIter(const HalfEdgeHandle& _ref_h,
        const PolyhedralMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h),
cur_index_(0) {

	if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h >= BaseIter::mesh()->incident_hfs_per_he_.size()) {

        BaseIter::valid(false);
        return;
    }
    if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[_ref_h].size()) {

    	BaseIter::valid(false);
    	return;
    }
    if((unsigned int)((BaseIter::mesh()->incident_hfs_per_he_[_ref_h])[cur_index_]) >=
    		BaseIter::mesh()->incident_cell_per_hf_.size()) {

    	BaseIter::valid(false);
    	return;
    }

    BaseIter::cur_handle(BaseIter::mesh()->incident_cell_per_hf_[((BaseIter::mesh()->incident_hfs_per_he_[_ref_h])[cur_index_])]);
}

template <class VecT>
HalfedgeCellIter<VecT>& HalfedgeCellIter<VecT>::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
		return *this;
	}
	BaseIter::cur_handle(BaseIter::mesh()->incident_cell_per_hf_[(
			(BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()])[cur_index_])]);
	return *this;
}

template <class VecT>
HalfedgeCellIter<VecT>& HalfedgeCellIter<VecT>::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()].size()) {

    	BaseIter::valid(false);
    	return *this;
    }
    if((unsigned int)((BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()])[cur_index_]) >=
    		BaseIter::mesh()->incident_cell_per_hf_.size()) {

    	BaseIter::valid(false);
    	return *this;
    }
    BaseIter::cur_handle(BaseIter::mesh()->incident_cell_per_hf_[(
    		(BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()])[cur_index_])]);
	return *this;
}

////================================================================================================
//// CellVertexIter
////================================================================================================

template <class VecT>
CellVertexIter<VecT>::CellVertexIter(const CellHandle& _ref_h,
        const PolyhedralMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h) {

    typename std::vector<HalfFaceHandle>::const_iterator hf_iter = BaseIter::mesh()->cell(_ref_h).halffaces().begin();
    for(; hf_iter != BaseIter::mesh()->cell(_ref_h).halffaces().end(); ++hf_iter) {
        std::vector<HalfEdgeHandle> hes = BaseIter::mesh()->halfface(*hf_iter).halfedges();
        for(typename std::vector<HalfEdgeHandle>::const_iterator he_iter = hes.begin(); he_iter != hes.end(); ++he_iter) {
            incident_vertices_.insert(BaseIter::mesh()->halfedge(*he_iter).to_vertex());
        }
    }

    v_iter_ = incident_vertices_.begin();
    BaseIter::valid(v_iter_ != incident_vertices_.end());

    if(BaseIter::valid()) {
    	BaseIter::cur_handle(*v_iter_);
    }
}

template <class VecT>
CellVertexIter<VecT>& CellVertexIter<VecT>::operator--() {

	--v_iter_;
	if(v_iter_ >= incident_vertices_.begin()) {
		BaseIter::cur_handle(*v_iter_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

template <class VecT>
CellVertexIter<VecT>& CellVertexIter<VecT>::operator++() {

	++v_iter_;
	if(v_iter_ != incident_vertices_.end()) {
		BaseIter::cur_handle(*v_iter_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

////================================================================================================
//// CellCellIter
////================================================================================================

template <class VecT>
CellCellIter<VecT>::CellCellIter(const CellHandle& _ref_h,
        const PolyhedralMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h) {

    if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

	typename std::vector<HalfFaceHandle>::const_iterator hf_iter = BaseIter::mesh()->cell(_ref_h).halffaces().begin();
	for(; hf_iter != BaseIter::mesh()->cell(_ref_h).halffaces().end(); ++hf_iter) {

		HalfFaceHandle opp_hf = BaseIter::mesh()->opposite_halfface_handle(*hf_iter);
		CellHandle ch = BaseIter::mesh()->incident_cell_per_hf_[opp_hf];
		if(ch != PolyhedralMesh<VecT>::InvalidCellHandle) {
			adjacent_cells_.insert(ch);
		}
	}

	c_iter_ = adjacent_cells_.begin();
	BaseIter::valid(c_iter_ != adjacent_cells_.end());
	if(BaseIter::valid()) {
		BaseIter::cur_handle(*c_iter_);
	}
}

template <class VecT>
CellCellIter<VecT>& CellCellIter<VecT>::operator--() {

	--c_iter_;
	if(c_iter_ >= adjacent_cells_.begin()) {
		BaseIter::cur_handle(*c_iter_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

template <class VecT>
CellCellIter<VecT>& CellCellIter<VecT>::operator++() {

	++c_iter_;
	if(c_iter_ != adjacent_cells_.end()) {
		BaseIter::cur_handle(*c_iter_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

////================================================================================================
//// BoundaryFaceIter
////================================================================================================

template <class VecT>
BoundaryFaceIter<VecT>::BoundaryFaceIter(const PolyhedralMesh<VecT>* _mesh) :
BaseIter(_mesh, PolyhedralMesh<VecT>::InvalidFaceHandle) {

	if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

	bf_it_ = BaseIter::mesh()->boundary_faces_.begin();
	BaseIter::valid(bf_it_ != BaseIter::mesh()->boundary_faces_.end());
	if(BaseIter::valid()) {
		BaseIter::cur_handle(*bf_it_);
	}
}

template <class VecT>
BoundaryFaceIter<VecT>& BoundaryFaceIter<VecT>::operator--() {

	--bf_it_;
	if(bf_it_ >= BaseIter::mesh()->boundary_faces_.begin()) {
		BaseIter::cur_handle(*bf_it_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

template <class VecT>
BoundaryFaceIter<VecT>& BoundaryFaceIter<VecT>::operator++() {

	++bf_it_;
	if(bf_it_ != BaseIter::mesh()->boundary_faces_.end()) {
		BaseIter::cur_handle(*bf_it_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

////================================================================================================
//// VertexIter
////================================================================================================

template <class VecT>
VertexIter<VecT>::VertexIter(const PolyhedralMesh<VecT>* _mesh, const VertexHandle& _vh) :
BaseIter(_mesh, PolyhedralMesh<VecT>::InvalidVertexHandle, _vh),
cur_index_(_vh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->vertices_.size()) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(typename PolyhedralMesh<VecT>::VertexHandle(cur_index_));
	}
}

template <class VecT>
VertexIter<VecT>& VertexIter<VecT>::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::VertexHandle(cur_index_));
	return *this;
}

template <class VecT>
VertexIter<VecT>& VertexIter<VecT>::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->vertices_.size()) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::VertexHandle(cur_index_));
	return *this;
}

////================================================================================================
//// EdgeIter
////================================================================================================

template <class VecT>
EdgeIter<VecT>::EdgeIter(const PolyhedralMesh<VecT>* _mesh, const EdgeHandle& _eh) :
BaseIter(_mesh, PolyhedralMesh<VecT>::InvalidEdgeHandle, _eh),
cur_index_(_eh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size()) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(typename PolyhedralMesh<VecT>::EdgeHandle(cur_index_));
	}
}

template <class VecT>
EdgeIter<VecT>& EdgeIter<VecT>::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::EdgeHandle(cur_index_));
	return *this;
}

template <class VecT>
EdgeIter<VecT>& EdgeIter<VecT>::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size()) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::EdgeHandle(cur_index_));
	return *this;
}

////================================================================================================
//// HalfEdgeIter
////================================================================================================

template <class VecT>
HalfEdgeIter<VecT>::HalfEdgeIter(const PolyhedralMesh<VecT>* _mesh, const HalfEdgeHandle& _heh) :
BaseIter(_mesh, PolyhedralMesh<VecT>::InvalidHalfEdgeHandle, _heh),
cur_index_(_heh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size() * 2) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(typename PolyhedralMesh<VecT>::HalfEdgeHandle(cur_index_));
	}
}

template <class VecT>
HalfEdgeIter<VecT>& HalfEdgeIter<VecT>::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::HalfEdgeHandle(cur_index_));
	return *this;
}

template <class VecT>
HalfEdgeIter<VecT>& HalfEdgeIter<VecT>::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size() * 2) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::HalfEdgeHandle(cur_index_));
	return *this;
}

////================================================================================================
//// FaceIter
////================================================================================================

template <class VecT>
FaceIter<VecT>::FaceIter(const PolyhedralMesh<VecT>* _mesh, const FaceHandle& _fh) :
BaseIter(_mesh, PolyhedralMesh<VecT>::InvalidFaceHandle, _fh),
cur_index_(_fh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size()) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(typename PolyhedralMesh<VecT>::FaceHandle(cur_index_));
	}
}

template <class VecT>
FaceIter<VecT>& FaceIter<VecT>::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::FaceHandle(cur_index_));
	return *this;
}

template <class VecT>
FaceIter<VecT>& FaceIter<VecT>::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size()) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::FaceHandle(cur_index_));
	return *this;
}

////================================================================================================
//// HalfFaceIter
////================================================================================================

template <class VecT>
HalfFaceIter<VecT>::HalfFaceIter(const PolyhedralMesh<VecT>* _mesh, const HalfFaceHandle& _hfh) :
BaseIter(_mesh, PolyhedralMesh<VecT>::InvalidHalfFaceHandle, _hfh),
cur_index_(_hfh.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size() * 2) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(typename PolyhedralMesh<VecT>::HalfFaceHandle(cur_index_));
	}
}

template <class VecT>
HalfFaceIter<VecT>& HalfFaceIter<VecT>::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::HalfFaceHandle(cur_index_));
	return *this;
}

template <class VecT>
HalfFaceIter<VecT>& HalfFaceIter<VecT>::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size() * 2) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::HalfFaceHandle(cur_index_));
	return *this;
}

////================================================================================================
//// CellIter
////================================================================================================

template <class VecT>
CellIter<VecT>::CellIter(const PolyhedralMesh<VecT>* _mesh, const CellHandle& _ch) :
BaseIter(_mesh, PolyhedralMesh<VecT>::InvalidCellHandle, _ch),
cur_index_(_ch.idx()) {

	if((unsigned int)cur_index_ >= BaseIter::mesh()->cells_.size()) {
		BaseIter::valid(false);
	}
	if(BaseIter::valid()) {
		BaseIter::cur_handle(typename PolyhedralMesh<VecT>::CellHandle(cur_index_));
	}
}

template <class VecT>
CellIter<VecT>& CellIter<VecT>::operator--() {

	--cur_index_;
	if(cur_index_ < 0) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::CellHandle(cur_index_));
	return *this;
}

template <class VecT>
CellIter<VecT>& CellIter<VecT>::operator++() {

	++cur_index_;
	if((unsigned int)cur_index_ >= BaseIter::mesh()->cells_.size()) {
		BaseIter::valid(false);
	}
	BaseIter::cur_handle(typename PolyhedralMesh<VecT>::CellHandle(cur_index_));
	return *this;
}

} // Namespace OpenVolumeMesh
