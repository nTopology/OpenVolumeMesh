/*
 * OpenHexMeshIteratorsT.cc
 *
 *  Created on: Jun 20, 2011
 *      Author: kremer
 */

#define OPENHEXMESHITERATORST_CC

#include "OpenHexMeshIterators.hh"

//================================================================================================
// CellSheetCellIter
//================================================================================================

template <class VecT>
CellSheetCellIter<VecT>::CellSheetCellIter(const CellHandle& _ref_h,
        const unsigned char _orthDir, const OpenHexMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h) {

    if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

	// First off, get all surrounding cells
	std::vector<HalfFaceHandle> halffaces = _mesh->cell(_ref_h).halffaces();
	for(typename std::vector<HalfFaceHandle>::const_iterator hf_it = halffaces.begin();
			hf_it != halffaces.end(); ++hf_it) {
		// Add those, that are perpendicular to the specified _orthDir
		if(_mesh->orientation(*hf_it, _ref_h) != _orthDir &&
				_mesh->orientation(*hf_it, _ref_h) != OpenHexMesh<VecT>::opposite_orientation(_orthDir)) {
			CellHandle ch = _mesh->incident_cell_per_hf_[_mesh->opposite_halfface_handle(*hf_it)];
			if(ch != OpenVolumeMesh<VecT>::InvalidCellHandle) {
				neighb_sheet_cell_hs_.insert(ch);
			}
		}
	}

	cur_it_ = neighb_sheet_cell_hs_.begin();
	BaseIter::valid(cur_it_ != neighb_sheet_cell_hs_.end());
	if(BaseIter::valid()) {
		BaseIter::cur_handle(*cur_it_);
	}
}

template <class VecT>
CellSheetCellIter<VecT>& CellSheetCellIter<VecT>::operator--() {

    --cur_it_;
    if(cur_it_ >= neighb_sheet_cell_hs_.begin()) {
		BaseIter::cur_handle(*cur_it_);
	} else {
		BaseIter::valid(false);
	}
    return *this;
}

template <class VecT>
CellSheetCellIter<VecT>& CellSheetCellIter<VecT>::operator++() {

    ++cur_it_;
	if(cur_it_ != neighb_sheet_cell_hs_.end()) {
		BaseIter::cur_handle(*cur_it_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

//================================================================================================
// HalfFaceSheetHalfFaceIter
//================================================================================================

template <class VecT>
HalfFaceSheetHalfFaceIter<VecT>::HalfFaceSheetHalfFaceIter(const HalfFaceHandle& _ref_h,
        const OpenHexMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h) {

	if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

	/*
	 * Each halfface uniquely belongs to either a cell
	 * or the boundary. If the halfface belongs
	 * to a cell, it suffices to determine the local axis
	 * the halfface represents w.r.t. the cell and to
	 * iterate over all neighboring cells orthogonal to
	 * this direction. We have to find those halffaces
	 * of the neighboring cells that contain exactly one
	 * of the initial halfface's opposite halfedges.
	 */

	if(_mesh->is_boundary(_ref_h)) {
		std::cerr << "HalfFaceSheetHalfFaceIter: HalfFace is boundary!" << std::endl;
		BaseIter::valid(false);
        return;
	}

	typename OpenHexMesh<VecT>::CellHandle ch = _mesh->incident_cell(_ref_h);
	unsigned char orientation = _mesh->orientation(_ref_h, ch);
	std::vector<HalfEdgeHandle> hes_v = _mesh->halfface(_ref_h).opposite().halfedges();
	std::set<HalfEdgeHandle> hes;
	hes.insert(hes_v.begin(), hes_v.end());

	for(typename OpenHexMesh<VecT>::CellSheetCellIter csc_it = _mesh->csc_iter(ch, orientation);
			csc_it.valid(); ++csc_it) {

		std::vector<HalfFaceHandle> hfs = _mesh->cell(*csc_it).halffaces();
		for(typename std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin();
				hf_it != hfs.end(); ++hf_it) {

			std::vector<HalfEdgeHandle> hf_hes = _mesh->halfface(*hf_it).halfedges();
			for(typename std::vector<HalfEdgeHandle>::const_iterator he_it = hf_hes.begin();
					he_it != hf_hes.end(); ++he_it) {

				if(hes.count(*he_it) > 0) {
					// Found halfface that lies on the same sheet
					adjacent_halffaces_.push_back(*hf_it);
					common_edges_.push_back(_mesh->edge_handle(*he_it));
					break;
				}
			}
		}
	}

	cur_it_ = adjacent_halffaces_.begin();
	edge_it_ = common_edges_.begin();
    BaseIter::valid(cur_it_ != adjacent_halffaces_.end());
    if(BaseIter::valid()) {
    	BaseIter::cur_handle(*cur_it_);
    }
}

template <class VecT>
HalfFaceSheetHalfFaceIter<VecT>& HalfFaceSheetHalfFaceIter<VecT>::operator--() {

    --cur_it_;
    --edge_it_;
    if(cur_it_ >= adjacent_halffaces_.begin()) {
		BaseIter::cur_handle(*cur_it_);
	} else {
		BaseIter::valid(false);
	}
    return *this;
}

template <class VecT>
HalfFaceSheetHalfFaceIter<VecT>& HalfFaceSheetHalfFaceIter<VecT>::operator++() {

    ++cur_it_;
    ++edge_it_;
	if(cur_it_ != adjacent_halffaces_.end()) {
		BaseIter::cur_handle(*cur_it_);
	} else {
		BaseIter::valid(false);
	}
	return *this;
}

//================================================================================================
// OutsideNeighborHalfFaceIter
//================================================================================================

template <class VecT>
OutsideNeighborHalfFaceIter<VecT>::OutsideNeighborHalfFaceIter(const HalfFaceHandle& _ref_h,
        const OpenHexMesh<VecT>* _mesh) :
BaseIter(_mesh, _ref_h) {

    if(!_mesh->has_bottom_up_adjacencies()) {
        std::cerr << "This iterator needs bottom-up adjacencies!" << std::endl;
        BaseIter::valid(false);
        return;
    }

    // Go over all incident halfedges
    std::vector<HalfEdgeHandle> halfedges = _mesh->halfface(_ref_h).halfedges();
    for(typename std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();
            he_it != halfedges.end(); ++he_it) {

        // Get outside halffaces
        typename OpenVolumeMesh<VecT>::HalfEdgeHalfFaceIter hehf_it = _mesh->hehf_iter(_mesh->opposite_halfedge_handle(*he_it));
        for(; hehf_it.valid(); ++hehf_it) {

            if(_mesh->is_boundary(*hehf_it)) {
                neighbor_halffaces_.push_back(*hehf_it);
                common_edges_.push_back(_mesh->edge_handle(*he_it));
            }
        }
    }

    cur_it_ = neighbor_halffaces_.begin();
    edge_it_ = common_edges_.begin();
    BaseIter::valid(cur_it_ != neighbor_halffaces_.end());
    if(BaseIter::valid()) {
        BaseIter::cur_handle(*cur_it_);
    }
}

template <class VecT>
OutsideNeighborHalfFaceIter<VecT>& OutsideNeighborHalfFaceIter<VecT>::operator--() {

    --cur_it_;
    --edge_it_;
    if(cur_it_ >= neighbor_halffaces_.begin()) {
        BaseIter::cur_handle(*cur_it_);
    } else {
        BaseIter::valid(false);
    }
    return *this;
}

template <class VecT>
OutsideNeighborHalfFaceIter<VecT>& OutsideNeighborHalfFaceIter<VecT>::operator++() {

    ++cur_it_;
    ++edge_it_;
    if(cur_it_ != neighbor_halffaces_.end()) {
        BaseIter::cur_handle(*cur_it_);
    } else {
        BaseIter::valid(false);
    }
    return *this;
}
