/*
 * OpenHexMeshIterators.hh
 *
 *  Created on: 27.06.2011
 *      Author: mike
 */

#ifndef OPENHEXMESHITERATORS_HH_
#define OPENHEXMESHITERATORS_HH_

//#include "../OpenHexMesh.hh"
template <class VecT>
class OpenHexMesh;

#include "../../OpenVolumeMesh/Utils/Iterators.hh"

template <class VecT>
class CellSheetCellIter : public BaseIterator<VecT,
	typename OpenVolumeMesh<VecT>::CellHandle,
	typename OpenVolumeMesh<VecT>::CellHandle> {
private:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::CellHandle,
			typename OpenVolumeMesh<VecT>::CellHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::CellHandle CellHandle;
	typedef typename OpenVolumeMesh<VecT>::HalfFaceHandle HalfFaceHandle;
public:
	CellSheetCellIter(const CellHandle& _ref_h, const unsigned char _orthDir,
			const OpenHexMesh<VecT>* _mesh);

	CellSheetCellIter& operator=(const CellSheetCellIter& _c) {
		BaseIter::operator=(_c);
		neighb_sheet_cell_hs_ = _c.neighb_sheet_cell_hs_;
		cur_it_ = neighb_sheet_cell_hs_.begin();
		return *this;
	}

	// Post increment/decrement operator
	CellSheetCellIter operator++(int) {
		CellSheetCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellSheetCellIter operator--(int) {
		CellSheetCellIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellSheetCellIter operator+(int _n) {
		CellSheetCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellSheetCellIter operator-(int _n) {
		CellSheetCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellSheetCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellSheetCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellSheetCellIter& operator++();
	CellSheetCellIter& operator--();

private:
	std::set<CellHandle> neighb_sheet_cell_hs_;
	typename std::set<CellHandle>::const_iterator cur_it_;
};

template <class VecT>
class HalfFaceSheetHalfFaceIter : public BaseIterator<VecT,
	typename OpenVolumeMesh<VecT>::HalfFaceHandle,
	typename OpenVolumeMesh<VecT>::HalfFaceHandle> {
private:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::HalfFaceHandle,
			typename OpenVolumeMesh<VecT>::HalfFaceHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::HalfFaceHandle HalfFaceHandle;
	typedef typename OpenVolumeMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;
	typedef typename OpenVolumeMesh<VecT>::EdgeHandle     EdgeHandle;
public:
	HalfFaceSheetHalfFaceIter(const HalfFaceHandle& _ref_h,
			const OpenHexMesh<VecT>* _mesh);
	HalfFaceSheetHalfFaceIter& operator=(const HalfFaceSheetHalfFaceIter& _c) {
		BaseIter::operator=(_c);
		adjacent_halffaces_ = _c.adjacent_halffaces_;
		cur_it_ = adjacent_halffaces_.begin();
		return *this;
	}

	// Post increment/decrement operator
	HalfFaceSheetHalfFaceIter operator++(int) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator--(int) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator+(int _n) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfFaceSheetHalfFaceIter operator-(int _n) {
		HalfFaceSheetHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfFaceSheetHalfFaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfFaceSheetHalfFaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfFaceSheetHalfFaceIter& operator++();
	HalfFaceSheetHalfFaceIter& operator--();

	const EdgeHandle& common_edge() const { return *edge_it_; }

private:
	std::vector<HalfFaceHandle> adjacent_halffaces_;
	typename std::vector<HalfFaceHandle>::const_iterator cur_it_;
	std::vector<EdgeHandle> common_edges_;
    typename std::vector<EdgeHandle>::const_iterator edge_it_;
};

template <class VecT>
class OutsideNeighborHalfFaceIter : public BaseIterator<VecT,
    typename OpenVolumeMesh<VecT>::HalfFaceHandle,
    typename OpenVolumeMesh<VecT>::HalfFaceHandle> {
private:
    typedef BaseIterator<VecT,
            typename OpenVolumeMesh<VecT>::HalfFaceHandle,
            typename OpenVolumeMesh<VecT>::HalfFaceHandle> BaseIter;
    typedef typename OpenVolumeMesh<VecT>::HalfFaceHandle HalfFaceHandle;
    typedef typename OpenVolumeMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;
    typedef typename OpenVolumeMesh<VecT>::EdgeHandle     EdgeHandle;
public:
    OutsideNeighborHalfFaceIter(const HalfFaceHandle& _ref_h,
            const OpenHexMesh<VecT>* _mesh);
    OutsideNeighborHalfFaceIter& operator=(const OutsideNeighborHalfFaceIter& _c) {
        BaseIter::operator=(_c);
        neighbor_halffaces_ = _c.adjacent_halffaces_;
        cur_it_ = neighbor_halffaces_.begin();
        return *this;
    }

    // Post increment/decrement operator
    OutsideNeighborHalfFaceIter operator++(int) {
        OutsideNeighborHalfFaceIter cpy = *this;
        ++(*this);
        return cpy;
    }
    OutsideNeighborHalfFaceIter operator--(int) {
        OutsideNeighborHalfFaceIter cpy = *this;
        --(*this);
        return cpy;
    }
    OutsideNeighborHalfFaceIter operator+(int _n) {
        OutsideNeighborHalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    OutsideNeighborHalfFaceIter operator-(int _n) {
        OutsideNeighborHalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    OutsideNeighborHalfFaceIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    OutsideNeighborHalfFaceIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    const EdgeHandle& common_edge() const { return *edge_it_; }

    OutsideNeighborHalfFaceIter& operator++();
    OutsideNeighborHalfFaceIter& operator--();

private:
    std::vector<HalfFaceHandle> neighbor_halffaces_;
    std::vector<EdgeHandle> common_edges_;
    typename std::vector<HalfFaceHandle>::const_iterator cur_it_;
    typename std::vector<EdgeHandle>::const_iterator edge_it_;
};

#if defined(INCLUDE_TEMPLATES) && !defined(OPENHEXMESHITERATORST_CC)
#include "OpenHexMeshIteratorsT.cc"
#endif

#endif /* OPENHEXMESHITERATORS_HH_ */
