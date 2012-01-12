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
 *   $Revision$                                                          *
 *   $Date$                   *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef HEXAHEDRALMESHITERATORS_HH
#define HEXAHEDRALMESHITERATORS_HH

#include "../../PolyhedralMesh/Utils/Iterators.hh"

namespace OpenVolumeMesh {

template <class VecT>
class HexahedralMesh;

template <class VecT>
class CellSheetCellIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::CellHandle,
	typename PolyhedralMesh<VecT>::CellHandle> {
private:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::CellHandle,
			typename PolyhedralMesh<VecT>::CellHandle>    BaseIter;
	typedef typename PolyhedralMesh<VecT>::CellHandle     CellHandle;
	typedef typename PolyhedralMesh<VecT>::HalfFaceHandle HalfFaceHandle;
public:
	CellSheetCellIter(const CellHandle& _ref_h, const unsigned char _orthDir,
			const HexahedralMesh<VecT>* _mesh);

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
	typename PolyhedralMesh<VecT>::HalfFaceHandle,
	typename PolyhedralMesh<VecT>::HalfFaceHandle> {
private:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::HalfFaceHandle,
			typename PolyhedralMesh<VecT>::HalfFaceHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::HalfFaceHandle HalfFaceHandle;
	typedef typename PolyhedralMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;
	typedef typename PolyhedralMesh<VecT>::EdgeHandle     EdgeHandle;
public:
	HalfFaceSheetHalfFaceIter(const HalfFaceHandle& _ref_h,
			const HexahedralMesh<VecT>* _mesh);
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
    typename PolyhedralMesh<VecT>::HalfFaceHandle,
    typename PolyhedralMesh<VecT>::HalfFaceHandle> {
private:
    typedef BaseIterator<VecT,
            typename PolyhedralMesh<VecT>::HalfFaceHandle,
            typename PolyhedralMesh<VecT>::HalfFaceHandle> BaseIter;
    typedef typename PolyhedralMesh<VecT>::HalfFaceHandle HalfFaceHandle;
    typedef typename PolyhedralMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;
    typedef typename PolyhedralMesh<VecT>::EdgeHandle     EdgeHandle;
public:
    OutsideNeighborHalfFaceIter(const HalfFaceHandle& _ref_h,
            const HexahedralMesh<VecT>* _mesh);
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

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(HEXAHEDRALMESHITERATORST_CC)
#include "IteratorsT.cc"
#endif

#endif /* HEXAHEDRALMESHITERATORS_HH */
