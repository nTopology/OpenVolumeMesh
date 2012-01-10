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

#ifndef ITERATORS_HH_
#define ITERATORS_HH_

#include <vector>
#include <iterator>

namespace OpenVolumeMesh {

// Forward declaration
template <class VecT>
class PolyhedralMesh;

template <class VecT,
class IH /*  Input handle type */,
class OH /* Output handle type */>
class BaseIterator {
public:

	// STL compliance
	typedef std::input_iterator_tag	iterator_category;
	typedef int						          distance_type;
	typedef OH  					          value_type;
	typedef OH* 					          pointer;
	typedef OH& 					          reference;

	BaseIterator(const PolyhedralMesh<VecT>* _mesh, const IH& _ih = -1, const OH& _ch = -1) :
		valid_(true), cur_handle_(_ch), ref_handle_(_ih), mesh_(_mesh) {}

	// STL compliance (needs to have default constructor)
	BaseIterator() : valid_(false), mesh_(0) {}
	virtual ~BaseIterator() {}
	bool operator== (const BaseIterator& _c) const {
		return (this->cur_handle_ == _c.cur_handle() &&
				this->ref_handle_ == _c.ref_handle() &&
				this->mesh_ == _c.mesh());
	}
	bool operator!= (const BaseIterator& _c) const {
		return !this->operator==(_c);
	}

	const OH* operator->() const {
		return &cur_handle_;
	}

	const OH& operator*() const {
		return cur_handle_;
	}

	BaseIterator& operator=(const BaseIterator& _c) {
		this->valid_ = _c.valid();
		this->cur_handle_ = _c.cur_handle();
		this->ref_handle_ = _c.ref_handle();
		this->mesh_ = _c.mesh();
		return *this;
	}

	// Increment/decrement operators to be overloaded
//	virtual BaseIterator& operator++() 	  	 = 0;
//	virtual BaseIterator  operator++(int)	 = 0;
//	virtual BaseIterator& operator--() 	  	 = 0;
//	virtual BaseIterator  operator--(int)	 = 0;
//
//	virtual BaseIterator  operator+ (int _n) = 0;
//	virtual BaseIterator& operator+=(int _n) = 0;
//	virtual BaseIterator  operator- (int _n) = 0;
//	virtual BaseIterator& operator-=(int _n) = 0;

//	operator int() const {
//		return cur_handle_;
//	}

	operator bool() const {
		return valid_;
	}

	void valid(bool _valid) {
		valid_ = _valid;
	}
	bool valid() const {
		return valid_;
	}
	void cur_handle(const OH& _h) {
		cur_handle_ = _h;
	}
	const OH& cur_handle() const {
		return cur_handle_;
	}
	const IH& ref_handle() const {
		return ref_handle_;
	}
	const PolyhedralMesh<VecT>* mesh() const {
		return mesh_;
	}

private:

	bool valid_;
	OH cur_handle_;
	IH ref_handle_;
	const PolyhedralMesh<VecT>* mesh_;
};

//===========================================================================

template <class VecT>
class VertexOHalfedgeIter :
	public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::VertexHandle,
	typename PolyhedralMesh<VecT>::HalfEdgeHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::VertexHandle,
			typename PolyhedralMesh<VecT>::HalfEdgeHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::VertexHandle VertexHandle;

	VertexOHalfedgeIter(const VertexHandle& _vIdx,
			const PolyhedralMesh<VecT>* _mesh);

	// Post increment/decrement operator
	VertexOHalfedgeIter operator++(int) {
		VertexOHalfedgeIter cpy = *this;
		++(*this);
		return cpy;
	}
	VertexOHalfedgeIter operator--(int) {
		VertexOHalfedgeIter cpy = *this;
		--(*this);
		return cpy;
	}
	VertexOHalfedgeIter operator+(int _n) {
		VertexOHalfedgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	VertexOHalfedgeIter operator-(int _n) {
		VertexOHalfedgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	VertexOHalfedgeIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	VertexOHalfedgeIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	VertexOHalfedgeIter& operator++();
	VertexOHalfedgeIter& operator--();

private:

	int cur_index_;
};

//===========================================================================

template <class VecT> class HalfEdgeHalfFaceIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::HalfEdgeHandle,
	typename PolyhedralMesh<VecT>::HalfFaceHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::HalfEdgeHandle,
			typename PolyhedralMesh<VecT>::HalfFaceHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;

	HalfEdgeHalfFaceIter(const HalfEdgeHandle& _heIdx, const PolyhedralMesh<VecT>* _mesh);

	// Post increment/decrement operator
	HalfEdgeHalfFaceIter operator++(int) {
		HalfEdgeHalfFaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfEdgeHalfFaceIter operator--(int) {
		HalfEdgeHalfFaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfEdgeHalfFaceIter operator+(int _n) {
		HalfEdgeHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfEdgeHalfFaceIter operator-(int _n) {
		HalfEdgeHalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfEdgeHalfFaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfEdgeHalfFaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfEdgeHalfFaceIter& operator++();
	HalfEdgeHalfFaceIter& operator--();

private:
	int cur_index_;
};

//===========================================================================

template <class VecT> class VertexCellIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::VertexHandle,
	typename PolyhedralMesh<VecT>::CellHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::VertexHandle,
			typename PolyhedralMesh<VecT>::CellHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::VertexHandle VertexHandle;
	typedef typename PolyhedralMesh<VecT>::CellHandle CellHandle;
	typedef typename PolyhedralMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;
	typedef typename PolyhedralMesh<VecT>::HalfFaceHandle HalfFaceHandle;

	VertexCellIter(const VertexHandle& _vIdx, const PolyhedralMesh<VecT>* _mesh);
	VertexCellIter& operator=(const VertexCellIter& _c) {
		BaseIter::operator=(_c);
		cells_ = cells_;
		cell_iter_ = cells_.begin();
		return *this;
	}

	// Post increment/decrement operator
	VertexCellIter operator++(int) {
		VertexCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	VertexCellIter operator--(int) {
		VertexCellIter cpy = *this;
		--(*this);
		return cpy;
	}
	VertexCellIter operator+(int _n) {
		VertexCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	VertexCellIter operator-(int _n) {
		VertexCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	VertexCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	VertexCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	VertexCellIter& operator++();
	VertexCellIter& operator--();

private:
	std::set<CellHandle> cells_;
	typename std::set<CellHandle>::const_iterator cell_iter_;
};

template <class VecT> class HalfedgeCellIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::HalfEdgeHandle,
	typename PolyhedralMesh<VecT>::CellHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::HalfEdgeHandle,
			typename PolyhedralMesh<VecT>::CellHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;

	HalfedgeCellIter(const HalfEdgeHandle& _heIdx, const PolyhedralMesh<VecT>* _mesh);

	// Post increment/decrement operator
	HalfedgeCellIter operator++(int) {
		HalfedgeCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfedgeCellIter operator--(int) {
		HalfedgeCellIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfedgeCellIter operator+(int _n) {
		HalfedgeCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfedgeCellIter operator-(int _n) {
		HalfedgeCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfedgeCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfedgeCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfedgeCellIter& operator++();
	HalfedgeCellIter& operator--();

private:
	int cur_index_;
};

//===========================================================================

template <class VecT> class CellVertexIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::CellHandle,
	typename PolyhedralMesh<VecT>::VertexHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::CellHandle,
			typename PolyhedralMesh<VecT>::VertexHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::CellHandle CellHandle;
	typedef typename PolyhedralMesh<VecT>::VertexHandle VertexHandle;
	typedef typename PolyhedralMesh<VecT>::HalfFaceHandle HalfFaceHandle;
	typedef typename PolyhedralMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;

	CellVertexIter(const CellHandle& _cIdx, const PolyhedralMesh<VecT>* _mesh);
	CellVertexIter& operator=(const CellVertexIter& _c) {
		BaseIter::operator=(_c);
		incident_vertices_ = _c.incident_vertices_;
		v_iter_ = incident_vertices_.begin();
		return *this;
	}

	// Post increment/decrement operator
	CellVertexIter operator++(int) {
		CellVertexIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellVertexIter operator--(int) {
		CellVertexIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellVertexIter operator+(int _n) {
		CellVertexIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellVertexIter operator-(int _n) {
		CellVertexIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellVertexIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellVertexIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellVertexIter& operator++();
	CellVertexIter& operator--();

private:
	std::set<VertexHandle> incident_vertices_;
	typename std::set<VertexHandle>::const_iterator v_iter_;
};

//===========================================================================

template <class VecT> class CellCellIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::CellHandle,
	typename PolyhedralMesh<VecT>::CellHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::CellHandle,
			typename PolyhedralMesh<VecT>::CellHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::CellHandle CellHandle;
	typedef typename PolyhedralMesh<VecT>::HalfFaceHandle HalfFaceHandle;

	CellCellIter(const CellHandle& _cIdx, const PolyhedralMesh<VecT>* _mesh);
	CellCellIter& operator=(const CellCellIter& _c) {
		BaseIter::operator=(_c);
		adjacent_cells_ = _c.adjacent_cells_;
		c_iter_ = adjacent_cells_.begin();
		return *this;
	}

	// Post increment/decrement operator
	CellCellIter operator++(int) {
		CellCellIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellCellIter operator--(int) {
		CellCellIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellCellIter operator+(int _n) {
		CellCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellCellIter operator-(int _n) {
		CellCellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellCellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellCellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellCellIter& operator++();
	CellCellIter& operator--();

private:
	std::set<CellHandle> adjacent_cells_;
	typename std::set<CellHandle>::const_iterator c_iter_;
};

//===========================================================================

template <class VecT> class BoundaryFaceIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::FaceHandle,
	typename PolyhedralMesh<VecT>::FaceHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::FaceHandle,
			typename PolyhedralMesh<VecT>::FaceHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::FaceHandle FaceHandle;

	BoundaryFaceIter(const PolyhedralMesh<VecT>* _mesh);

	// Post increment/decrement operator
	BoundaryFaceIter operator++(int) {
		BoundaryFaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	BoundaryFaceIter operator--(int) {
		BoundaryFaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	BoundaryFaceIter operator+(int _n) {
		BoundaryFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	BoundaryFaceIter operator-(int _n) {
		BoundaryFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	BoundaryFaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	BoundaryFaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	BoundaryFaceIter& operator++();
	BoundaryFaceIter& operator--();

private:
	typename std::vector<FaceHandle>::const_iterator bf_it_;
};

//===========================================================================

template <class VecT> class VertexIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::VertexHandle,
	typename PolyhedralMesh<VecT>::VertexHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::VertexHandle,
			typename PolyhedralMesh<VecT>::VertexHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::VertexHandle VertexHandle;

	VertexIter(const PolyhedralMesh<VecT>* _mesh, const VertexHandle& _vh = VertexHandle(0));

	// Post increment/decrement operator
	VertexIter operator++(int) {
		VertexIter cpy = *this;
		++(*this);
		return cpy;
	}
	VertexIter operator--(int) {
		VertexIter cpy = *this;
		--(*this);
		return cpy;
	}
	VertexIter operator+(int _n) {
		VertexIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	VertexIter operator-(int _n) {
		VertexIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	VertexIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	VertexIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	VertexIter& operator++();
	VertexIter& operator--();

private:
	int cur_index_;
};

//===========================================================================

template <class VecT> class EdgeIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::EdgeHandle,
	typename PolyhedralMesh<VecT>::EdgeHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::EdgeHandle,
			typename PolyhedralMesh<VecT>::EdgeHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::EdgeHandle EdgeHandle;

	EdgeIter(const PolyhedralMesh<VecT>* _mesh, const EdgeHandle& _eh = EdgeHandle(0));

	// Post increment/decrement operator
	EdgeIter operator++(int) {
		EdgeIter cpy = *this;
		++(*this);
		return cpy;
	}
	EdgeIter operator--(int) {
		EdgeIter cpy = *this;
		--(*this);
		return cpy;
	}
	EdgeIter operator+(int _n) {
		EdgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	EdgeIter operator-(int _n) {
		EdgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	EdgeIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	EdgeIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	EdgeIter& operator++();
	EdgeIter& operator--();

private:
	int cur_index_;
};

//===========================================================================

template <class VecT> class HalfEdgeIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::HalfEdgeHandle,
	typename PolyhedralMesh<VecT>::HalfEdgeHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::HalfEdgeHandle,
			typename PolyhedralMesh<VecT>::HalfEdgeHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;

	HalfEdgeIter(const PolyhedralMesh<VecT>* _mesh, const HalfEdgeHandle& _heh = HalfEdgeHandle(0));

	// Post increment/decrement operator
	HalfEdgeIter operator++(int) {
		HalfEdgeIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfEdgeIter operator--(int) {
		HalfEdgeIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfEdgeIter operator+(int _n) {
		HalfEdgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfEdgeIter operator-(int _n) {
		HalfEdgeIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfEdgeIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfEdgeIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfEdgeIter& operator++();
	HalfEdgeIter& operator--();

private:
	int cur_index_;
};

//===========================================================================

template <class VecT> class FaceIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::FaceHandle,
	typename PolyhedralMesh<VecT>::FaceHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::FaceHandle,
			typename PolyhedralMesh<VecT>::FaceHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::FaceHandle FaceHandle;

	FaceIter(const PolyhedralMesh<VecT>* _mesh, const FaceHandle& _fh = FaceHandle(0));

	// Post increment/decrement operator
	FaceIter operator++(int) {
		FaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	FaceIter operator--(int) {
		FaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	FaceIter operator+(int _n) {
		FaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	FaceIter operator-(int _n) {
		FaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	FaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	FaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	FaceIter& operator++();
	FaceIter& operator--();

private:
	int cur_index_;
};

//===========================================================================

template <class VecT> class HalfFaceIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::HalfFaceHandle,
	typename PolyhedralMesh<VecT>::HalfFaceHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::HalfFaceHandle,
			typename PolyhedralMesh<VecT>::HalfFaceHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::HalfFaceHandle HalfFaceHandle;

	HalfFaceIter(const PolyhedralMesh<VecT>* _mesh, const HalfFaceHandle& _hfh = HalfFaceHandle(0));

	// Post increment/decrement operator
	HalfFaceIter operator++(int) {
		HalfFaceIter cpy = *this;
		++(*this);
		return cpy;
	}
	HalfFaceIter operator--(int) {
		HalfFaceIter cpy = *this;
		--(*this);
		return cpy;
	}
	HalfFaceIter operator+(int _n) {
		HalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	HalfFaceIter operator-(int _n) {
		HalfFaceIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	HalfFaceIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	HalfFaceIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	HalfFaceIter& operator++();
	HalfFaceIter& operator--();

private:
	int cur_index_;
};

//===========================================================================

template <class VecT> class CellIter : public BaseIterator<VecT,
	typename PolyhedralMesh<VecT>::CellHandle,
	typename PolyhedralMesh<VecT>::CellHandle> {
public:
	typedef BaseIterator<VecT,
			typename PolyhedralMesh<VecT>::CellHandle,
			typename PolyhedralMesh<VecT>::CellHandle> BaseIter;
	typedef typename PolyhedralMesh<VecT>::CellHandle CellHandle;

	CellIter(const PolyhedralMesh<VecT>* _mesh, const CellHandle& _ch = CellHandle(0));

	// Post increment/decrement operator
	CellIter operator++(int) {
		CellIter cpy = *this;
		++(*this);
		return cpy;
	}
	CellIter operator--(int) {
		CellIter cpy = *this;
		--(*this);
		return cpy;
	}
	CellIter operator+(int _n) {
		CellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			++cpy;
		}
		return cpy;
	}
	CellIter operator-(int _n) {
		CellIter cpy = *this;
		for(int i = 0; i < _n; ++i) {
			--cpy;
		}
		return cpy;
	}
	CellIter& operator+=(int _n) {
		for(int i = 0; i < _n; ++i) {
			++(*this);
		}
		return *this;
	}
	CellIter& operator-=(int _n) {
		for(int i = 0; i < _n; ++i) {
			--(*this);
		}
		return *this;
	}

	CellIter& operator++();
	CellIter& operator--();

private:
	int cur_index_;
};

//===========================================================================

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(ITERATORST_CC)
#include "IteratorsT.cc"
#endif

#endif /* ITERATORS_HH_ */
