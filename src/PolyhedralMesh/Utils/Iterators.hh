/*
 * Iterators.hh
 *
 *  Created on: 24.06.2011
 *      Author: mike
 */

#ifndef ITERATORS_HH_
#define ITERATORS_HH_

#include <vector>
#include <iterator>

// Forward declaration
template <class VecT>
class OpenVolumeMesh;

template <class VecT,
class IH /*  Input handle type */,
class OH /* Output handle type */>
class BaseIterator {
public:

	// STL compliance
	typedef std::input_iterator_tag	iterator_category;
	typedef int						distance_type;
	typedef OH  					value_type;
	typedef OH* 					pointer;
	typedef OH& 					reference;

	BaseIterator(const OpenVolumeMesh<VecT>* _mesh, const IH& _ih, const OH& _ch = -1) :
		valid_(true), cur_handle_(_ch), ref_handle_(_ih), mesh_(_mesh) {}

	// STL compliance (needs to have default constructor)
	BaseIterator() : valid_(false) {}
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
	const OpenVolumeMesh<VecT>* mesh() const {
		return mesh_;
	}

private:

	bool valid_;
	OH cur_handle_;
	IH ref_handle_;
	const OpenVolumeMesh<VecT>* mesh_;
};

//===========================================================================

template <class VecT>
class VertexOHalfedgeIter :
	public BaseIterator<VecT,
	typename OpenVolumeMesh<VecT>::VertexHandle,
	typename OpenVolumeMesh<VecT>::HalfEdgeHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::VertexHandle,
			typename OpenVolumeMesh<VecT>::HalfEdgeHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::VertexHandle VertexHandle;

	VertexOHalfedgeIter(const VertexHandle& _vIdx,
			const OpenVolumeMesh<VecT>* _mesh);

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
	typename OpenVolumeMesh<VecT>::HalfEdgeHandle,
	typename OpenVolumeMesh<VecT>::HalfFaceHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::HalfEdgeHandle,
			typename OpenVolumeMesh<VecT>::HalfFaceHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;

	HalfEdgeHalfFaceIter(const HalfEdgeHandle& _heIdx, const OpenVolumeMesh<VecT>* _mesh);

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
	typename OpenVolumeMesh<VecT>::VertexHandle,
	typename OpenVolumeMesh<VecT>::CellHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::VertexHandle,
			typename OpenVolumeMesh<VecT>::CellHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::VertexHandle VertexHandle;
	typedef typename OpenVolumeMesh<VecT>::CellHandle CellHandle;
	typedef typename OpenVolumeMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;
	typedef typename OpenVolumeMesh<VecT>::HalfFaceHandle HalfFaceHandle;

	VertexCellIter(const VertexHandle& _vIdx, const OpenVolumeMesh<VecT>* _mesh);
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
	typename OpenVolumeMesh<VecT>::HalfEdgeHandle,
	typename OpenVolumeMesh<VecT>::CellHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::HalfEdgeHandle,
			typename OpenVolumeMesh<VecT>::CellHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;

	HalfedgeCellIter(const HalfEdgeHandle& _heIdx, const OpenVolumeMesh<VecT>* _mesh);

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
	typename OpenVolumeMesh<VecT>::CellHandle,
	typename OpenVolumeMesh<VecT>::VertexHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::CellHandle,
			typename OpenVolumeMesh<VecT>::VertexHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::CellHandle CellHandle;
	typedef typename OpenVolumeMesh<VecT>::VertexHandle VertexHandle;
	typedef typename OpenVolumeMesh<VecT>::HalfFaceHandle HalfFaceHandle;
	typedef typename OpenVolumeMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;

	CellVertexIter(const CellHandle& _cIdx, const OpenVolumeMesh<VecT>* _mesh);
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
	typename OpenVolumeMesh<VecT>::CellHandle,
	typename OpenVolumeMesh<VecT>::CellHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::CellHandle,
			typename OpenVolumeMesh<VecT>::CellHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::CellHandle CellHandle;
	typedef typename OpenVolumeMesh<VecT>::HalfFaceHandle HalfFaceHandle;

	CellCellIter(const CellHandle& _cIdx, const OpenVolumeMesh<VecT>* _mesh);
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
	typename OpenVolumeMesh<VecT>::FaceHandle,
	typename OpenVolumeMesh<VecT>::FaceHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::FaceHandle,
			typename OpenVolumeMesh<VecT>::FaceHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::FaceHandle FaceHandle;

	BoundaryFaceIter(const OpenVolumeMesh<VecT>* _mesh);

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
	typename OpenVolumeMesh<VecT>::VertexHandle,
	typename OpenVolumeMesh<VecT>::VertexHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::VertexHandle,
			typename OpenVolumeMesh<VecT>::VertexHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::VertexHandle VertexHandle;

	VertexIter(const OpenVolumeMesh<VecT>* _mesh, const VertexHandle& _vh = VertexHandle(0));

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
	typename OpenVolumeMesh<VecT>::EdgeHandle,
	typename OpenVolumeMesh<VecT>::EdgeHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::EdgeHandle,
			typename OpenVolumeMesh<VecT>::EdgeHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::EdgeHandle EdgeHandle;

	EdgeIter(const OpenVolumeMesh<VecT>* _mesh, const EdgeHandle& _eh = EdgeHandle(0));

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
	typename OpenVolumeMesh<VecT>::HalfEdgeHandle,
	typename OpenVolumeMesh<VecT>::HalfEdgeHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::HalfEdgeHandle,
			typename OpenVolumeMesh<VecT>::HalfEdgeHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::HalfEdgeHandle HalfEdgeHandle;

	HalfEdgeIter(const OpenVolumeMesh<VecT>* _mesh, const HalfEdgeHandle& _heh = HalfEdgeHandle(0));

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
	typename OpenVolumeMesh<VecT>::FaceHandle,
	typename OpenVolumeMesh<VecT>::FaceHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::FaceHandle,
			typename OpenVolumeMesh<VecT>::FaceHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::FaceHandle FaceHandle;

	FaceIter(const OpenVolumeMesh<VecT>* _mesh, const FaceHandle& _fh = FaceHandle(0));

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
	typename OpenVolumeMesh<VecT>::HalfFaceHandle,
	typename OpenVolumeMesh<VecT>::HalfFaceHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::HalfFaceHandle,
			typename OpenVolumeMesh<VecT>::HalfFaceHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::HalfFaceHandle HalfFaceHandle;

	HalfFaceIter(const OpenVolumeMesh<VecT>* _mesh, const HalfFaceHandle& _hfh = HalfFaceHandle(0));

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
	typename OpenVolumeMesh<VecT>::CellHandle,
	typename OpenVolumeMesh<VecT>::CellHandle> {
public:
	typedef BaseIterator<VecT,
			typename OpenVolumeMesh<VecT>::CellHandle,
			typename OpenVolumeMesh<VecT>::CellHandle> BaseIter;
	typedef typename OpenVolumeMesh<VecT>::CellHandle CellHandle;

	CellIter(const OpenVolumeMesh<VecT>* _mesh, const CellHandle& _ch = CellHandle(0));

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

#if defined(INCLUDE_TEMPLATES) && !defined(ITERATORST_CC)
#include "IteratorsT.cc"
#endif

#endif /* ITERATORS_HH_ */