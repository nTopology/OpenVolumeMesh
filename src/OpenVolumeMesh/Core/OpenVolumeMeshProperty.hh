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
 *   $Revision$                                                          *
 *   $Date$                   *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef OPENVOLUMEMESHPROPERTY_HH
#define OPENVOLUMEMESHPROPERTY_HH

//== INCLUDES =================================================================

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "OpenVolumeMeshBaseProperty.hh"
#include "OpenVolumeMeshHandle.hh"

namespace OpenVolumeMesh {

//== CLASS DEFINITION =========================================================

/** \class OpenVolumeMeshPropertyT
 *
 *  \brief Default property class for any type T.
 *
 *  The default property class for any type T.
 *
 *  The property supports persistence if T is a "fundamental" type:
 *  - integer fundamental types except bool:
 *    char, short, int, long, long long (__int64 for MS VC++) and
 *    their unsigned companions.
 *  - float fundamentals except <tt>long double</tt>:
 *    float, double
 */

template<class T>
class OpenVolumeMeshPropertyT: public OpenVolumeMeshBaseProperty {
public:

	typedef T 										Value;
	typedef std::vector<T> 				vector_type;
	typedef T 										value_type;
	typedef typename vector_type::reference 		reference;
	typedef typename vector_type::const_reference 	const_reference;

public:

	/// Default constructor
	OpenVolumeMeshPropertyT(const std::string& _name = "<unknown>") :
		OpenVolumeMeshBaseProperty(_name) {
	}

	/// Copy constructor
	OpenVolumeMeshPropertyT(const OpenVolumeMeshPropertyT& _rhs) :
		OpenVolumeMeshBaseProperty(_rhs), data_(_rhs.data_) {
	}

public:
	// inherited from OpenVolumeMeshBaseProperty
	virtual void reserve(size_t _n) {
		data_.reserve(_n);
	}
	virtual void resize(size_t _n) {
		data_.resize(_n, Value());
	}
	virtual void clear() {
		data_.clear();
		vector_type().swap(data_);
	}
	virtual void push_back() {
		data_.push_back(T());
	}
	virtual void swap(size_t _i0, size_t _i1) {
		std::swap(data_[_i0], data_[_i1]);
	}

	void delete_element(size_t _idx) {
	    data_.erase(data_.begin() + _idx);
	}

	virtual void stats(std::ostream& _ostr) const {
		for(typename vector_type::const_iterator it = data_.begin();
			it != data_.end(); ++it) {
				_ostr << *it << " ";
		}
	}

public:

	virtual size_t n_elements() const {
		return data_.size();
	}
	virtual size_t element_size() const {
		return sizeof(T);
	}

#ifndef DOXY_IGNORE_THIS
	struct plus {
		size_t operator ()(size_t _b, const T& /*_v*/) {
			return _b + sizeof(T);
		}
	};
#endif

	virtual size_t size_of() const {
		if (element_size() != OpenVolumeMeshBaseProperty::UnknownSize)
			return this->OpenVolumeMeshBaseProperty::size_of(n_elements());
		return std::accumulate(data_.begin(), data_.end(), 0, plus());
	}

	virtual size_t size_of(size_t _n_elem) const {
		return this->OpenVolumeMeshBaseProperty::size_of(_n_elem);
	}

public:
	// data access interface

	/// Get pointer to array (does not work for T==bool)
	const T* data() const {

		if (data_.empty())
			return 0;

		return &data_[0];
	}

	/// Get reference to property vector (be careful, improper usage, e.g. resizing, may crash)
	vector_type& data_vector() {

		return data_;
	}

	/// Access the i'th element. No range check is performed!
	reference operator[](int _idx) {
		assert(size_t(_idx) < data_.size());
		return data_[_idx];
	}

	/// Const access to the i'th element. No range check is performed!
	const_reference operator[](int _idx) const {
		assert(size_t(_idx) < data_.size());
		return data_[_idx];
	}

	/// Make a copy of self.
	OpenVolumeMeshPropertyT<T>* clone() const {
		OpenVolumeMeshPropertyT<T>* p = new OpenVolumeMeshPropertyT<T> (*this);
		return p;
	}

private:

	vector_type data_;
};

//-----------------------------------------------------------------------------


/** \class OpenVolumeMeshPropertyT<bool>

 Property specialization for bool type.

 */
template<>
class OpenVolumeMeshPropertyT<bool> : public OpenVolumeMeshBaseProperty {
public:

	typedef std::vector<bool> 				vector_type;
	typedef bool 							value_type;
	typedef vector_type::reference 			reference;
	typedef vector_type::const_reference 	const_reference;

public:

	OpenVolumeMeshPropertyT(const std::string& _name = "<unknown>") :
		OpenVolumeMeshBaseProperty(_name) {
	}

public:
	// inherited from OpenVolumeMeshBaseProperty

	virtual void reserve(size_t _n) {
		data_.reserve(_n);
	}
	virtual void resize(size_t _n) {
		data_.resize(_n);
	}
	virtual void clear() {
		data_.clear();
		vector_type().swap(data_);
	}
	virtual void push_back() {
		data_.push_back(bool());
	}
	virtual void swap(size_t _i0, size_t _i1) {
		bool t(data_[_i0]);
		data_[_i0] = data_[_i1];
		data_[_i1] = t;
	}

	void delete_element(size_t _idx) {
        data_.erase(data_.begin() + _idx);
    }

public:

//	virtual void set_persistent(bool _yn) {
//		check_and_set_persistent<bool> (_yn);
//	}

	virtual size_t n_elements() const {
		return data_.size();
	}
	virtual size_t element_size() const {
		return OpenVolumeMeshBaseProperty::UnknownSize;
	}
	virtual size_t size_of() const {
		return size_of(n_elements());
	}
	virtual size_t size_of(size_t _n_elem) const {
		return _n_elem / 8 + ((_n_elem % 8) != 0);
	}

	virtual void stats(std::ostream& _ostr) const {
		for(vector_type::const_iterator it = data_.begin();
			it != data_.end(); ++it) {
				_ostr << (*it ? "true" : "false") << " ";
		}
	}

public:

	/// Access the i'th element. No range check is performed!
	reference operator[](int _idx) {
		assert(size_t(_idx) < data_.size());
		return data_[_idx];
	}

	/// Const access to the i'th element. No range check is performed!
	const_reference operator[](int _idx) const {
		assert(size_t(_idx) < data_.size());
		return data_[_idx];
	}

	/// Make a copy of self.
	OpenVolumeMeshPropertyT<bool>* clone() const {
		OpenVolumeMeshPropertyT<bool>* p = new OpenVolumeMeshPropertyT<bool> (
				*this);
		return p;
	}

private:

	vector_type data_;
};

//-----------------------------------------------------------------------------


/** \class OpenVolumeMeshPropertyT<std::string>

 Property specialization for std::string type.
 */
template<>
class OpenVolumeMeshPropertyT<std::string> : public OpenVolumeMeshBaseProperty {
public:

	typedef std::string 					Value;
	typedef std::vector<std::string> 		vector_type;
	typedef std::string 					value_type;
	typedef vector_type::reference 			reference;
	typedef vector_type::const_reference 	const_reference;

public:

	OpenVolumeMeshPropertyT(const std::string& _name = "<unknown>") :
		OpenVolumeMeshBaseProperty(_name) {
	}

public:
	// inherited from OpenVolumeMeshBaseProperty

	virtual void reserve(size_t _n) {
		data_.reserve(_n);
	}
	virtual void resize(size_t _n) {
		data_.resize(_n);
	}
	virtual void clear() {
		data_.clear();
		vector_type().swap(data_);
	}
	virtual void push_back() {
		data_.push_back(std::string());
	}
	virtual void swap(size_t _i0, size_t _i1) {
		std::swap(data_[_i0], data_[_i1]);
	}

	virtual void delete_element(size_t _idx) {
        data_.erase(data_.begin() + _idx);
    }

public:

//	virtual void set_persistent(bool _yn) {
//		check_and_set_persistent<std::string> (_yn);
//	}

	virtual size_t n_elements() const {
		return data_.size();
	}
	virtual size_t element_size() const {
		return OpenVolumeMeshBaseProperty::UnknownSize;
	}
	virtual size_t size_of() const {
		return sizeof(data_);
	}

	virtual size_t size_of(size_t /* _n_elem */) const {
		return OpenVolumeMeshBaseProperty::UnknownSize;
	}

	virtual void stats(std::ostream& _ostr) const {
		for(vector_type::const_iterator it = data_.begin();
			it != data_.end(); ++it) {
				_ostr << *it << " ";
		}
	}

public:

	const value_type* data() const {
		if (data_.empty())
			return 0;

		return (value_type*) &data_[0];
	}

	/// Access the i'th element. No range check is performed!
	reference operator[](int _idx) {
		assert(size_t(_idx) < data_.size());
		return ((value_type*) &data_[0])[_idx];
	}

	/// Const access the i'th element. No range check is performed!
	const_reference operator[](int _idx) const {
		assert(size_t(_idx) < data_.size());
		return ((value_type*) &data_[0])[_idx];
	}

	OpenVolumeMeshPropertyT<value_type>* clone() const {
		OpenVolumeMeshPropertyT<value_type>* p = new OpenVolumeMeshPropertyT<
				value_type> (*this);
		return p;
	}

private:

	vector_type data_;

};

/// Base property handle.
template<class T>
struct OpenVolumeMeshBasePropHandleT: public OpenVolumeMeshHandle {

	typedef T 										Value;
	typedef std::vector<T> 							vector_type;
	typedef T 										value_type;
	typedef typename vector_type::reference 		reference;
	typedef typename vector_type::const_reference 	const_reference;

	explicit OpenVolumeMeshBasePropHandleT(int _idx = -1) :
		OpenVolumeMeshHandle(_idx) {
	}
};

/** \ingroup mesh_property_handle_group
 *  Handle representing a vertex property
 */
template<class T>
struct VPropHandleT: public OpenVolumeMeshBasePropHandleT<T> {
	typedef T Value;
	typedef T value_type;

	explicit VPropHandleT(int _idx = -1) :
		OpenVolumeMeshBasePropHandleT<T> (_idx) {
	}
	explicit VPropHandleT(const OpenVolumeMeshBasePropHandleT<T>& _b) :
		OpenVolumeMeshBasePropHandleT<T> (_b) {
	}
};

/** \ingroup mesh_property_handle_group
 *  Handle representing a halfedge property
 */
template<class T>
struct HEPropHandleT: public OpenVolumeMeshBasePropHandleT<T> {
	typedef T Value;
	typedef T value_type;

	explicit HEPropHandleT(int _idx = -1) :
		OpenVolumeMeshBasePropHandleT<T> (_idx) {
	}
	explicit HEPropHandleT(const OpenVolumeMeshBasePropHandleT<T>& _b) :
		OpenVolumeMeshBasePropHandleT<T> (_b) {
	}
};

/** \ingroup mesh_property_handle_group
 *  Handle representing an edge property
 */
template<class T>
struct EPropHandleT: public OpenVolumeMeshBasePropHandleT<T> {
	typedef T Value;
	typedef T value_type;

	explicit EPropHandleT(int _idx = -1) :
		OpenVolumeMeshBasePropHandleT<T> (_idx) {
	}
	explicit EPropHandleT(const OpenVolumeMeshBasePropHandleT<T>& _b) :
		OpenVolumeMeshBasePropHandleT<T> (_b) {
	}
};

/** \ingroup mesh_property_handle_group
 *  Handle representing a halfface property
 */
template<class T>
struct HFPropHandleT: public OpenVolumeMeshBasePropHandleT<T> {
	typedef T Value;
	typedef T value_type;

	explicit HFPropHandleT(int _idx = -1) :
		OpenVolumeMeshBasePropHandleT<T> (_idx) {
	}
	explicit HFPropHandleT(const OpenVolumeMeshBasePropHandleT<T>& _b) :
		OpenVolumeMeshBasePropHandleT<T> (_b) {
	}
};

/** \ingroup mesh_property_handle_group
 *  Handle representing a face property
 */
template<class T>
struct FPropHandleT: public OpenVolumeMeshBasePropHandleT<T> {
	typedef T Value;
	typedef T value_type;

	explicit FPropHandleT(int _idx = -1) :
		OpenVolumeMeshBasePropHandleT<T> (_idx) {
	}
	explicit FPropHandleT(const OpenVolumeMeshBasePropHandleT<T>& _b) :
		OpenVolumeMeshBasePropHandleT<T> (_b) {
	}
};

/** \ingroup mesh_property_handle_group
 *  Handle representing a cell property
 */
template<class T>
struct CPropHandleT: public OpenVolumeMeshBasePropHandleT<T> {
	typedef T Value;
	typedef T value_type;

	explicit CPropHandleT(int _idx = -1) :
		OpenVolumeMeshBasePropHandleT<T> (_idx) {
	}
	explicit CPropHandleT(const OpenVolumeMeshBasePropHandleT<T>& _b) :
		OpenVolumeMeshBasePropHandleT<T> (_b) {
	}
};

/** \ingroup mesh_property_handle_group
 *  Handle representing a mesh property
 */
template<class T>
struct MPropHandleT: public OpenVolumeMeshBasePropHandleT<T> {
	typedef T Value;
	typedef T value_type;

	explicit MPropHandleT(int _idx = -1) :
		OpenVolumeMeshBasePropHandleT<T> (_idx) {
	}
	explicit MPropHandleT(const OpenVolumeMeshBasePropHandleT<T>& _b) :
		OpenVolumeMeshBasePropHandleT<T> (_b) {
	}
};

} // Namespace OpenVolumeMesh

//=============================================================================
#endif // OPENVOLUMEMESHPROPERTY_HH defined
//=============================================================================
