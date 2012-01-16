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

#ifndef OPENVOLUMEMESHPROPERTYCONTAINER_HH
#define OPENVOLUMEMESHPROPERTYCONTAINER_HH

#include <algorithm>
#include <cassert>

#include "OpenVolumeMeshProperty.hh"

namespace OpenVolumeMesh {

//== FORWARDDECLARATIONS ======================================================
class BaseKernel;

//== CLASS DEFINITION =========================================================
/// A a container for properties.
class OpenVolumeMeshPropertyContainer {
public:

	//-------------------------------------------------- constructor / destructor

	OpenVolumeMeshPropertyContainer() {
	}

	virtual ~OpenVolumeMeshPropertyContainer() {
		std::for_each(properties_.begin(), properties_.end(), Delete());
	}

	//------------------------------------------------------------- info / access

	typedef std::vector<OpenVolumeMeshBaseProperty*> Properties;
	const Properties& properties() const {
		return properties_;
	}
	size_t size() const {
		return properties_.size();
	}

	//--------------------------------------------------------- copy / assignment

	OpenVolumeMeshPropertyContainer(const OpenVolumeMeshPropertyContainer& _rhs) {
		operator=(_rhs);
	}

	OpenVolumeMeshPropertyContainer& operator=(const OpenVolumeMeshPropertyContainer& _rhs) {
		// The assignment below relies on all previous OpenVolumeMeshBaseProperty* elements having been deleted
		std::for_each(properties_.begin(), properties_.end(), Delete());
		properties_ = _rhs.properties_;
		Properties::iterator p_it = properties_.begin(), p_end =
				properties_.end();
		for (; p_it != p_end; ++p_it)
			if (*p_it)
				*p_it = (*p_it)->clone();
		return *this;
	}

	//--------------------------------------------------------- manage properties

	template<class T>
	OpenVolumeMeshBasePropHandleT<T> add(const T&, const std::string& _name = "<unknown>") {
		Properties::iterator p_it = properties_.begin(), p_end =
				properties_.end();
		int idx = 0;
		for (; p_it != p_end && *p_it != NULL; ++p_it, ++idx) {
		}
		if (p_it == p_end)
			properties_.push_back(NULL);
		properties_[idx] = new OpenVolumeMeshPropertyT<T> (_name);
		return OpenVolumeMeshBasePropHandleT<T> (idx);
	}

	template<class T>
	OpenVolumeMeshBasePropHandleT<T> handle(const T&, const std::string& _name) const {
		Properties::const_iterator p_it = properties_.begin();
		for (int idx = 0; p_it != properties_.end(); ++p_it, ++idx) {
			if (*p_it != NULL && (*p_it)->name() == _name) {
				return OpenVolumeMeshBasePropHandleT<T> (idx);
			}
		}
		return OpenVolumeMeshBasePropHandleT<T> ();
	}

	OpenVolumeMeshBaseProperty* property(const std::string& _name) const {
		Properties::const_iterator p_it = properties_.begin();
		for (int idx = 0; p_it != properties_.end(); ++p_it, ++idx) {
			if (*p_it != NULL && (*p_it)->name() == _name) //skip deleted properties
			{
				return *p_it;
			}
		}
		return NULL;
	}

	template<class T> OpenVolumeMeshPropertyT<T>& property(OpenVolumeMeshBasePropHandleT<T> _h) {
		assert(_h.idx() >= 0 && _h.idx() < (int) properties_.size());
		assert(properties_[_h.idx()] != NULL);
		OpenVolumeMeshPropertyT<T>* p = dynamic_cast<OpenVolumeMeshPropertyT<T>*> (properties_[_h.idx()]);
		assert(p != NULL);
		return *p;
	}

	template<class T> const OpenVolumeMeshPropertyT<T>& property(OpenVolumeMeshBasePropHandleT<T> _h) const {
		assert(_h.idx() >= 0 && _h.idx() < (int) properties_.size());
		assert(properties_[_h.idx()] != NULL);
		OpenVolumeMeshPropertyT<T>* p = dynamic_cast<OpenVolumeMeshPropertyT<T>*> (properties_[_h.idx()]);
		assert(p != NULL);
		return *p;
	}

	template<class T> void remove(OpenVolumeMeshBasePropHandleT<T> _h) {
		assert(_h.idx() >= 0 && _h.idx() < (int) properties_.size());
		delete properties_[_h.idx()];
		properties_[_h.idx()] = NULL;
	}

	void clear() {
		// Clear properties vector:
		// Replaced the old version with new one
		// which performs a swap to clear values and
		// deallocate memory.

		// Old version (changed 22.07.09) {
		// std::for_each(properties_.begin(), properties_.end(), Delete());
		// }

		std::for_each(properties_.begin(), properties_.end(), ClearAll());
		properties_.clear();
	}

	//---------------------------------------------------- synchronize properties

	void reserve(size_t _n) const {
		std::for_each(properties_.begin(), properties_.end(), Reserve(_n));
	}

	void resize(size_t _n) const {
		std::for_each(properties_.begin(), properties_.end(), Resize(_n));
	}

	void swap(size_t _i0, size_t _i1) const {
		std::for_each(properties_.begin(), properties_.end(), Swap(_i0, _i1));
	}

protected:
	// generic add/get

	size_t _add(OpenVolumeMeshBaseProperty* _bp) {
		Properties::iterator p_it = properties_.begin(), p_end =
				properties_.end();
		size_t idx = 0;
		for (; p_it != p_end && *p_it != NULL; ++p_it, ++idx) {
		};
		if (p_it == p_end)
			properties_.push_back(NULL);
		properties_[idx] = _bp;
		return idx;
	}

	OpenVolumeMeshBaseProperty& _property(size_t _idx) {
		assert(_idx < properties_.size());
		assert(properties_[_idx] != NULL);
		OpenVolumeMeshBaseProperty *p = properties_[_idx];
		assert(p != NULL);
		return *p;
	}

	const OpenVolumeMeshBaseProperty& _property(size_t _idx) const {
		assert(_idx < properties_.size());
		assert(properties_[_idx] != NULL);
		OpenVolumeMeshBaseProperty *p = properties_[_idx];
		assert(p != NULL);
		return *p;
	}

	typedef Properties::iterator iterator;
	typedef Properties::const_iterator const_iterator;
	iterator begin() {
		return properties_.begin();
	}
	iterator end() {
		return properties_.end();
	}
	const_iterator begin() const {
		return properties_.begin();
	}
	const_iterator end() const {
		return properties_.end();
	}

	friend class OpenVolumeMeshBaseKernel;

private:

	//-------------------------------------------------- synchronization functors

#ifndef DOXY_IGNORE_THIS
	struct Reserve {
		Reserve(size_t _n) :
			n_(_n) {
		}
		void operator()(OpenVolumeMeshBaseProperty* _p) const {
			if (_p)
				_p->reserve(n_);
		}
		size_t n_;
	};

	struct Resize {
		Resize(size_t _n) :
			n_(_n) {
		}
		void operator()(OpenVolumeMeshBaseProperty* _p) const {
			if (_p)
				_p->resize(n_);
		}
		size_t n_;
	};

	struct ClearAll {
		ClearAll() {
		}
		void operator()(OpenVolumeMeshBaseProperty* _p) const {
			if (_p)
				_p->clear();
		}
	};

	struct Swap {
		Swap(size_t _i0, size_t _i1) :
			i0_(_i0), i1_(_i1) {
		}
		void operator()(OpenVolumeMeshBaseProperty* _p) const {
			if (_p)
				_p->swap(i0_, i1_);
		}
		size_t i0_, i1_;
	};

	struct Delete {
		Delete() {
		}
		void operator()(OpenVolumeMeshBaseProperty* _p) const {
			if (_p)
				delete _p;
			_p = NULL;
		}
	};
#endif

	Properties properties_;
};

} // Namespace OpenVolumeMesh

#endif //OPENVOLUMEMESHPROPERTYCONTAINER_HH
