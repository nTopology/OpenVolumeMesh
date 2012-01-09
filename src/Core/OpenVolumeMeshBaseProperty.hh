#ifndef OPENVOLUMEMESHBASEPROPERTY_HH
#define OPENVOLUMEMESHBASEPROPERTY_HH

#include <string>

//== CLASS DEFINITION =========================================================

/** \class OpenVolumeMeshBaseProperty

 Abstract class defining the basic interface of a dynamic property.

 **/

class OpenVolumeMeshBaseProperty {
public:

	/// Indicates an error when a size is returned by a member.
	static const size_t UnknownSize = size_t(-1);

public:

	OpenVolumeMeshBaseProperty(const std::string& _name = "<unknown>") :
		name_(_name), persistent_(false) {
	}

	OpenVolumeMeshBaseProperty(const OpenVolumeMeshBaseProperty& _rhs) :
		name_(_rhs.name_), persistent_(_rhs.persistent_) {
	}

	virtual ~OpenVolumeMeshBaseProperty() {}

public:

	/// Reserve memory for n elements.
	virtual void reserve(size_t _n) = 0;

	/// Resize storage to hold n elements.
	virtual void resize(size_t _n) = 0;

	/// Clear all elements and free memory.
	virtual void clear() = 0;

	/// Extend the number of elements by one.
	virtual void push_back() = 0;

	/// Let two elements swap their storage place.
	virtual void swap(size_t _i0, size_t _i1) = 0;

	/// Erase an element of the vector
	virtual void delete_element(size_t _idx) = 0;

	/// Return a deep copy of self.
	virtual OpenVolumeMeshBaseProperty* clone() const = 0;

public:

	/// Return the name of the property
	const std::string& name() const {
		return name_;
	}

	virtual void stats(std::ostream& _ostr) const {

		_ostr << name() << (persistent() ? ", persistent " : "") << "\n";
	}

public:
	// I/O support

	/// Returns true if the persistent flag is enabled else false.
	bool persistent() const {
		return persistent_;
	}

	/// Enable or disable persistency. Self must be a named property to enable
	/// persistency.
	void set_persistent(bool _yn) {
		persistent_ = _yn;
	}

	/// Number of elements in property
	virtual size_t n_elements() const = 0;

	/// Size of one element in bytes or UnknownSize if not known.
	virtual size_t element_size() const = 0;

	/// Return size of property in bytes
	virtual size_t size_of() const {
		return size_of(n_elements());
	}

	/// Estimated size of property if it has _n_elem elements.
	/// The member returns UnknownSize if the size cannot be estimated.
	virtual size_t size_of(size_t _n_elem) const {
		return (element_size() != UnknownSize) ? (_n_elem * element_size())
				: UnknownSize;
	}

private:

	std::string name_;
	bool persistent_;
};

#endif //OPENVOLUMEMESHBASEPROPERTY_HH

