//=============================================================================
//
//  CLASS OpenVolumeMeshBaseKernel
//
//=============================================================================


#ifndef OPENVOLUMEMESHBASEKERNEL_HH
#define OPENVOLUMEMESHBASEKERNEL_HH

//== INCLUDES =================================================================

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
// --------------------
#include "OpenVolumeMeshProperty.hh"
#include "OpenVolumeMeshPropertyContainer.hh"
#include "OpenVolumeMeshHandle.hh"

//== CLASS DEFINITION =========================================================

/// This class provides the basic property management like adding/removing
/// properties and access to properties.
/// All operations provided by %OpenVolumeMeshBaseKernel need at least a property handle
/// (VPropHandleT, HEPropHandleT, EPropHandleT, HFPropHandleT, FPropHandleT, CPropHandleT, MPropHandleT).
/// which keeps the data type of the property, too.

class OpenVolumeMeshBaseKernel {
public:
	//-------------------------------------------- constructor / destructor

	OpenVolumeMeshBaseKernel() {
	}

	virtual ~OpenVolumeMeshBaseKernel() {
		vprops_.clear();
		eprops_.clear();
		heprops_.clear();
		fprops_.clear();
		hfprops_.clear();
		cprops_.clear();
		mprops_.clear();
	}

public:
	//-------------------------------------------------- add new properties

	template<class T>
	void add_property(VPropHandleT<T>& _ph,
			const std::string& _name = "<vprop>") {
		_ph = VPropHandleT<T> (vprops_.add(T(), _name));
		vprops_.resize(n_vertices());
	}

	template<class T>
	void add_property(HEPropHandleT<T>& _ph,
			const std::string& _name = "<heprop>") {
		_ph = HEPropHandleT<T> (heprops_.add(T(), _name));
		heprops_.resize(n_halfedges());
	}

	template<class T>
	void add_property(EPropHandleT<T>& _ph,
			const std::string& _name = "<eprop>") {
		_ph = EPropHandleT<T> (eprops_.add(T(), _name));
		eprops_.resize(n_edges());
	}

	template<class T>
	void add_property(FPropHandleT<T>& _ph,
			const std::string& _name = "<fprop>") {
		_ph = FPropHandleT<T> (fprops_.add(T(), _name));
		fprops_.resize(n_faces());
	}

	template<class T>
	void add_property(HFPropHandleT<T>& _ph,
			const std::string& _name = "<hfprop>") {
		_ph = HFPropHandleT<T> (hfprops_.add(T(), _name));
		hfprops_.resize(n_halffaces());
	}

	template<class T>
	void add_property(CPropHandleT<T>& _ph,
			const std::string& _name = "<cprop>") {
		_ph = CPropHandleT<T> (cprops_.add(T(), _name));
		cprops_.resize(n_cells());
	}

	template<class T>
	void add_property(MPropHandleT<T>& _ph,
			const std::string& _name = "<mprop>") {
		_ph = MPropHandleT<T> (mprops_.add(T(), _name));
		mprops_.resize(1);
	}


public:
	//--------------------------------------------------- remove properties

	template<typename T>
	void remove_property(VPropHandleT<T>& _ph) {
		if (_ph.is_valid())
			vprops_.remove(_ph);
		_ph.reset();
	}

	template<typename T>
	void remove_property(HEPropHandleT<T>& _ph) {
		if (_ph.is_valid())
			heprops_.remove(_ph);
		_ph.reset();
	}

	template<typename T>
	void remove_property(EPropHandleT<T>& _ph) {
		if (_ph.is_valid())
			eprops_.remove(_ph);
		_ph.reset();
	}

	template<typename T>
	void remove_property(FPropHandleT<T>& _ph) {
		if (_ph.is_valid())
			fprops_.remove(_ph);
		_ph.reset();
	}

	template<typename T>
	void remove_property(HFPropHandleT<T>& _ph) {
		if (_ph.is_valid())
			hfprops_.remove(_ph);
		_ph.reset();
	}

	template<typename T>
	void remove_property(CPropHandleT<T>& _ph) {
		if (_ph.is_valid())
			cprops_.remove(_ph);
		_ph.reset();
	}

	template<typename T>
	void remove_property(MPropHandleT<T>& _ph) {
		if (_ph.is_valid())
			mprops_.remove(_ph);
		_ph.reset();
	}

	//--------------------------------------------------- remove single elements

	void remove_vprop_element(size_t _idx) {
	    for(prop_iterator p_it = vprops_.begin(); p_it != vprops_.end(); ++p_it) {
	        (*p_it)->delete_element(_idx);
	    }
	}

    void remove_eprop_element(size_t _idx) {
        for(prop_iterator p_it = eprops_.begin(); p_it != eprops_.end(); ++p_it) {
            (*p_it)->delete_element(_idx);
        }
    }

    void remove_heprop_element(size_t _idx) {
        for(prop_iterator p_it = heprops_.begin(); p_it != heprops_.end(); ++p_it) {
            (*p_it)->delete_element(_idx);
        }
    }

    void remove_fprop_element(size_t _idx) {
        for(prop_iterator p_it = fprops_.begin(); p_it != fprops_.end(); ++p_it) {
            (*p_it)->delete_element(_idx);
        }
    }

    void remove_hfprop_element(size_t _idx) {
        for(prop_iterator p_it = hfprops_.begin(); p_it != hfprops_.end(); ++p_it) {
            (*p_it)->delete_element(_idx);
        }
    }

    void remove_cprop_element(size_t _idx) {
        for(prop_iterator p_it = cprops_.begin(); p_it != cprops_.end(); ++p_it) {
            (*p_it)->delete_element(_idx);
        }
    }

public:
	//------------------------------------------------ get handle from name

	template<class T>
	bool get_property_handle(VPropHandleT<T>& _ph, const std::string& _name) const {
		return (_ph = VPropHandleT<T> (vprops_.handle(T(), _name))).is_valid();
	}

	template<class T>
	bool get_property_handle(HEPropHandleT<T>& _ph, const std::string& _name) const {
		return (_ph = HEPropHandleT<T> (heprops_.handle(T(), _name))).is_valid();
	}

	template<class T>
	bool get_property_handle(EPropHandleT<T>& _ph, const std::string& _name) const {
		return (_ph = EPropHandleT<T> (eprops_.handle(T(), _name))).is_valid();
	}

	template<class T>
	bool get_property_handle(FPropHandleT<T>& _ph, const std::string& _name) const {
		return (_ph = FPropHandleT<T> (fprops_.handle(T(), _name))).is_valid();
	}

	template<class T>
	bool get_property_handle(HFPropHandleT<T>& _ph, const std::string& _name) const {
		return (_ph = HFPropHandleT<T> (hfprops_.handle(T(), _name))).is_valid();
	}

	template<class T>
	bool get_property_handle(CPropHandleT<T>& _ph, const std::string& _name) const {
		return (_ph = CPropHandleT<T> (cprops_.handle(T(), _name))).is_valid();
	}

	template<class T>
	bool get_property_handle(MPropHandleT<T>& _ph, const std::string& _name) const {
		return (_ph = MPropHandleT<T> (mprops_.handle(T(), _name))).is_valid();
	}


public:
	//--------------------------------------------------- access properties

	template<class T>
	OpenVolumeMeshPropertyT<T>& property(VPropHandleT<T> _ph) {
		return vprops_.property(_ph);
	}
	template<class T>
	const OpenVolumeMeshPropertyT<T>& property(VPropHandleT<T> _ph) const {
		return vprops_.property(_ph);
	}

	template<class T>
	OpenVolumeMeshPropertyT<T>& property(HEPropHandleT<T> _ph) {
		return heprops_.property(_ph);
	}
	template<class T>
	const OpenVolumeMeshPropertyT<T>& property(HEPropHandleT<T> _ph) const {
		return heprops_.property(_ph);
	}

	template<class T>
	OpenVolumeMeshPropertyT<T>& property(EPropHandleT<T> _ph) {
		return eprops_.property(_ph);
	}
	template<class T>
	const OpenVolumeMeshPropertyT<T>& property(EPropHandleT<T> _ph) const {
		return eprops_.property(_ph);
	}

	template<class T>
	OpenVolumeMeshPropertyT<T>& property(FPropHandleT<T> _ph) {
		return fprops_.property(_ph);
	}
	template<class T>
	const OpenVolumeMeshPropertyT<T>& property(FPropHandleT<T> _ph) const {
		return fprops_.property(_ph);
	}

	template<class T>
	OpenVolumeMeshPropertyT<T>& property(HFPropHandleT<T> _ph) {
		return hfprops_.property(_ph);
	}
	template<class T>
	const OpenVolumeMeshPropertyT<T>& property(HFPropHandleT<T> _ph) const {
		return hfprops_.property(_ph);
	}

	template<class T>
	OpenVolumeMeshPropertyT<T>& property(CPropHandleT<T> _ph) {
		return cprops_.property(_ph);
	}
	template<class T>
	const OpenVolumeMeshPropertyT<T>& property(CPropHandleT<T> _ph) const {
		return cprops_.property(_ph);
	}

	template<class T>
	OpenVolumeMeshPropertyT<T>& mproperty(MPropHandleT<T> _ph) {
		return mprops_.property(_ph);
	}
	template<class T>
	const OpenVolumeMeshPropertyT<T>& mproperty(MPropHandleT<T> _ph) const {
		return mprops_.property(_ph);
	}


public:
	//-------------------------------------------- access property elements

	template<class T>
	typename VPropHandleT<T>::reference property(VPropHandleT<T> _ph,
			VertexHandle _vh) {
		return vprops_.property(_ph)[_vh.idx()];
	}

	template<class T>
	typename VPropHandleT<T>::const_reference property(VPropHandleT<T> _ph,
			VertexHandle _vh) const {
		return vprops_.property(_ph)[_vh.idx()];
	}

	template<class T>
	typename HEPropHandleT<T>::reference property(HEPropHandleT<T> _ph,
			HalfEdgeHandle _hh) {
		return heprops_.property(_ph)[_hh.idx()];
	}

	template<class T>
	typename HEPropHandleT<T>::const_reference property(HEPropHandleT<T> _ph,
			HalfEdgeHandle _hh) const {
		return heprops_.property(_ph)[_hh.idx()];
	}

	template<class T>
	typename EPropHandleT<T>::reference property(EPropHandleT<T> _ph,
			EdgeHandle _eh) {
		return eprops_.property(_ph)[_eh.idx()];
	}

	template<class T>
	typename EPropHandleT<T>::const_reference property(EPropHandleT<T> _ph,
			EdgeHandle _eh) const {
		return eprops_.property(_ph)[_eh.idx()];
	}

	template<class T>
	typename FPropHandleT<T>::reference property(FPropHandleT<T> _ph,
			FaceHandle _fh) {
		return fprops_.property(_ph)[_fh.idx()];
	}

	template<class T>
	typename FPropHandleT<T>::const_reference property(FPropHandleT<T> _ph,
			FaceHandle _fh) const {
		return fprops_.property(_ph)[_fh.idx()];
	}

	template<class T>
	typename HFPropHandleT<T>::reference property(HFPropHandleT<T> _ph,
			HalfFaceHandle _hfh) {
		return hfprops_.property(_ph)[_hfh.idx()];
	}

	template<class T>
	typename HFPropHandleT<T>::const_reference property(HFPropHandleT<T> _ph,
			HalfFaceHandle _hfh) const {
		return hfprops_.property(_ph)[_hfh.idx()];
	}

	template<class T>
	typename CPropHandleT<T>::reference property(CPropHandleT<T> _ph,
			CellHandle _ch) {
		return cprops_.property(_ph)[_ch.idx()];
	}

	template<class T>
	typename CPropHandleT<T>::const_reference property(CPropHandleT<T> _ph,
			CellHandle _ch) const {
		return cprops_.property(_ph)[_ch.idx()];
	}

	template<class T>
	typename MPropHandleT<T>::reference property(MPropHandleT<T> _ph) {
		return mprops_.property(_ph)[0];
	}

	template<class T>
	typename MPropHandleT<T>::const_reference property(MPropHandleT<T> _ph) const {
		return mprops_.property(_ph)[0];
	}

	//@}

protected:
	//------------------------------------------------- low-level access

public:
	// used by non-native kernel and MeshIO, should be protected

	size_t n_vprops(void) const {
		return vprops_.size();
	}

	size_t n_eprops(void) const {
		return eprops_.size();
	}

	size_t n_heprops(void) const {
		return heprops_.size();
	}

	size_t n_fprops(void) const {
		return fprops_.size();
	}

	size_t n_hfprops(void) const {
		return hfprops_.size();
	}

	size_t n_cprops(void) const {
		return cprops_.size();
	}

	size_t n_mprops(void) const {
		return mprops_.size();
	}

	OpenVolumeMeshBaseProperty* _get_vprop(const std::string& _name) {
		return vprops_.property(_name);
	}

	OpenVolumeMeshBaseProperty* _get_eprop(const std::string& _name) {
		return eprops_.property(_name);
	}

	OpenVolumeMeshBaseProperty* _get_heprop(const std::string& _name) {
		return heprops_.property(_name);
	}

	OpenVolumeMeshBaseProperty* _get_fprop(const std::string& _name) {
		return fprops_.property(_name);
	}

	OpenVolumeMeshBaseProperty* _get_hfprop(const std::string& _name) {
		return hfprops_.property(_name);
	}

	OpenVolumeMeshBaseProperty* _get_cprop(const std::string& _name) {
		return cprops_.property(_name);
	}

	OpenVolumeMeshBaseProperty* _get_mprop(const std::string& _name) {
		return mprops_.property(_name);
	}

	const OpenVolumeMeshBaseProperty* _get_vprop(const std::string& _name) const {
		return vprops_.property(_name);
	}

	const OpenVolumeMeshBaseProperty* _get_eprop(const std::string& _name) const {
		return eprops_.property(_name);
	}

	const OpenVolumeMeshBaseProperty* _get_heprop(const std::string& _name) const {
		return heprops_.property(_name);
	}

	const OpenVolumeMeshBaseProperty* _get_fprop(const std::string& _name) const {
		return fprops_.property(_name);
	}

	const OpenVolumeMeshBaseProperty* _get_hfprop(const std::string& _name) const {
		return hfprops_.property(_name);
	}

	const OpenVolumeMeshBaseProperty* _get_cprop(const std::string& _name) const {
		return cprops_.property(_name);
	}

	const OpenVolumeMeshBaseProperty* _get_mprop(const std::string& _name) const {
		return mprops_.property(_name);
	}

	OpenVolumeMeshBaseProperty& _vprop(size_t _idx) {
		return vprops_._property(_idx);
	}
	OpenVolumeMeshBaseProperty& _eprop(size_t _idx) {
		return eprops_._property(_idx);
	}
	OpenVolumeMeshBaseProperty& _heprop(size_t _idx) {
		return heprops_._property(_idx);
	}
	OpenVolumeMeshBaseProperty& _fprop(size_t _idx) {
		return fprops_._property(_idx);
	}
	OpenVolumeMeshBaseProperty& _hfprop(size_t _idx) {
		return hfprops_._property(_idx);
	}
	OpenVolumeMeshBaseProperty& _cprop(size_t _idx) {
		return cprops_._property(_idx);
	}
	OpenVolumeMeshBaseProperty& _mprop(size_t _idx) {
		return mprops_._property(_idx);
	}

	const OpenVolumeMeshBaseProperty& _vprop(size_t _idx) const {
		return vprops_._property(_idx);
	}
	const OpenVolumeMeshBaseProperty& _eprop(size_t _idx) const {
		return eprops_._property(_idx);
	}
	const OpenVolumeMeshBaseProperty& _heprop(size_t _idx) const {
		return heprops_._property(_idx);
	}
	const OpenVolumeMeshBaseProperty& _fprop(size_t _idx) const {
		return fprops_._property(_idx);
	}
	const OpenVolumeMeshBaseProperty& _hfprop(size_t _idx) const {
		return hfprops_._property(_idx);
	}
	const OpenVolumeMeshBaseProperty& _cprop(size_t _idx) const {
		return cprops_._property(_idx);
	}
	const OpenVolumeMeshBaseProperty& _mprop(size_t _idx) const {
		return mprops_._property(_idx);
	}

	size_t _add_vprop(OpenVolumeMeshBaseProperty* _bp) {
		return vprops_._add(_bp);
	}
	size_t _add_eprop(OpenVolumeMeshBaseProperty* _bp) {
		return eprops_._add(_bp);
	}
	size_t _add_heprop(OpenVolumeMeshBaseProperty* _bp) {
		return heprops_._add(_bp);
	}
	size_t _add_fprop(OpenVolumeMeshBaseProperty* _bp) {
		return fprops_._add(_bp);
	}
	size_t _add_hfprop(OpenVolumeMeshBaseProperty* _bp) {
		return hfprops_._add(_bp);
	}
	size_t _add_cprop(OpenVolumeMeshBaseProperty* _bp) {
		return cprops_._add(_bp);
	}
	size_t _add_mprop(OpenVolumeMeshBaseProperty* _bp) {
		return mprops_._add(_bp);
	}

protected:
	// low-level access non-public

	OpenVolumeMeshBaseProperty& _vprop(OpenVolumeMeshHandle _h) {
		return vprops_._property(_h.idx());
	}
	OpenVolumeMeshBaseProperty& _eprop(OpenVolumeMeshHandle _h) {
		return eprops_._property(_h.idx());
	}
	OpenVolumeMeshBaseProperty& _heprop(OpenVolumeMeshHandle _h) {
		return heprops_._property(_h.idx());
	}
	OpenVolumeMeshBaseProperty& _fprop(OpenVolumeMeshHandle _h) {
		return fprops_._property(_h.idx());
	}
	OpenVolumeMeshBaseProperty& _hfprop(OpenVolumeMeshHandle _h) {
		return hfprops_._property(_h.idx());
	}
	OpenVolumeMeshBaseProperty& _cprop(OpenVolumeMeshHandle _h) {
		return cprops_._property(_h.idx());
	}
	OpenVolumeMeshBaseProperty& _mprop(OpenVolumeMeshHandle _h) {
		return mprops_._property(_h.idx());
	}

	const OpenVolumeMeshBaseProperty& _vprop(OpenVolumeMeshHandle _h) const {
		return vprops_._property(_h.idx());
	}
	const OpenVolumeMeshBaseProperty& _eprop(OpenVolumeMeshHandle _h) const {
		return eprops_._property(_h.idx());
	}
	const OpenVolumeMeshBaseProperty& _heprop(OpenVolumeMeshHandle _h) const {
		return heprops_._property(_h.idx());
	}
	const OpenVolumeMeshBaseProperty& _fprop(OpenVolumeMeshHandle _h) const {
		return fprops_._property(_h.idx());
	}
	const OpenVolumeMeshBaseProperty& _hfprop(OpenVolumeMeshHandle _h) const {
		return hfprops_._property(_h.idx());
	}
	const OpenVolumeMeshBaseProperty& _cprop(OpenVolumeMeshHandle _h) const {
		return cprops_._property(_h.idx());
	}
	const OpenVolumeMeshBaseProperty& _mprop(OpenVolumeMeshHandle _h) const {
		return mprops_._property(_h.idx());
	}

public:
	//----------------------------------------------------- element numbers


	virtual unsigned int n_vertices() const {
		return 0;
	}
	virtual unsigned int n_halfedges() const {
		return 0;
	}
	virtual unsigned int n_edges() const {
		return 0;
	}
	virtual unsigned int n_faces() const {
		return 0;
	}
	virtual unsigned int n_halffaces() const {
		return 0;
	}
	virtual unsigned int n_cells() const {
		return 0;
	}

protected:
	//------------------------------------------- synchronize properties

	void vprops_reserve(unsigned int _n) const {
		vprops_.reserve(_n);
	}
	void vprops_resize(unsigned int _n) const {
		vprops_.resize(_n);
	}
	void vprops_clear() {
		vprops_.clear();
	}
	void vprops_swap(unsigned int _i0, unsigned int _i1) const {
		vprops_.swap(_i0, _i1);
	}

	void heprops_reserve(unsigned int _n) const {
		heprops_.reserve(_n);
	}
	void heprops_resize(unsigned int _n) const {
		heprops_.resize(_n);
	}
	void heprops_clear() {
		heprops_.clear();
	}
	void heprops_swap(unsigned int _i0, unsigned int _i1) const {
		heprops_.swap(_i0, _i1);
	}

	void eprops_reserve(unsigned int _n) const {
		eprops_.reserve(_n);
	}
	void eprops_resize(unsigned int _n) const {
		eprops_.resize(_n);
	}
	void eprops_clear() {
		eprops_.clear();
	}
	void eprops_swap(unsigned int _i0, unsigned int _i1) const {
		eprops_.swap(_i0, _i1);
	}

	void fprops_reserve(unsigned int _n) const {
		fprops_.reserve(_n);
	}
	void fprops_resize(unsigned int _n) const {
		fprops_.resize(_n);
	}
	void fprops_clear() {
		fprops_.clear();
	}
	void fprops_swap(unsigned int _i0, unsigned int _i1) const {
		fprops_.swap(_i0, _i1);
	}

	void hfprops_reserve(unsigned int _n) const {
		hfprops_.reserve(_n);
	}
	void hfprops_resize(unsigned int _n) const {
		hfprops_.resize(_n);
	}
	void hfprops_clear() {
		hfprops_.clear();
	}
	void hfprops_swap(unsigned int _i0, unsigned int _i1) const {
		hfprops_.swap(_i0, _i1);
	}

	void cprops_reserve(unsigned int _n) const {
		cprops_.reserve(_n);
	}
	void cprops_resize(unsigned int _n) const {
		cprops_.resize(_n);
	}
	void cprops_clear() {
		cprops_.clear();
	}
	void cprops_swap(unsigned int _i0, unsigned int _i1) const {
		cprops_.swap(_i0, _i1);
	}

	void mprops_resize(unsigned int _n) const {
		mprops_.resize(_n);
	}
	void mprops_clear() {
		mprops_.clear();
	}

public:

	void property_stats(std::ostream& _ostr = std::clog) const;

	void vprop_stats(std::string& _string) const;
	void heprop_stats(std::string& _string) const;
	void eprop_stats(std::string& _string) const;
	void fprop_stats(std::string& _string) const;
	void hfprop_stats(std::string& _string) const;
	void cprop_stats(std::string& _string) const;
	void mprop_stats(std::string& _string) const;

	void vprop_stats(std::ostream& _ostr = std::clog) const;
	void heprop_stats(std::ostream& _ostr = std::clog) const;
	void eprop_stats(std::ostream& _ostr = std::clog) const;
	void fprop_stats(std::ostream& _ostr = std::clog) const;
	void hfprop_stats(std::ostream& _ostr = std::clog) const;
	void cprop_stats(std::ostream& _ostr = std::clog) const;
	void mprop_stats(std::ostream& _ostr = std::clog) const;

public:

	typedef OpenVolumeMeshPropertyContainer::iterator prop_iterator;
	typedef OpenVolumeMeshPropertyContainer::const_iterator const_prop_iterator;

	prop_iterator vprops_begin() {
		return vprops_.begin();
	}
	prop_iterator vprops_end() {
		return vprops_.end();
	}
	const_prop_iterator vprops_begin() const {
		return vprops_.begin();
	}
	const_prop_iterator vprops_end() const {
		return vprops_.end();
	}

	prop_iterator eprops_begin() {
		return eprops_.begin();
	}
	prop_iterator eprops_end() {
		return eprops_.end();
	}
	const_prop_iterator eprops_begin() const {
		return eprops_.begin();
	}
	const_prop_iterator eprops_end() const {
		return eprops_.end();
	}

	prop_iterator heprops_begin() {
		return heprops_.begin();
	}
	prop_iterator heprops_end() {
		return heprops_.end();
	}
	const_prop_iterator heprops_begin() const {
		return heprops_.begin();
	}
	const_prop_iterator heprops_end() const {
		return heprops_.end();
	}

	prop_iterator fprops_begin() {
		return fprops_.begin();
	}
	prop_iterator fprops_end() {
		return fprops_.end();
	}
	const_prop_iterator fprops_begin() const {
		return fprops_.begin();
	}
	const_prop_iterator fprops_end() const {
		return fprops_.end();
	}

	prop_iterator hfprops_begin() {
		return hfprops_.begin();
	}
	prop_iterator hfprops_end() {
		return hfprops_.end();
	}
	const_prop_iterator hfprops_begin() const {
		return hfprops_.begin();
	}
	const_prop_iterator hfprops_end() const {
		return hfprops_.end();
	}

	prop_iterator cprops_begin() {
		return cprops_.begin();
	}
	prop_iterator cprops_end() {
		return cprops_.end();
	}
	const_prop_iterator cprops_begin() const {
		return cprops_.begin();
	}
	const_prop_iterator cprops_end() const {
		return cprops_.end();
	}

	prop_iterator mprops_begin() {
		return mprops_.begin();
	}
	prop_iterator mprops_end() {
		return mprops_.end();
	}
	const_prop_iterator mprops_begin() const {
		return mprops_.begin();
	}
	const_prop_iterator mprops_end() const {
		return mprops_.end();
	}

private:

	OpenVolumeMeshPropertyContainer vprops_;
	OpenVolumeMeshPropertyContainer eprops_;
	OpenVolumeMeshPropertyContainer heprops_;
	OpenVolumeMeshPropertyContainer fprops_;
	OpenVolumeMeshPropertyContainer hfprops_;
	OpenVolumeMeshPropertyContainer cprops_;
	OpenVolumeMeshPropertyContainer mprops_;
};

//=============================================================================
#endif // OPENVOLUMEMESHBASEKERNEL_HH defined
//=============================================================================
