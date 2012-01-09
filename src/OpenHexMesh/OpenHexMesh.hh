/*
 * OpenHexMesh.hh
 *
 *  Created on: Jun 10, 2011
 *      Author: kremer
 */

#ifndef OPENHEXMESH_HH_
#define OPENHEXMESH_HH_

#include <set>

#include "../OpenVolumeMesh/OpenVolumeMesh.hh"
#include "Utils/OpenHexMeshIterators.hh"

/*
 * HexMesh data structure basing on OpenVolumeMesh
 *
 * Define hexahedra to span a 3D coordinate system:
 *
 * X-axis: From X-Front to X-Back face
 * Y-axis: From Y-Front to Y-Back face
 * Z-axis: From Z-Front to Z-Back face
 */

template <typename VecT>
class OpenHexMesh : public OpenVolumeMesh<VecT> {
public:

    typedef typename OpenVolumeMesh<VecT>::VertexHandle     VertexHandle;
    typedef typename OpenVolumeMesh<VecT>::EdgeHandle       EdgeHandle;
    typedef typename OpenVolumeMesh<VecT>::HalfEdgeHandle   HalfEdgeHandle;
    typedef typename OpenVolumeMesh<VecT>::FaceHandle       FaceHandle;
    typedef typename OpenVolumeMesh<VecT>::HalfFaceHandle   HalfFaceHandle;
    typedef typename OpenVolumeMesh<VecT>::CellHandle       CellHandle;

    // Orientation constants
    static const unsigned char XF = 0;
    static const unsigned char XB = 1;
    static const unsigned char YF = 2;
    static const unsigned char YB = 3;
    static const unsigned char ZF = 4;
    static const unsigned char ZB = 5;
    static const unsigned char INVALID = 6;

    static inline unsigned char opposite_orientation(const unsigned char _d) {
        return (_d % 2 == 0 ? _d + 1 : _d - 1);
    }

    // Constructor
    OpenHexMesh();

    // Destructor
    ~OpenHexMesh();

    // Overridden function
    FaceHandle add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck = true);

    // Overridden function
    FaceHandle add_face(const std::vector<VertexHandle>& _vertices);

    /// Overridden function
    CellHandle add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck = true,
            bool _reorderFaces = false);

    // ======================= Specialized Iterators =============================

    friend class CellSheetCellIter<VecT>;
    friend class HalfFaceSheetHalfFaceIter<VecT>;
    friend class OutsideNeighborHalfFaceIter<VecT>;

    typedef class CellSheetCellIter<VecT> 		    CellSheetCellIter;
    typedef class HalfFaceSheetHalfFaceIter<VecT>   HalfFaceSheetHalfFaceIter;
    typedef class OutsideNeighborHalfFaceIter<VecT> OutsideNeighborHalfFaceIter;

    CellSheetCellIter csc_iter(const CellHandle& _ref_h, const unsigned char _orthDir) const {
        return CellSheetCellIter(_ref_h, _orthDir, this);
    }

    HalfFaceSheetHalfFaceIter hfshf_iter(const HalfFaceHandle& _ref_h) const {
        return HalfFaceSheetHalfFaceIter(_ref_h, this);
    }

    OutsideNeighborHalfFaceIter onhf_iter(const HalfFaceHandle& _ref_h) const {
        return OutsideNeighborHalfFaceIter(_ref_h, this);
    }

    // ======================= Connectivity functions =============================

    inline HalfFaceHandle opposite_halfface_handle_in_cell(const HalfFaceHandle& _hfh, const CellHandle& _ch) {

        assert((unsigned int)_ch < OpenVolumeMesh<VecT>::cells_.size());

        if(orientation(_hfh, _ch) == XF) return xback_halfface(_ch);
        if(orientation(_hfh, _ch) == XB) return xfront_halfface(_ch);
        if(orientation(_hfh, _ch) == YF) return yback_halfface(_ch);
        if(orientation(_hfh, _ch) == YB) return yfront_halfface(_ch);
        if(orientation(_hfh, _ch) == ZF) return zback_halfface(_ch);
        if(orientation(_hfh, _ch) == ZB) return zfront_halfface(_ch);

        return OpenVolumeMesh<VecT>::InvalidHalfFaceHandle;
    }

    inline HalfFaceHandle xfront_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch < OpenVolumeMesh<VecT>::cells_.size());

        return cell(_ch).halffaces()[XF];
    }

    inline HalfFaceHandle xback_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch < OpenVolumeMesh<VecT>::cells_.size());

        return cell(_ch).halffaces()[XB];
    }

    inline HalfFaceHandle yfront_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch < OpenVolumeMesh<VecT>::cells_.size());

        return cell(_ch).halffaces()[YF];
    }

    inline HalfFaceHandle yback_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch < OpenVolumeMesh<VecT>::cells_.size());

        return cell(_ch).halffaces()[YB];
    }

    inline HalfFaceHandle zfront_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch < OpenVolumeMesh<VecT>::cells_.size());

        return cell(_ch).halffaces()[ZF];
    }

    inline HalfFaceHandle zback_halfface(const CellHandle& _ch) const {

        assert((unsigned int)_ch < OpenVolumeMesh<VecT>::cells_.size());

        return cell(_ch).halffaces()[ZB];
    }

    unsigned char orientation(const HalfFaceHandle& _hfh, const CellHandle& _ch) const {

        assert((unsigned int)_ch < OpenVolumeMesh<VecT>::cells_.size());

        std::vector<HalfFaceHandle> halffaces = cell(_ch).halffaces();
        for(unsigned int i = 0; i < halffaces.size(); ++i) {
            if(halffaces[i] == _hfh) return (unsigned char)i;
        }

        return INVALID;
    }

    static inline unsigned char orthogonal_orientation(const unsigned char _o1, const unsigned char _o2) {

        if(_o1 == XF && _o2 == YF) return ZF;
        if(_o1 == XF && _o2 == YB) return ZB;
        if(_o1 == XF && _o2 == ZF) return YB;
        if(_o1 == XF && _o2 == ZB) return YF;
        if(_o1 == XB && _o2 == YF) return ZB;
        if(_o1 == XB && _o2 == YB) return ZF;
        if(_o1 == XB && _o2 == ZF) return YF;
        if(_o1 == XB && _o2 == ZB) return YB;

        if(_o1 == YF && _o2 == XF) return ZB;
        if(_o1 == YF && _o2 == XB) return ZF;
        if(_o1 == YF && _o2 == ZF) return XF;
        if(_o1 == YF && _o2 == ZB) return XB;
        if(_o1 == YB && _o2 == XF) return ZF;
        if(_o1 == YB && _o2 == XB) return ZB;
        if(_o1 == YB && _o2 == ZF) return XB;
        if(_o1 == YB && _o2 == ZB) return XF;

        if(_o1 == ZF && _o2 == YF) return XB;
        if(_o1 == ZF && _o2 == YB) return XF;
        if(_o1 == ZF && _o2 == XF) return YF;
        if(_o1 == ZF && _o2 == XB) return YB;
        if(_o1 == ZB && _o2 == YF) return XF;
        if(_o1 == ZB && _o2 == YB) return XB;
        if(_o1 == ZB && _o2 == XF) return YB;
        if(_o1 == ZB && _o2 == XB) return YF;

        return INVALID;

    }

    inline HalfFaceHandle get_oriented_halfface(const unsigned char _o, const CellHandle& _ch) const {

        if(_o == XF) return xfront_halfface(_ch);
        if(_o == XB) return xback_halfface(_ch);
        if(_o == YF) return yfront_halfface(_ch);
        if(_o == YB) return yback_halfface(_ch);
        if(_o == ZF) return zfront_halfface(_ch);
        if(_o == ZB) return zback_halfface(_ch);
        return OpenVolumeMesh<VecT>::InvalidHalfFaceHandle;
    }

    HalfFaceHandle adjacent_halfface_on_sheet(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh) const {

        if(!OpenVolumeMesh<VecT>::has_bottom_up_adjacencies()) {
            std::cerr << "No bottom-up adjacencies computed so far, could not get adjacent halfface on sheet!" << std::endl;
            return OpenVolumeMesh<VecT>::InvalidHalfFaceHandle;
        }

        HalfFaceHandle n_hf = _hfh;
        HalfEdgeHandle n_he = _heh;

        // Try the 1st way
        while(true) {
            n_hf = OpenVolumeMesh<VecT>::adjacent_halfface_in_cell(n_hf, n_he);
            if(n_hf == OpenVolumeMesh<VecT>::InvalidHalfFaceHandle) break;
            n_hf = OpenVolumeMesh<VecT>::opposite_halfface_handle(n_hf);
            if(n_hf == OpenVolumeMesh<VecT>::InvalidHalfFaceHandle) break;
            HalfEdgeHandle o_he = OpenVolumeMesh<VecT>::opposite_halfedge_handle(n_he);
            if(o_he == OpenVolumeMesh<VecT>::InvalidHalfEdgeHandle) break;
            n_hf = OpenVolumeMesh<VecT>::adjacent_halfface_in_cell(n_hf, o_he);
            if(n_hf == OpenVolumeMesh<VecT>::InvalidHalfFaceHandle) break;
            else return n_hf;
        }

        n_hf = OpenVolumeMesh<VecT>::opposite_halfface_handle(_hfh);
        n_he = OpenVolumeMesh<VecT>::opposite_halfedge_handle(_heh);

        // Try the 2nd way
        while(true) {
            n_hf = OpenVolumeMesh<VecT>::adjacent_halfface_in_cell(n_hf, n_he);
            if(n_hf == OpenVolumeMesh<VecT>::InvalidHalfFaceHandle) break;
            n_hf = OpenVolumeMesh<VecT>::opposite_halfface_handle(n_hf);
            if(n_hf == OpenVolumeMesh<VecT>::InvalidHalfFaceHandle) break;
            HalfEdgeHandle o_he = OpenVolumeMesh<VecT>::opposite_halfedge_handle(n_he);
            if(o_he == OpenVolumeMesh<VecT>::InvalidHalfEdgeHandle) break;
            n_hf = OpenVolumeMesh<VecT>::adjacent_halfface_in_cell(n_hf, o_he);
            if(n_hf == OpenVolumeMesh<VecT>::InvalidHalfFaceHandle) break;
            else return OpenVolumeMesh<VecT>::opposite_halfface_handle(n_hf);
        }

        return OpenVolumeMesh<VecT>::InvalidHalfFaceHandle;
    }

    HalfFaceHandle adjacent_halfface_on_surface(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh) const {

        for(typename OpenVolumeMesh<VecT>::HalfEdgeHalfFaceIter hehf_it = OpenVolumeMesh<VecT>::hehf_iter(_heh);
                hehf_it.valid(); ++hehf_it) {
            if(*hehf_it == _hfh) continue;
            if(OpenVolumeMesh<VecT>::is_boundary(*hehf_it)) {
                return *hehf_it;
            }
            if(OpenVolumeMesh<VecT>::is_boundary(OpenVolumeMesh<VecT>::opposite_halfface_handle(*hehf_it))) {
                return OpenVolumeMesh<VecT>::opposite_halfface_handle(*hehf_it);
            }
        }
        return OpenVolumeMesh<VecT>::InvalidHalfFaceHandle;
    }

    HalfFaceHandle neighboring_outside_halfface(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh) const {

        if(!OpenVolumeMesh<VecT>::has_bottom_up_adjacencies()) {
            std::cerr << "No bottom-up adjacencies computed so far, could not get neighboring outside halfface!" << std::endl;
            return OpenVolumeMesh<VecT>::InvalidHalfFaceHandle;
        }

        for(typename OpenVolumeMesh<VecT>::HalfEdgeHalfFaceIter hehf_it = OpenVolumeMesh<VecT>::hehf_iter(_heh);
                hehf_it; ++hehf_it) {
            if(*hehf_it == _hfh) continue;
            if(OpenVolumeMesh<VecT>::is_boundary(*hehf_it)) return *hehf_it;
            if(OpenVolumeMesh<VecT>::is_boundary(OpenVolumeMesh<VecT>::opposite_halfface_handle(*hehf_it)))
                return OpenVolumeMesh<VecT>::opposite_halfface_handle(*hehf_it);
        }

        return OpenVolumeMesh<VecT>::InvalidHalfFaceHandle;
    }

private:

    const HalfFaceHandle& get_adjacent_halfface(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh,
            const std::vector<HalfFaceHandle>& _halffaces) const;

};

#if defined(INCLUDE_TEMPLATES) && !defined(OPENHEXMESHT_CC)
#include "OpenHexMeshT.cc"
#endif

#endif /* OPENHEXMESH_HH_ */
