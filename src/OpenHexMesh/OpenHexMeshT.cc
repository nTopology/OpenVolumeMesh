/*
 * OpenHexMeshT.cc
 *
 *  Created on: Jun 10, 2011
 *      Author: kremer
 */

#define OPENVOLUMEMESHT_CC

#include "OpenHexMesh.hh"

template <typename VecT>
OpenHexMesh<VecT>::OpenHexMesh() {

}

//========================================================================================

template <typename VecT>
OpenHexMesh<VecT>::~OpenHexMesh() {

}

//========================================================================================

template <typename VecT>
typename OpenHexMesh<VecT>::FaceHandle
OpenHexMesh<VecT>::add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck) {

    if(_halfedges.size() != 4) {
        std::cerr << "Face valence is not four! Aborting." << std::endl;
        return OpenVolumeMesh<VecT>::InvalidFaceHandle;
    }

    return OpenVolumeMesh<VecT>::add_face(_halfedges, _topologyCheck);
}

//========================================================================================

template <typename VecT>
typename OpenHexMesh<VecT>::FaceHandle
OpenHexMesh<VecT>::add_face(const std::vector<VertexHandle>& _vertices) {

    if(_vertices.size() != 4) {
        std::cerr << "Face valence is not four! Aborting." << std::endl;
        return OpenVolumeMesh<VecT>::InvalidFaceHandle;
    }

    return OpenVolumeMesh<VecT>::add_face(_vertices);
}

//========================================================================================

template <typename VecT>
typename OpenHexMesh<VecT>::CellHandle
OpenHexMesh<VecT>::add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck, bool _reorderFaces) {

    if(_halffaces.size() != 6) {
        std::cerr << "Cell valence is not six! Aborting." << std::endl;
        return OpenVolumeMesh<VecT>::InvalidCellHandle;
    }
    for(typename std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
            it != _halffaces.end(); ++it) {
        if(halfface(*it).halfedges().size() != 4) {
            std::cerr << "Incident face does not have valence four! Aborting." << std::endl;
            return OpenVolumeMesh<VecT>::InvalidCellHandle;
        }
    }

    // Ordering array (see below for details)
    const int orderTop[] = {2, 4, 3, 5};
    const int orderBot[] = {3, 4, 2, 5};

    // Create new halffaces vector
    std::vector<HalfFaceHandle> ordered_halffaces;

    // The user wants the faces to be reordered
    if(_reorderFaces) {

        ordered_halffaces.resize(6, OpenVolumeMesh<VecT>::InvalidHalfFaceHandle);

        // Create top side
        ordered_halffaces[0] = _halffaces[0];

        // Go over all incident halfedges
        std::vector<HalfEdgeHandle> hes = halfface(ordered_halffaces[0]).halfedges();
        unsigned int idx = 0;
        for(typename std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin();
                he_it != hes.end(); ++he_it) {

            HalfFaceHandle ahfh = get_adjacent_halfface(ordered_halffaces[0], *he_it, _halffaces);
            if(ahfh == OpenVolumeMesh<VecT>::InvalidHalfFaceHandle) {
                std::cerr << "The current halfface is invalid!" << std::endl;
                continue;
            }
            ordered_halffaces[orderTop[idx]] = ahfh;
            ++idx;
        }

        // Now set bottom-halfface
        HalfFaceHandle cur_hf = ordered_halffaces[0];
        HalfEdgeHandle cur_he = *(halfface(cur_hf).halfedges().begin());
        cur_hf = get_adjacent_halfface(cur_hf, cur_he, _halffaces);
        cur_he = opposite_halfedge_handle(cur_he);
        cur_he = next_halfedge_in_halfface(cur_he, cur_hf);
        cur_he = next_halfedge_in_halfface(cur_he, cur_hf);
        cur_hf = get_adjacent_halfface(cur_hf, cur_he, _halffaces);

        if(cur_hf != OpenVolumeMesh<VecT>::InvalidHalfFaceHandle) {
            ordered_halffaces[1] = cur_hf;
        } else {
            std::cerr << "The current halfface is invalid!" << std::endl;
        }

    } else {
        // Just copy the original ones
        ordered_halffaces = _halffaces;
    }

    // Now check if faces are in right order
    /*
     * The test works as follows: Test for both the first and second face in the list,
     * whether the following order holds (clockwise):
     *
     * View from above, outside.
     *           ____
     *          | 4  |
     *      ____|____|____
     *     | 5  | 1  | 6  |
     *     |____|____|____|
     *          | 3  |
     *          |____|
     *
     * View from below, outside.
     *           ____
     *          | 3  |
     *      ____|____|____
     *     | 5  | 2  | 6  |
     *     |____|____|____|
     *          | 4  |
     *          |____|
     */

    HalfFaceHandle hfhTop = ordered_halffaces[0];
    HalfFaceHandle hfhBot = ordered_halffaces[1];

    std::vector<HalfEdgeHandle> halfedgesTop = halfface(ordered_halffaces[0]).halfedges();
    std::vector<HalfEdgeHandle> halfedgesBot = halfface(ordered_halffaces[1]).halfedges();

    int offsetTop = -1;
    int offsetBot = -1;

    // Traverse halfedges top
    for(typename std::vector<HalfEdgeHandle>::const_iterator it = halfedgesTop.begin();
            it != halfedgesTop.end(); ++it) {

        HalfFaceHandle ahfh = get_adjacent_halfface(hfhTop, *it, ordered_halffaces);

        if(offsetTop == -1) {
            if(ahfh == ordered_halffaces[2])       offsetTop = 0;
            else if(ahfh == ordered_halffaces[4])  offsetTop = 1;
            else if(ahfh == ordered_halffaces[3])  offsetTop = 2;
            else if(ahfh == ordered_halffaces[5])  offsetTop = 3;
        } else {
            offsetTop = (offsetTop + 1) % 4;
            if(ahfh != ordered_halffaces[orderTop[offsetTop]]) {
                std::cerr << "Faces not in right order!" << std::endl;
                return OpenVolumeMesh<VecT>::InvalidCellHandle;
            }
        }
    }

    if(offsetTop == -1) {
        std::cerr << "Faces not in right order!" << std::endl;
        return OpenVolumeMesh<VecT>::InvalidCellHandle;
    }

    // Traverse halfedges bottom
    for(typename std::vector<HalfEdgeHandle>::const_iterator it = halfedgesBot.begin();
            it != halfedgesBot.end(); ++it) {

        HalfFaceHandle ahfh = get_adjacent_halfface(hfhBot, *it, ordered_halffaces);

        if(offsetBot == -1) {
            if(ahfh == ordered_halffaces[3])       offsetBot = 0;
            else if(ahfh == ordered_halffaces[4])  offsetBot = 1;
            else if(ahfh == ordered_halffaces[2])  offsetBot = 2;
            else if(ahfh == ordered_halffaces[5])  offsetBot = 3;
        } else {
            offsetBot = (offsetBot + 1) % 4;
            if(ahfh != ordered_halffaces[orderBot[offsetBot]]) {
                std::cerr << "Faces not in right order!" << std::endl;
                return OpenVolumeMesh<VecT>::InvalidCellHandle;
            }
        }
    }

    if(offsetBot == -1) {
        std::cerr << "Faces not in right order!" << std::endl;
        return OpenVolumeMesh<VecT>::InvalidCellHandle;
    }

    return OpenVolumeMesh<VecT>::add_cell(ordered_halffaces, _topologyCheck);
}

//========================================================================================

template <typename VecT>
const typename OpenHexMesh<VecT>::HalfFaceHandle&
OpenHexMesh<VecT>::get_adjacent_halfface(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh,
        const std::vector<HalfFaceHandle>& _halffaces) const {

    // Search for halfface that is incident to the opposite
    // halfedge of _heh
    HalfEdgeHandle o_he = opposite_halfedge_handle(_heh);

    for(typename std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
            it != _halffaces.end(); ++it) {
        if(*it == _hfh) continue;
        std::vector<HalfEdgeHandle> halfedges = halfface(*it).halfedges();
        for(typename std::vector<HalfEdgeHandle>::const_iterator h_it = halfedges.begin();
                h_it != halfedges.end(); ++h_it) {
            if(*h_it == o_he) return *it;
        }
    }

    return OpenVolumeMesh<VecT>::InvalidHalfFaceHandle;
}
