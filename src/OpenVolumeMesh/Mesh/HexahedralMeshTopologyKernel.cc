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
 *   $Revision$                                                         *
 *   $Date$                    *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#include "HexahedralMeshTopologyKernel.hh"

namespace OpenVolumeMesh {


HexahedralMeshTopologyKernel::HexahedralMeshTopologyKernel() {

}

//========================================================================================


HexahedralMeshTopologyKernel::~HexahedralMeshTopologyKernel() {

}

//========================================================================================


FaceHandle HexahedralMeshTopologyKernel::add_face(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck) {

    if(_halfedges.size() != 4) {
        std::cerr << "Face valence is not four! Aborting." << std::endl;
        return TopologyKernel::InvalidFaceHandle;
    }

    return TopologyKernel::add_face(_halfedges, _topologyCheck);
}

//========================================================================================


FaceHandle
HexahedralMeshTopologyKernel::add_face(const std::vector<VertexHandle>& _vertices) {

    if(_vertices.size() != 4) {
        std::cerr << "Face valence is not four! Aborting." << std::endl;
        return TopologyKernel::InvalidFaceHandle;
    }

    return TopologyKernel::add_face(_vertices);
}

//========================================================================================


CellHandle
HexahedralMeshTopologyKernel::add_cell(const std::vector<HalfFaceHandle>& _halffaces, bool _topologyCheck, bool _reorderFaces) {

    if(_halffaces.size() != 6) {
        std::cerr << "Cell valence is not six! Aborting." << std::endl;
        return TopologyKernel::InvalidCellHandle;
    }
    for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
            it != _halffaces.end(); ++it) {
        if(TopologyKernel::halfface(*it).halfedges().size() != 4) {
            std::cerr << "Incident face does not have valence four! Aborting." << std::endl;
            return TopologyKernel::InvalidCellHandle;
        }
    }

    // Ordering array (see below for details)
    const int orderTop[] = {2, 4, 3, 5};
    const int orderBot[] = {3, 4, 2, 5};

    // Create new halffaces vector
    std::vector<HalfFaceHandle> ordered_halffaces;

    // The user wants the faces to be reordered
    if(_reorderFaces) {

        ordered_halffaces.resize(6, TopologyKernel::InvalidHalfFaceHandle);

        // Create top side
        ordered_halffaces[0] = _halffaces[0];

        // Go over all incident halfedges
        std::vector<HalfEdgeHandle> hes = TopologyKernel::halfface(ordered_halffaces[0]).halfedges();
        unsigned int idx = 0;
        for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin();
                he_it != hes.end(); ++he_it) {

            HalfFaceHandle ahfh = get_adjacent_halfface(ordered_halffaces[0], *he_it, _halffaces);
            if(ahfh == TopologyKernel::InvalidHalfFaceHandle) {
                std::cerr << "The current halfface is invalid!" << std::endl;
                continue;
            }
            ordered_halffaces[orderTop[idx]] = ahfh;
            ++idx;
        }

        // Now set bottom-halfface
        HalfFaceHandle cur_hf = ordered_halffaces[0];
        HalfEdgeHandle cur_he = *(TopologyKernel::halfface(cur_hf).halfedges().begin());
        cur_hf = get_adjacent_halfface(cur_hf, cur_he, _halffaces);
        cur_he = TopologyKernel::opposite_halfedge_handle(cur_he);
        cur_he = TopologyKernel::next_halfedge_in_halfface(cur_he, cur_hf);
        cur_he = TopologyKernel::next_halfedge_in_halfface(cur_he, cur_hf);
        cur_hf = get_adjacent_halfface(cur_hf, cur_he, _halffaces);

        if(cur_hf != TopologyKernel::InvalidHalfFaceHandle) {
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

    std::vector<HalfEdgeHandle> halfedgesTop = TopologyKernel::halfface(ordered_halffaces[0]).halfedges();
    std::vector<HalfEdgeHandle> halfedgesBot = TopologyKernel::halfface(ordered_halffaces[1]).halfedges();

    int offsetTop = -1;
    int offsetBot = -1;

    // Traverse halfedges top
    for(std::vector<HalfEdgeHandle>::const_iterator it = halfedgesTop.begin();
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
                return TopologyKernel::InvalidCellHandle;
            }
        }
    }

    if(offsetTop == -1) {
        std::cerr << "Faces not in right order!" << std::endl;
        return TopologyKernel::InvalidCellHandle;
    }

    // Traverse halfedges bottom
    for(std::vector<HalfEdgeHandle>::const_iterator it = halfedgesBot.begin();
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
                return TopologyKernel::InvalidCellHandle;
            }
        }
    }

    if(offsetBot == -1) {
        std::cerr << "Faces not in right order!" << std::endl;
        return TopologyKernel::InvalidCellHandle;
    }

    return TopologyKernel::add_cell(ordered_halffaces, _topologyCheck);
}

//========================================================================================

CellHandle
HexahedralMeshTopologyKernel::add_cell(const std::vector<VertexHandle>& _vertices) {

    assert(_vertices.size() == 8);

    if(!TopologyKernel::has_full_bottom_up_adjacencies()) {
        std::cerr << "Error: This function needs bottom-up adjacencies to be enabled!" << std::endl;
        return CellHandle(-1);
    }

    if(_vertices.size() != 8) {
        std::cerr << "The number of vertices is not eight!" << std::endl;
        return CellHandle(-1);
    }

    HalfFaceHandle hf0, hf1, hf2, hf3, hf4, hf5, hf6, hf7;

    std::vector<VertexHandle> vs;

    // Half-face XF
    vs.push_back(_vertices[2]);
    vs.push_back(_vertices[1]);
    vs.push_back(_vertices[0]);
    hf0 = TopologyKernel::halfface(vs); vs.clear();

    // Half-face XB
    vs.push_back(_vertices[6]);
    vs.push_back(_vertices[5]);
    vs.push_back(_vertices[4]);
    hf1 = TopologyKernel::halfface(vs); vs.clear();

    // Half-face YF
    vs.push_back(_vertices[1]);
    vs.push_back(_vertices[2]);
    vs.push_back(_vertices[6]);
    hf2 = TopologyKernel::halfface(vs); vs.clear();

    // Half-face YB
    vs.push_back(_vertices[5]);
    vs.push_back(_vertices[3]);
    vs.push_back(_vertices[0]);
    hf3 = TopologyKernel::halfface(vs); vs.clear();

    // Half-face ZF
    vs.push_back(_vertices[7]);
    vs.push_back(_vertices[4]);
    vs.push_back(_vertices[0]);
    hf4 = TopologyKernel::halfface(vs); vs.clear();

    // Half-face ZB
    vs.push_back(_vertices[3]);
    vs.push_back(_vertices[5]);
    vs.push_back(_vertices[6]);
    hf5 = TopologyKernel::halfface(vs);

    if(!hf0.is_valid()) {

        vs.clear();
        vs.push_back(_vertices[3]); vs.push_back(_vertices[2]);
        vs.push_back(_vertices[1]); vs.push_back(_vertices[0]);
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf0 = halfface_handle(fh, 0);
    }

    if(!hf1.is_valid()) {

        vs.clear();
        vs.push_back(_vertices[7]); vs.push_back(_vertices[6]);
        vs.push_back(_vertices[5]); vs.push_back(_vertices[4]);
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf1 = halfface_handle(fh, 0);
    }

    if(!hf2.is_valid()) {

        vs.clear();
        vs.push_back(_vertices[1]); vs.push_back(_vertices[2]);
        vs.push_back(_vertices[6]); vs.push_back(_vertices[7]);
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf2 = halfface_handle(fh, 0);
    }

    if(!hf3.is_valid()) {

        vs.clear();
        vs.push_back(_vertices[4]); vs.push_back(_vertices[5]);
        vs.push_back(_vertices[3]); vs.push_back(_vertices[0]);
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf3 = halfface_handle(fh, 0);
    }

    if(!hf4.is_valid()) {

        vs.clear();
        vs.push_back(_vertices[1]); vs.push_back(_vertices[7]);
        vs.push_back(_vertices[4]); vs.push_back(_vertices[0]);
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf4 = halfface_handle(fh, 0);
    }

    if(!hf5.is_valid()) {

        vs.clear();
        vs.push_back(_vertices[2]); vs.push_back(_vertices[3]);
        vs.push_back(_vertices[5]); vs.push_back(_vertices[6]);
        FaceHandle fh = TopologyKernel::add_face(vs);
        hf5 = halfface_handle(fh, 0);
    }

    assert(hf0.valid()); assert(hf1.valid()); assert(hf2.valid());
    assert(hf3.valid()); assert(hf4.valid()); assert(hf5.valid());

    std::vector<HalfFaceHandle> hfs;
    hfs.push_back(hf0); hfs.push_back(hf1); hfs.push_back(hf2);
    hfs.push_back(hf3); hfs.push_back(hf4); hfs.push_back(hf5);

    return TopologyKernel::add_cell(hfs,false);
}

//========================================================================================

const HalfFaceHandle&
HexahedralMeshTopologyKernel::get_adjacent_halfface(const HalfFaceHandle& _hfh, const HalfEdgeHandle& _heh,
        const std::vector<HalfFaceHandle>& _halffaces) const {

    // Search for halfface that is incident to the opposite
    // halfedge of _heh
    HalfEdgeHandle o_he = TopologyKernel::opposite_halfedge_handle(_heh);

    for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin();
            it != _halffaces.end(); ++it) {
        if(*it == _hfh) continue;
        std::vector<HalfEdgeHandle> halfedges = TopologyKernel::halfface(*it).halfedges();
        for(std::vector<HalfEdgeHandle>::const_iterator h_it = halfedges.begin();
                h_it != halfedges.end(); ++h_it) {
            if(*h_it == o_he) return *it;
        }
    }

    return TopologyKernel::InvalidHalfFaceHandle;
}

} // Namespace OpenVolumeMesh
