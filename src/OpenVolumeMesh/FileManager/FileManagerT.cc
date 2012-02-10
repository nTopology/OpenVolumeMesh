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

#define FILEMANAGERT_CC

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <typeinfo>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include "FileManager.hh"

namespace OpenVolumeMesh {

namespace IO {

//==================================================

template <class MeshT>
bool FileManager::readFile(const std::string& _filename, MeshT& _mesh,
    bool _topologyCheck, bool _computeBottomUpAdjacencies) const {

    std::ifstream iff(_filename.c_str(), std::ios::in);

    if(!iff.good()) {
        std::cerr << "Error: Could not open file " << _filename << " for reading!" << std::endl;
        iff.close();
        return false;
    }

    std::stringstream sstr;
    std::string line;
    std::string s_tmp;
    unsigned int c = 0u;
    typedef typename MeshT::PointT Point;
    Point v = Point(0.0, 0.0, 0.0);
    unsigned int v1 = 0; unsigned int v2 = 0;

    /*
     * Header
     */

    bool header_found = true;

    // Get first line
    getCleanLine(iff, line);
    sstr.str(line);

    // Check header
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "OVM") {
        //iff.close();
        header_found = false;
        std::cerr << "The specified file might not be in OpenVolumeMesh format!" << std::endl;
        //return false;
    }

    // Get ASCII/BINARY string
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp == "BINARY") {
        iff.close();
        std::cerr << "Binary files are not supported at the moment!" << std::endl;
        return false;
    }

    /*
     * Vertices
     */
    if(!header_found) {
        sstr.clear();
        sstr.str(line);
    } else {
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
    }

    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "VERTICES") {
        iff.close();
        std::cerr << "No vertex section defined!" << std::endl;
        return false;
    } else {

        // Read in number of vertices
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        // Read in vertices
        for(unsigned int i = 0u; i < c; ++i) {

            getCleanLine(iff, line);
            sstr.clear();
            sstr.str(line);
            sstr >> v[0];
            sstr >> v[1];
            sstr >> v[2];
            _mesh.add_vertex(v);
        }
    }

    /*
     * Edges
     */
    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "EDGES") {
        iff.close();
        std::cerr << "No edge section defined!" << std::endl;
        return false;
    } else {

        // Read in number of edges
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        // Read in edges
        for(unsigned int i = 0u; i < c; ++i) {

            getCleanLine(iff, line);
            sstr.clear();
            sstr.str(line);
            sstr >> v1;
            sstr >> v2;
            _mesh.add_edge(VertexHandle(v1), VertexHandle(v2));
        }
    }

    /*
     * Faces
     */
    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "FACES") {
        iff.close();
        std::cerr << "No face section defined!" << std::endl;
        return false;
    } else {

        // Read in number of faces
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        // Read in faces
        for(unsigned int i = 0u; i < c; ++i) {

            getCleanLine(iff, line);
            sstr.clear();
            sstr.str(line);

            std::vector<HalfEdgeHandle> hes;

            // Get face valence
            unsigned int val = 0u;
            sstr >> val;

            // Read half-edge indices
            for(unsigned int e = 0; e < val; ++e) {
                sstr >> v1;
                hes.push_back(HalfEdgeHandle(v1));
            }

            _mesh.add_face(hes, _topologyCheck);
        }
    }

    /*
     * Cells
     */
    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "POLYHEDRA") {
        iff.close();
        std::cerr << "No polyhedra section defined!" << std::endl;
        return false;
    } else {

        // Read in number of cells
        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        // Read in cells
        for(unsigned int i = 0u; i < c; ++i) {

            getCleanLine(iff, line);
            sstr.clear();
            sstr.str(line);

            std::vector<HalfFaceHandle> hfs;

            // Get cell valence
            unsigned int val = 0u;
            sstr >> val;

            // Read half-face indices
            for(unsigned int f = 0; f < val; ++f) {
                sstr >> v1;
                hfs.push_back(HalfFaceHandle(v1));
            }

            _mesh.add_cell(hfs, _topologyCheck);
        }
    }

    iff.close();

    // Compute top-down-adjacencies
    if(_computeBottomUpAdjacencies) {
        _mesh.update_adjacencies();
    }

    std::cerr << "######## openvolumemesh info #########" << std::endl;
    std::cerr << "#vertices: " << _mesh.n_vertices() << std::endl;
    std::cerr << "#edges:    " << _mesh.n_edges() << std::endl;
    std::cerr << "#faces:    " << _mesh.n_faces() << std::endl;
    std::cerr << "#cells:    " << _mesh.n_cells() << std::endl;
    std::cerr << "######################################" << std::endl;

    return true;
}

//==================================================

template<class MeshT>
bool FileManager::writeFile(const std::string& _filename, const MeshT& _mesh) const {

    typedef typename MeshT::Face Face;
    std::ofstream off(_filename.c_str(), std::ios::out);

    if(!off.good()) {
        std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
        off.close();
        return false;
    }

    // Write header
    off << "OVM ASCII" << std::endl;

    unsigned int n_vertices(_mesh.n_vertices());
    off << "Vertices" << std::endl;
    off << n_vertices << std::endl;

    typedef typename MeshT::PointT Point;

    // write vertices
    for(VertexIter v_it = _mesh.v_iter(); v_it; ++v_it) {

        Point v = _mesh.vertex(*v_it);
        off << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    unsigned int n_edges(_mesh.n_edges());
    off << "Edges" << std::endl;
    off << n_edges << std::endl;

    // write edges
    for(EdgeIter e_it = _mesh.e_iter(); e_it; ++e_it) {

        VertexHandle from_vertex = _mesh.edge(*e_it).from_vertex();
        VertexHandle to_vertex = _mesh.edge(*e_it).to_vertex();
        off << from_vertex << " " << to_vertex << std::endl;
    }

    unsigned int n_faces(_mesh.n_faces());
    off << "Faces" << std::endl;
    off << n_faces << std::endl;

    // write faces
    for(FaceIter f_it = _mesh.f_iter(); f_it; ++f_it) {

        off << _mesh.face(*f_it).halfedges().size() << " ";

        std::vector<HalfEdgeHandle> halfedges = _mesh.face(*f_it).halfedges();

        for(typename std::vector<HalfEdgeHandle>::const_iterator it = halfedges.begin(); it
                != halfedges.end(); ++it) {

            off << it->idx();

            if((it + 1) != halfedges.end())
                off << " ";
        }

        off << std::endl;
    }

    unsigned int n_cells(_mesh.n_cells());
    off << "Polyhedra" << std::endl;
    off << n_cells << std::endl;

    for(CellIter c_it = _mesh.c_iter(); c_it; ++c_it) {

        off << _mesh.cell(*c_it).halffaces().size() << " ";

        std::vector<HalfFaceHandle> halffaces = _mesh.cell(*c_it).halffaces();

        for(typename std::vector<HalfFaceHandle>::const_iterator it = halffaces.begin(); it
                != halffaces.end(); ++it) {

            off << it->idx();

            if((it + 1) != halffaces.end())
                off << " ";
        }

        off << std::endl;
    }

    off.close();

    return true;
}

//==================================================

} // Namespace IO

} // Namespace FileManager
