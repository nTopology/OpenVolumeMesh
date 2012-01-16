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

#define FILEMANAGERT_CC

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <typeinfo>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/PolyhedralMesh/PolyhedralMesh.hh>

#include "FileManager.hh"

namespace OpenVolumeMesh {

namespace IO {

//==================================================

FileManager::FileManager() {

}

//==================================================

FileManager::~FileManager() {

}

//==================================================

void FileManager::trimString(std::string& _string) const {

    // Trim Both leading and trailing spaces
    size_t start = _string.find_first_not_of(" \t\r\n");
    size_t end = _string.find_last_not_of(" \t\r\n");

    if((std::string::npos == start) || (std::string::npos == end)) {
        _string = "";
    } else {
        _string = _string.substr(start, end - start + 1);
    }
}

//==================================================

void FileManager::extractQuotedText(std::string& _string) const {

    // Trim Both leading and trailing quote marks
    size_t start = _string.find_first_of("\""); ++start;
    size_t end = _string.find_last_not_of("\"");

    if((std::string::npos == start) || (std::string::npos == end)) {
        _string = "";
    } else {
        _string = _string.substr(start, end - start + 1);
    }
}

//==================================================

bool FileManager::getCleanLine(std::istream& _ifs, std::string& _string, bool _skipEmptyLines) const {

    // While we are not at the end of the file
    while(true) {

        // Get the current line:
        std::getline(_ifs, _string);

        // Remove whitespace at beginning and end
        trimString(_string);

        // Check if string is not empty ( otherwise we continue
        if(_string.size() != 0) {

            // Check if string is a comment ( starting with # )
            if(_string[0] != '#') {
                return true;
            }

        } else {
            if(!_skipEmptyLines)
                return true;
        }

        if(_ifs.eof()) {
            std::cerr << "End of file reached while searching for input!" << std::endl;
            return false;
        }
    }

    return false;
}

//==================================================

template <class MeshT>
bool FileManager::readFile(const std::string& _filename, MeshT& _mesh,
    bool _topologyCheck, bool _computeBottomUpAdjacencies, bool _computeFaceNormals) const {

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
    typedef OpenVolumeMesh::Geometry::Vec3d Vec3d;
    Vec3d v = Vec3d(0.0, 0.0, 0.0);
    unsigned int v1 = 0; unsigned int v2 = 0;

    /*
     * Header
     */

    // Get first line
    getCleanLine(iff, line);
    sstr.str(line);

    // Check header
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "OVM") {
        iff.close();
        std::cerr << "The specified file is not in OpenVolumeMesh format!" << std::endl;
        return false;
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
    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);

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

    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);

    parseVertexProperties<MeshT, typename MeshT::VertexIter>
        (iff, sstr, "VERTEX_PROPERTY", _mesh, _mesh.vertices_begin(), _mesh.vertices_end());

    /*
     * Edges
     */
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
            _mesh.add_edge(typename MeshT::VertexHandle(v1), typename MeshT::VertexHandle(v2));
        }
    }

    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);

    parseEdgeProperties<MeshT, typename MeshT::EdgeIter>
        (iff, sstr, "EDGE_PROPERTY", _mesh, _mesh.edges_begin(), _mesh.edges_end());
    parseHalfEdgeProperties<MeshT, typename MeshT::HalfEdgeIter>
        (iff, sstr, "HALFEDGE_PROPERTY", _mesh, _mesh.halfedges_begin(), _mesh.halfedges_end());

    /*
     * Faces
     */
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

            std::vector<typename MeshT::HalfEdgeHandle> hes;

            // Get face valence
            unsigned int val = 0u;
            sstr >> val;

            // Read half-edge indices
            for(unsigned int e = 0; e < val; ++e) {
                sstr >> v1;
                hes.push_back(typename MeshT::HalfEdgeHandle(v1));
            }

            _mesh.add_face(hes, _topologyCheck);
        }
    }

    getCleanLine(iff, line);
    sstr.clear();
    sstr.str(line);

    parseFaceProperties<MeshT, typename MeshT::FaceIter>
        (iff, sstr, "FACE_PROPERTY", _mesh, _mesh.faces_begin(), _mesh.faces_end());
    parseHalfFaceProperties<MeshT, typename MeshT::HalfFaceIter>
        (iff, sstr, "HALFFACE_PROPERTY", _mesh, _mesh.halffaces_begin(), _mesh.halffaces_end());

    /*
     * Cells
     */
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

            std::vector<typename MeshT::HalfFaceHandle> hfs;

            // Get cell valence
            unsigned int val = 0u;
            sstr >> val;

            // Read half-face indices
            for(unsigned int f = 0; f < val; ++f) {
                sstr >> v1;
                hfs.push_back(typename MeshT::HalfFaceHandle(v1));
            }

            _mesh.add_cell(hfs, _topologyCheck);
        }
    }

    if(!iff.eof()) {

        getCleanLine(iff, line);
        sstr.clear();
        sstr.str(line);

        parseCellProperties<MeshT, typename MeshT::CellIter>
            (iff, sstr, "POLYHEDRON_PROPERTY", _mesh, _mesh.cells_begin(), _mesh.cells_end());
    }

    iff.close();

    // Compute top-down-adjacencies
    if(_computeBottomUpAdjacencies) {
        _mesh.update_adjacencies();
    }

    // Compute face normals
    if(_computeFaceNormals) {
        _mesh.request_face_normals();
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

template<class MeshT, typename IterT>
void FileManager::parseVertexProperties(std::ifstream& _iff, std::stringstream& _sstr, const std::string& _identifier, MeshT& _mesh,
                                  const IterT& _begin, const IterT& _end) const {

    if(_iff.eof()) return;

    std::string line = _sstr.str();
    std::string s_tmp;

    // Check for properties
    _sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);

    // No property found
    if(s_tmp != _identifier) {
        // Reset stream
        _sstr.clear();
        _sstr.str(line);
        return;
    }

    while(s_tmp == _identifier && !_iff.eof()) {

        s_tmp = _sstr.str();
        extractQuotedText(s_tmp);

        std::string type;
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        _sstr >> type;

        // Add and initialize property
        if(type == "int") {
            initializeProperty<MeshT, VPropHandleT<int>, int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned int") {
            initializeProperty<MeshT, VPropHandleT<unsigned int>, unsigned int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "float") {
            initializeProperty<MeshT, VPropHandleT<float>, float, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "double") {
            initializeProperty<MeshT, VPropHandleT<double>, double, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "char") {
            initializeProperty<MeshT, VPropHandleT<char>, char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned char") {
            initializeProperty<MeshT, VPropHandleT<unsigned char>, unsigned char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "bool") {
            initializeProperty<MeshT, VPropHandleT<bool>, bool, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "string") {
            initializeProperty<MeshT, VPropHandleT<std::string>, std::string, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else {
            std::cerr << "Data type '" << type << "' not recognized!" << std::endl;

            // Skip number of entities lines
            for(IterT it = _begin; it != _end; ++it) {
                getCleanLine(_iff, line);
            }
        }

        // Get next line
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        // Check for properties
        _sstr >> s_tmp;
        std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    }

    _sstr.clear();
    _sstr.str(line);
}

//==================================================

template<class MeshT, typename IterT>
void FileManager::parseEdgeProperties(std::ifstream& _iff, std::stringstream& _sstr, const std::string& _identifier, MeshT& _mesh,
                                  const IterT& _begin, const IterT& _end) const {

    if(_iff.eof()) return;

    std::string line = _sstr.str();
    std::string s_tmp;

    // Check for properties
    _sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);

    // No property found
    if(s_tmp != _identifier) {
        // Reset stream
        _sstr.clear();
        _sstr.str(line);
        return;
    }

    while(s_tmp == _identifier && !_iff.eof()) {

        s_tmp = _sstr.str();
        extractQuotedText(s_tmp);

        std::string type;
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        _sstr >> type;

        // Add and initialize property
        if(type == "int") {
            initializeProperty<MeshT, EPropHandleT<int>, int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned int") {
            initializeProperty<MeshT, EPropHandleT<unsigned int>, unsigned int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "float") {
            initializeProperty<MeshT, EPropHandleT<float>, float, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "double") {
            initializeProperty<MeshT, EPropHandleT<double>, double, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "char") {
            initializeProperty<MeshT, EPropHandleT<char>, char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned char") {
            initializeProperty<MeshT, EPropHandleT<unsigned char>, unsigned char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "bool") {
            initializeProperty<MeshT, EPropHandleT<bool>, bool, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "string") {
            initializeProperty<MeshT, EPropHandleT<std::string>, std::string, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else {
            std::cerr << "Data type '" << type << "' not recognized!" << std::endl;

            // Skip number of entities lines
            for(IterT it = _begin; it != _end; ++it) {
                getCleanLine(_iff, line);
            }
        }

        // Get next line
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        // Check for properties
        _sstr >> s_tmp;
        std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    }

    _sstr.clear();
    _sstr.str(line);
}

//==================================================

template<class MeshT, typename IterT>
void FileManager::parseHalfEdgeProperties(std::ifstream& _iff, std::stringstream& _sstr, const std::string& _identifier, MeshT& _mesh,
                                  const IterT& _begin, const IterT& _end) const {

    if(_iff.eof()) return;

    std::string line = _sstr.str();
    std::string s_tmp;

    // Check for properties
    _sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);

    // No property found
    if(s_tmp != _identifier) {
        // Reset stream
        _sstr.clear();
        _sstr.str(line);

        return;
    }

    while(s_tmp == _identifier && !_iff.eof()) {

        s_tmp = _sstr.str();
        extractQuotedText(s_tmp);

        std::string type;
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        _sstr >> type;

        // Add and initialize property
        if(type == "int") {
            initializeProperty<MeshT, HEPropHandleT<int>, int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned int") {
            initializeProperty<MeshT, HEPropHandleT<unsigned int>, unsigned int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "float") {
            initializeProperty<MeshT, HEPropHandleT<float>, float, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "double") {
            initializeProperty<MeshT, HEPropHandleT<double>, double, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "char") {
            initializeProperty<MeshT, HEPropHandleT<char>, char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned char") {
            initializeProperty<MeshT, HEPropHandleT<unsigned char>, unsigned char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "bool") {
            initializeProperty<MeshT, HEPropHandleT<bool>, bool, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "string") {
            initializeProperty<MeshT, HEPropHandleT<std::string>, std::string, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else {
            std::cerr << "Data type '" << type << "' not recognized!" << std::endl;

            // Skip number of entities lines
            for(IterT it = _begin; it != _end; ++it) {
                getCleanLine(_iff, line);
            }
        }

        // Get next line
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        // Check for properties
        _sstr >> s_tmp;
        std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    }

    _sstr.clear();
    _sstr.str(line);
}

//==================================================

template<class MeshT, typename IterT>
void FileManager::parseFaceProperties(std::ifstream& _iff, std::stringstream& _sstr, const std::string& _identifier, MeshT& _mesh,
                                  const IterT& _begin, const IterT& _end) const {

    if(_iff.eof()) return;

    std::string line = _sstr.str();
    std::string s_tmp;

    // Check for properties
    _sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);

    // No property found
    if(s_tmp != _identifier) {
        // Reset stream
        _sstr.clear();
        _sstr.str(line);
        return;
    }

    while(s_tmp == _identifier && !_iff.eof()) {

        s_tmp = _sstr.str();
        extractQuotedText(s_tmp);

        std::string type;
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        _sstr >> type;

        // Add and initialize property
        if(type == "int") {
            initializeProperty<MeshT, FPropHandleT<int>, int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned int") {
            initializeProperty<MeshT, FPropHandleT<unsigned int>, unsigned int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "float") {
            initializeProperty<MeshT, FPropHandleT<float>, float, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "double") {
            initializeProperty<MeshT, FPropHandleT<double>, double, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "char") {
            initializeProperty<MeshT, FPropHandleT<char>, char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned char") {
            initializeProperty<MeshT, FPropHandleT<unsigned char>, unsigned char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "bool") {
            initializeProperty<MeshT, FPropHandleT<bool>, bool, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "string") {
            initializeProperty<MeshT, FPropHandleT<std::string>, std::string, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else {
            std::cerr << "Data type '" << type << "' not recognized!" << std::endl;

            // Skip number of entities lines
            for(IterT it = _begin; it != _end; ++it) {
                getCleanLine(_iff, line);
            }
        }

        // Get next line
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        // Check for properties
        _sstr >> s_tmp;
        std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    }

    _sstr.clear();
    _sstr.str(line);
}

//==================================================

template<class MeshT, typename IterT>
void FileManager::parseHalfFaceProperties(std::ifstream& _iff, std::stringstream& _sstr, const std::string& _identifier, MeshT& _mesh,
                                  const IterT& _begin, const IterT& _end) const {

    if(_iff.eof()) return;

    std::string line = _sstr.str();
    std::string s_tmp;

    // Check for properties
    _sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);

    // No property found
    if(s_tmp != _identifier) {
        // Reset stream
        _sstr.clear();
        _sstr.str(line);
        return;
    }

    while(s_tmp == _identifier && !_iff.eof()) {

        s_tmp = _sstr.str();
        extractQuotedText(s_tmp);

        std::string type;
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        _sstr >> type;

        // Add and initialize property
        if(type == "int") {
            initializeProperty<MeshT, HFPropHandleT<int>, int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned int") {
            initializeProperty<MeshT, HFPropHandleT<unsigned int>, unsigned int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "float") {
            initializeProperty<MeshT, HFPropHandleT<float>, float, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "double") {
            initializeProperty<MeshT, HFPropHandleT<double>, double, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "char") {
            initializeProperty<MeshT, HFPropHandleT<char>, char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned char") {
            initializeProperty<MeshT, HFPropHandleT<unsigned char>, unsigned char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "bool") {
            initializeProperty<MeshT, HFPropHandleT<bool>, bool, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "string") {
            initializeProperty<MeshT, HFPropHandleT<std::string>, std::string, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else {
            std::cerr << "Data type '" << type << "' not recognized!" << std::endl;

            // Skip number of entities lines
            for(IterT it = _begin; it != _end; ++it) {
                getCleanLine(_iff, line);
            }
        }

        // Get next line
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        // Check for properties
        _sstr >> s_tmp;
        std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    }

    _sstr.clear();
    _sstr.str(line);
}

//==================================================

template<class MeshT, typename IterT>
void FileManager::parseCellProperties(std::ifstream& _iff, std::stringstream& _sstr, const std::string& _identifier, MeshT& _mesh,
                                  const IterT& _begin, const IterT& _end) const {

    if(_iff.eof()) return;

    std::string line = _sstr.str();
    std::string s_tmp;

    // Check for properties
    _sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);

    // No property found
    if(s_tmp != _identifier) {
        // Reset stream
        _sstr.clear();
        _sstr.str(line);
        return;
    }

    while(s_tmp == _identifier && !_iff.eof()) {

        s_tmp = _sstr.str();
        extractQuotedText(s_tmp);

        std::string type;
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        _sstr >> type;

        // Add and initialize property
        if(type == "int") {
            initializeProperty<MeshT, CPropHandleT<int>, int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned int") {
            initializeProperty<MeshT, CPropHandleT<unsigned int>, unsigned int, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "float") {
            initializeProperty<MeshT, CPropHandleT<float>, float, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "double") {
            initializeProperty<MeshT, CPropHandleT<double>, double, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "char") {
            initializeProperty<MeshT, CPropHandleT<char>, char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "unsigned char") {
            initializeProperty<MeshT, CPropHandleT<unsigned char>, unsigned char, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "bool") {
            initializeProperty<MeshT, CPropHandleT<bool>, bool, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else if(type == "string") {
            initializeProperty<MeshT, CPropHandleT<std::string>, std::string, IterT>
            (_iff, _mesh, s_tmp, _begin, _end);
        } else {
            std::cerr << "Data type '" << type << "' not recognized!" << std::endl;

            // Skip number of entities lines
            for(IterT it = _begin; it != _end; ++it) {
                getCleanLine(_iff, line);
            }
        }

        // Get next line
        getCleanLine(_iff, line);
        _sstr.clear();
        _sstr.str(line);

        // Check for properties
        _sstr >> s_tmp;
        std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    }

    _sstr.clear();
    _sstr.str(line);
}

//==================================================

template<class MeshT, typename PropHandleT, typename PropT, typename IterT>
void FileManager::initializeProperty(std::ifstream& _iff, MeshT& _mesh, const std::string& _s_tmp,
                                     const IterT& _begin, const IterT& _end) const {

    std::string line;
    std::stringstream sstr;

    PropHandleT handle;
    _mesh.add_property(handle, _s_tmp);

    PropT c;

    IterT it = _begin;
    for(; it != _end; ++it) {

        getCleanLine(_iff, line);
        sstr.clear();
        sstr.str(line);
        sstr >> c;

        _mesh.property(handle, *it) = c;
    }
}

//==================================================

template <class MeshT>
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

  // write vertices
  for (typename MeshT::VertexIter v_it = _mesh.v_iter(); v_it; ++v_it) {

    Vec3d v = _mesh.vertex(*v_it).position();
    off << v[0] << " " << v[1] << " " << v[2] << std::endl;
  }

  // Write vertex props
  writePropertySection<typename MeshT::const_prop_iterator>
                      (off, "Vertex_Property", _mesh.vprops_begin(), _mesh.vprops_end());

  unsigned int n_edges(_mesh.n_edges());
  off << "Edges" << std::endl;
  off << n_edges << std::endl;

  // write edges
  for (typename MeshT::EdgeIter e_it = _mesh.e_iter(); e_it; ++e_it) {

    typename MeshT::VertexHandle from_vertex = _mesh.edge(*e_it).from_vertex();
    typename MeshT::VertexHandle to_vertex = _mesh.edge(*e_it).to_vertex();
    off << from_vertex << " " << to_vertex << std::endl;
  }

  // Write edge props
  writePropertySection<typename MeshT::const_prop_iterator>
                      (off, "Edge_Property", _mesh.eprops_begin(), _mesh.eprops_end());
  // Write half-edge props
  writePropertySection<typename MeshT::const_prop_iterator>
                      (off, "HalfEdge_Property", _mesh.heprops_begin(), _mesh.heprops_end());

  unsigned int n_faces(_mesh.n_faces());
  off << "Faces" << std::endl;
  off << n_faces << std::endl;

  // write faces
  for (typename MeshT::FaceIter f_it = _mesh.f_iter(); f_it; ++f_it) {

    off << _mesh.face(*f_it).halfedges().size() << " ";

    std::vector<typename MeshT::HalfEdgeHandle> halfedges = _mesh.face(*f_it).halfedges();

    for (typename std::vector<typename MeshT::HalfEdgeHandle>::const_iterator it = halfedges.begin(); it
        != halfedges.end(); ++it) {

      off << it->idx();

      if ((it + 1) != halfedges.end())
        off << " ";
    }

    off << std::endl;
  }

  // Write face props
  writePropertySection<typename MeshT::const_prop_iterator>
                      (off, "Face_Property", _mesh.fprops_begin(), _mesh.fprops_end());
  // Write half-face props
  writePropertySection<typename MeshT::const_prop_iterator>
                      (off, "HalfFace_Property", _mesh.hfprops_begin(), _mesh.hfprops_end());

  unsigned int n_cells(_mesh.n_cells());
  off << "Polyhedra" << std::endl;
  off << n_cells << std::endl;

  for (typename MeshT::CellIter c_it = _mesh.c_iter(); c_it; ++c_it) {

    off << _mesh.cell(*c_it).halffaces().size() << " ";

    std::vector<typename MeshT::HalfFaceHandle> halffaces = _mesh.cell(*c_it).halffaces();

    for (typename std::vector<typename MeshT::HalfFaceHandle>::const_iterator it = halffaces.begin(); it
        != halffaces.end(); ++it) {

      off << it->idx();

      if ((it + 1) != halffaces.end())
        off << " ";
    }

    off << std::endl;
  }

  // Write polyhedron props
  writePropertySection<typename MeshT::const_prop_iterator>
                      (off, "Polyhedron_Property", _mesh.cprops_begin(), _mesh.cprops_end());

  off.close();

  return true;
}

//==================================================

template<typename PropIterT>
void FileManager::writePropertySection(std::ofstream& _ofs, const std::string& _identifier,
                            const PropIterT& _begin, const PropIterT& _end) const {

    for(PropIterT it = _begin; it != _end; ++it) {

        // Write property header
        _ofs << _identifier << " \"" << (*it)->name() << "\"" << std::endl;

        typedef OpenVolumeMeshPropertyT<int> P_I;
        typedef OpenVolumeMeshPropertyT<unsigned int> P_UI;
        typedef OpenVolumeMeshPropertyT<float> P_F;
        typedef OpenVolumeMeshPropertyT<double> P_D;
        typedef OpenVolumeMeshPropertyT<char> P_C;
        typedef OpenVolumeMeshPropertyT<unsigned char> P_UC;
        typedef OpenVolumeMeshPropertyT<bool> P_B;
        typedef OpenVolumeMeshPropertyT<std::string> P_S;

        P_I*  prop_i =  dynamic_cast<P_I*> (*it);
        P_UI* prop_ui = dynamic_cast<P_UI*>(*it);
        P_F*  prop_f =  dynamic_cast<P_F*> (*it);
        P_D*  prop_d =  dynamic_cast<P_D*> (*it);
        P_C*  prop_c =  dynamic_cast<P_C*> (*it);
        P_UC* prop_uc = dynamic_cast<P_UC*>(*it);
        P_B*  prop_b =  dynamic_cast<P_B*> (*it);
        P_S*  prop_s =  dynamic_cast<P_S*> (*it);

        if(prop_i != NULL) {
            _ofs << "int" << std::endl;
            for(typename P_I::vector_type::const_iterator it = prop_i->data_vector().begin();
                    it != prop_i->data_vector().end(); ++it)
                _ofs << *it << std::endl;
        } else if(prop_ui != NULL) {
            _ofs << "unsigned int" << std::endl;
            for(typename P_UI::vector_type::const_iterator it = prop_ui->data_vector().begin();
                    it != prop_ui->data_vector().end(); ++it)
                _ofs << *it << std::endl;
        } else if(prop_f != NULL) {
            _ofs << "float" << std::endl;
            for(typename P_F::vector_type::const_iterator it = prop_f->data_vector().begin();
                    it != prop_f->data_vector().end(); ++it)
                _ofs << *it << std::endl;
        } else if(prop_d != NULL) {
            _ofs << "double" << std::endl;
            for(typename P_D::vector_type::const_iterator it = prop_d->data_vector().begin();
                    it != prop_d->data_vector().end(); ++it)
                _ofs << *it << std::endl;
        } else if(prop_c != NULL) {
            _ofs << "char" << std::endl;
            for(typename P_C::vector_type::const_iterator it = prop_c->data_vector().begin();
                    it != prop_c->data_vector().end(); ++it)
                _ofs << *it << std::endl;
        } else if(prop_uc != NULL) {
            _ofs << "unsigned char" << std::endl;
            for(typename P_UC::vector_type::const_iterator it = prop_uc->data_vector().begin();
                    it != prop_uc->data_vector().end(); ++it)
                _ofs << *it << std::endl;
        } else if(prop_b != NULL) {
            _ofs << "bool" << std::endl;
            for(unsigned int i = 0; i < prop_b->n_elements(); ++i)
                _ofs << (*prop_b)[i] << std::endl;
        } else if(prop_s != NULL) {
            _ofs << "string" << std::endl;
            for(unsigned int i = 0; i < prop_s->n_elements(); ++i)
                _ofs << (*prop_b)[i] << std::endl;
        } else {
            _ofs << "unknown" << std::endl << "0" << std::endl;
            continue;
        }
    }
}

//==================================================

bool FileManager::isHexahedralMesh(const std::string& _filename) const {

  std::ifstream iff(_filename.c_str(), std::ios::in);

  if(!iff.good()) {
    std::cerr << "Could not open file " << _filename << " for reading!" << std::endl;
    iff.close();
    return false;
  }

  std::string s;
  unsigned int n = 0u;

  // Skip until we find polyhedra section
  while (true || !iff.eof()) {
    iff >> s;
    if (s == "Polyhedra") {
      break;
    }
  }

  if (iff.eof()) {
    // Polyhedra section not found in file. Defaulting to polyhedral type.
    iff.close();
    return false;
  }

  // Read in number of cells
  iff >> n;
  unsigned int v = 0;
  char tmp[256];
  for (unsigned int i = 0; i < n; ++i) {
    iff >> v;
    iff.getline(tmp, 256);
    if (v != 6u) {
      iff.close();
      return false;
    }
  }
  iff.close();
  return true;
}

//==================================================

} // Namespace IO

} // Namespace FileManager
