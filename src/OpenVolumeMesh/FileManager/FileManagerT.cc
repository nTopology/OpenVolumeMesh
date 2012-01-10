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

#define FILEMANAGERT_CC

#include <OpenVolumeMesh/Geometry/VectorT.hh>

#include "FileManager.hh"

#include <fstream>
#include <vector>
#include <iostream>

namespace OpenVolumeMesh {

namespace IO {

//==================================================

FileManager::FileManager() {

}

//==================================================

FileManager::~FileManager() {

}

//==================================================

template <class MeshT>
bool FileManager::readFile(const std::string& _filename, MeshT& _mesh,
    bool _topologyCheck, bool _computeBottomUpAdjacencies, bool _computeFaceNormals) const {

  // found edges in file
  bool edges_in_file = false;

  typedef typename MeshT::Vertex Vertex;
  _mesh.clear();

  std::ifstream iff(_filename.c_str(), std::ios::in);

  std::string s;

  // read header
  iff >> s;
  if (s != "Vertices") {
    std::cerr << "Error reading OpenVolumeMesh file: Vertex section failed!" << std::endl;
    iff.close();
    return false;
  }

  // read vertices
  int nv = 0;
  iff >> nv;
  for (int i = 0; i < nv; ++i) {
    double x, y, z;
    iff >> x;
    iff >> y;
    iff >> z;

    Vec3d v(x, y, z);
    _mesh.add_vertex(v);
  }

  iff >> s;
  if (s != "Edges") {
    std::cerr << "No edges found!" << std::endl;
  } else {
    edges_in_file = true;
    int ne = 0;
    iff >> ne;
    for (int e = 0; e < ne; ++e) {
      int v1 = 0;
      int v2 = 0;
      iff >> v1;
      iff >> v2;
      _mesh.add_edge(typename MeshT::VertexHandle(v1), typename MeshT::VertexHandle(v2));
    }
  }

  if (edges_in_file) {
    iff >> s;
  }

  if (s != "Faces") {
    std::cerr << "Error reading OpenVolumeMesh file: Face section failed!" << std::endl;
    iff.close();
    return false;
  }

  // Read halfface indices
  int nf = 0;
  iff >> nf;
  for (int f = 0; f < nf; ++f) {
    int nfhe;
    iff >> nfhe;
    std::vector<typename MeshT::HalfEdgeHandle> hes;
    for (int he = 0; he < nfhe; ++he) {
      int i;
      iff >> i;
      hes.push_back(i);
    }

    _mesh.add_face(hes, _topologyCheck);
  }

  // Read faces and find the respective halffaces
  iff >> s;
  if (s != "Polyhedra") {
    std::cerr << "Error reading OpenVolumeMesh file: Polyhedra section failed!" << std::endl;
    iff.close();
    return false;
  }

  // Just read the specified halffaces
  int nc = 0;
  iff >> nc;
  for (int c = 0; c < nc; ++c) {

    int nhf;
    iff >> nhf;
    std::vector<typename MeshT::HalfFaceHandle> hfs;
    for (int hf = 0; hf < nhf; ++hf) {
      int i;
      iff >> i;
      hfs.push_back(i);
    }

    _mesh.add_cell(hfs, _topologyCheck);
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

template <class MeshT>
bool FileManager::writeFile(const std::string& _filename, const MeshT& _mesh) const {

  typedef typename MeshT::Face Face;
  std::ofstream off(_filename.c_str(), std::ios::out);

  if(!off.good()) {
    std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
    off.close();
    return false;
  }

  int n_vertices(_mesh.n_vertices());
  off << "Vertices" << std::endl;
  off << n_vertices << std::endl;

  // write vertices
  for (typename MeshT::VertexIter v_it = _mesh.v_iter(); v_it; ++v_it) {

    Vec3d v = _mesh.vertex(*v_it).position();
    off << v[0] << " " << v[1] << " " << v[2] << std::endl;
  }

  int n_edges(_mesh.n_edges());
  off << "Edges" << std::endl;
  off << n_edges << std::endl;

  // write edges
  for (typename MeshT::EdgeIter e_it = _mesh.e_iter(); e_it; ++e_it) {

    typename MeshT::VertexHandle from_vertex = _mesh.edge(*e_it).from_vertex();
    typename MeshT::VertexHandle to_vertex = _mesh.edge(*e_it).to_vertex();
    off << from_vertex << " " << to_vertex << std::endl;
  }

  int n_faces(_mesh.n_faces());
  off << "Faces" << std::endl;
  off << n_faces << std::endl;

  // write faces
  for (typename MeshT::FaceIter f_it = _mesh.f_iter(); f_it; ++f_it) {

    off << _mesh.face(*f_it).halfedges().size() << " ";

    std::vector<typename MeshT::HalfEdgeHandle> halfedges = _mesh.face(*f_it).halfedges();

    for (typename std::vector<typename MeshT::HalfEdgeHandle>::const_iterator it = halfedges.begin(); it
        != halfedges.end(); ++it) {

      off << *it;

      if ((it + 1) != halfedges.end())
        off << " ";
    }

    off << std::endl;
  }

  int n_cells(_mesh.n_cells());
  off << "Polyhedra" << std::endl;
  off << n_cells << std::endl;

  for (typename MeshT::CellIter c_it = _mesh.c_iter(); c_it; ++c_it) {

    off << _mesh.cell(*c_it).halffaces().size() << " ";

    std::vector<typename MeshT::HalfFaceHandle> halffaces = _mesh.cell(*c_it).halffaces();

    for (typename std::vector<typename MeshT::HalfFaceHandle>::const_iterator it = halffaces.begin(); it
        != halffaces.end(); ++it) {

      off << *it;

      if ((it + 1) != halffaces.end())
        off << " ";
    }

    off << std::endl;
  }

  off.close();

  return true;
}

//==================================================

bool FileManager::isHexahedralMesh(const std::string& _filename) const {

  std::ifstream iff(_filename.c_str(), std::ios::in);

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
