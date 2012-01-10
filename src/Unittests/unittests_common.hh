#ifndef INCLUDE_UNITTESTS_COMMON_HH
#define INCLUDE_UNITTESTS_COMMON_HH

#include <gtest/gtest.h>

#include <OpenVolumeMesh/PolyhedralMesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/HexahedralMesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

// Vector defines
typedef OpenVolumeMesh::Vec3d Vec3d;
typedef OpenVolumeMesh::Vec3f Vec3f;

/*
 * Simple test setting for polyhedral meshes
 */

typedef OpenVolumeMesh::PolyhedralMesh<OpenVolumeMesh::Vec3d> PolyhedralMesh;

class PolyhedralMeshBase: public testing::Test {

protected:

  typedef PolyhedralMesh::VertexHandle    VertexHandle;
  typedef PolyhedralMesh::HalfEdgeHandle  HalfEdgeHandle;
  typedef PolyhedralMesh::EdgeHandle      EdgeHandle;
  typedef PolyhedralMesh::HalfFaceHandle  HalfFaceHandle;
  typedef PolyhedralMesh::FaceHandle      FaceHandle;
  typedef PolyhedralMesh::CellHandle      CellHandle;

  // This function is called before each test is run
  virtual void SetUp()
  {

    // Do some initial stuff with the member data here...
  }

  // This function is called after all tests are through
  virtual void TearDown()
  {

    // Do some final stuff with the member data here...
  }

  // Generate a basic hexahedral mesh
  void generatePolyhedralMesh(PolyhedralMesh& _mesh);

  // This member will be accessible in all tests
  PolyhedralMesh mesh_;
};

/*
 * Simple test setting for hexahedral meshes
 */

typedef OpenVolumeMesh::HexahedralMesh<OpenVolumeMesh::Vec3d> HexahedralMesh;

class HexahedralMeshBase: public testing::Test {

protected:

  typedef HexahedralMesh::VertexHandle    VertexHandle;
  typedef HexahedralMesh::HalfEdgeHandle  HalfEdgeHandle;
  typedef HexahedralMesh::EdgeHandle      EdgeHandle;
  typedef HexahedralMesh::HalfFaceHandle  HalfFaceHandle;
  typedef HexahedralMesh::FaceHandle      FaceHandle;
  typedef HexahedralMesh::CellHandle      CellHandle;

  // This function is called before each test is run
  virtual void SetUp()
  {

    // Do some initial stuff with the member data here...
  }

  // This function is called after all tests are through
  virtual void TearDown()
  {

    // Do some final stuff with the member data here...
  }

  // Generate a basic hexahedral mesh
  void generateHexahedralMesh(HexahedralMesh& _mesh);

  // This member will be accessible in all tests
  HexahedralMesh mesh_;
};

// Printer class (for STL compliance test)
class Print {
public:
  Print(bool _mute = false) : mute_(_mute) {}
  void mute(bool _mute) { mute_ = _mute; }
  bool mute() const { return mute_; }
  void operator()(const OpenVolumeMesh::OpenVolumeMeshHandle& _h) const {
    if(!mute_) std::cerr << "Handle: " << _h.idx() << std::endl;
  }
private:
  bool mute_;
};

#endif // INCLUDE GUARD
