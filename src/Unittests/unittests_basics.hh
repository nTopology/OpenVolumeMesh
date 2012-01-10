#ifndef INCLUDE_UNITTESTS_DECIMATER_HH
#define INCLUDE_UNITTESTS_DECIMATER_HH

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

class SimplePolyhedralMesh : public PolyhedralMeshBase {

protected:

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
};

/*
 * ====================================================================
 * Define tests below
 * ====================================================================
 */

/*
 */
TEST_F(SimplePolyhedralMesh, CreateSimpleMesh) {

  /*
   * Add vertices
   */

  VertexHandle v0 = mesh_.add_vertex(Vec3d(-1.0, -1.0, -1.0));
  VertexHandle v1 = mesh_.add_vertex(Vec3d( 1.0, -1.0, -1.0));
  VertexHandle v2 = mesh_.add_vertex(Vec3d( 1.0,  1.0, -1.0));
  VertexHandle v3 = mesh_.add_vertex(Vec3d(-1.0,  1.0, -1.0));
  VertexHandle v4 = mesh_.add_vertex(Vec3d(-1.0, -1.0,  1.0));
  VertexHandle v5 = mesh_.add_vertex(Vec3d( 1.0, -1.0,  1.0));
  VertexHandle v6 = mesh_.add_vertex(Vec3d( 1.0,  1.0,  1.0));
  VertexHandle v7 = mesh_.add_vertex(Vec3d(-1.0,  1.0,  1.0));

  EXPECT_EQ(8u, mesh_.n_vertices()) << "The number of vertices is not correct!";

  /*
   * Add faces
   */

  std::vector<VertexHandle> fvertices;

  fvertices.push_back(v3);
  fvertices.push_back(v2);
  fvertices.push_back(v1);
  fvertices.push_back(v0);

  FaceHandle fh0 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v4);
  fvertices.push_back(v5);
  fvertices.push_back(v6);
  fvertices.push_back(v7);

  FaceHandle fh1 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v0);
  fvertices.push_back(v4);
  fvertices.push_back(v7);
  fvertices.push_back(v3);

  FaceHandle fh2 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v1);
  fvertices.push_back(v2);
  fvertices.push_back(v6);
  fvertices.push_back(v5);

  FaceHandle fh3 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v7);
  fvertices.push_back(v6);
  fvertices.push_back(v2);
  fvertices.push_back(v3);

  FaceHandle fh4 = mesh_.add_face(fvertices);

  fvertices.clear();

  fvertices.push_back(v0);
  fvertices.push_back(v1);
  fvertices.push_back(v5);
  fvertices.push_back(v4);

  FaceHandle fh5 = mesh_.add_face(fvertices);

  EXPECT_EQ(12u, mesh_.n_edges()) << "The number of edges is not correct!";
  EXPECT_EQ(6u, mesh_.n_faces())  << "The number of faces is not correct!";

  /*
   * Add cell
   */

  std::vector<HalfFaceHandle> chfaces;

  chfaces.push_back(mesh_.halfface_handle(fh0, 0));
  chfaces.push_back(mesh_.halfface_handle(fh1, 0));
  chfaces.push_back(mesh_.halfface_handle(fh2, 0));
  chfaces.push_back(mesh_.halfface_handle(fh3, 0));
  chfaces.push_back(mesh_.halfface_handle(fh4, 0));
  chfaces.push_back(mesh_.halfface_handle(fh5, 0));

  mesh_.add_cell(chfaces);

  EXPECT_EQ(1u, mesh_.n_cells())  << "The number of cells is not correct!";
}

#endif // INCLUDE GUARD
