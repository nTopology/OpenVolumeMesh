#ifndef INCLUDE_UNITTESTS_FILES_HH
#define INCLUDE_UNITTESTS_FILES_HH

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

TEST_F(PolyhedralMeshBase, LoadFile) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());
}

TEST_F(HexahedralMeshBase, LoadFile) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.isHexahedralMesh("Cylinder.ovm"));
  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());
}

TEST_F(PolyhedralMeshBase, LoadFileWithProps) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cube_with_props.ovm", mesh_));

  EXPECT_EQ(8u, mesh_.n_vertices());
  EXPECT_EQ(12u, mesh_.n_edges());
  EXPECT_EQ(6u, mesh_.n_faces());
  EXPECT_EQ(1u, mesh_.n_cells());
}

#endif // INCLUDE GUARD
