#ifndef INCLUDE_UNITTESTS_BASICS_HH
#define INCLUDE_UNITTESTS_BASICS_HH

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

/*
 * ====================================================================
 * Define tests below
 * ====================================================================
 */

/*
 */
TEST_F(PolyhedralMeshBase, CreateSimpleMesh) {

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

//===========================================================================

TEST_F(PolyhedralMeshBase, CreateSimpleMeshWithoutCells) {

    Vec3d p1(0.0, 0.0, 0.0);
    Vec3d p2(1.0, 0.0, 0.0);
    Vec3d p3(1.0, 1.0, 0.0);
    Vec3d p4(0.0, 1.0, 0.0);

    Vec3d p5(0.0, 0.0, 1.0);
    Vec3d p6(1.0, 0.0, 1.0);
    Vec3d p7(1.0, 1.0, 1.0);
    Vec3d p8(0.0, 1.0, 1.0);

    VertexHandle v1 = mesh_.add_vertex(p1);
    VertexHandle v2 = mesh_.add_vertex(p2);
    VertexHandle v3 = mesh_.add_vertex(p3);
    VertexHandle v4 = mesh_.add_vertex(p4);

    VertexHandle v5 = mesh_.add_vertex(p5);
    VertexHandle v6 = mesh_.add_vertex(p6);
    VertexHandle v7 = mesh_.add_vertex(p7);
    VertexHandle v8 = mesh_.add_vertex(p8);

    EXPECT_EQ(0, v1);
    EXPECT_EQ(1, v2);
    EXPECT_EQ(2, v3);
    EXPECT_EQ(3, v4);
    EXPECT_EQ(4, v5);
    EXPECT_EQ(5, v6);
    EXPECT_EQ(6, v7);
    EXPECT_EQ(7, v8);

    EdgeHandle e1 = mesh_.add_edge(v1, v2);
    EdgeHandle e2 = mesh_.add_edge(v2, v3);
    EdgeHandle e3 = mesh_.add_edge(v3, v4);
    EdgeHandle e4 = mesh_.add_edge(v4, v1);

    EdgeHandle e5 = mesh_.add_edge(v5, v6);
    EdgeHandle e6 = mesh_.add_edge(v6, v7);
    EdgeHandle e7 = mesh_.add_edge(v7, v8);
    EdgeHandle e8 = mesh_.add_edge(v8, v5);

    EXPECT_EQ(0, e1);
    EXPECT_EQ(1, e2);
    EXPECT_EQ(2, e3);
    EXPECT_EQ(3, e4);
    EXPECT_EQ(4, e5);
    EXPECT_EQ(5, e6);
    EXPECT_EQ(6, e7);
    EXPECT_EQ(7, e8);

    // Get halfedges
    HalfEdgeHandle h1 = mesh_.halfedge_handle(e1, 0u);
    HalfEdgeHandle h2 = mesh_.halfedge_handle(e2, 0u);
    HalfEdgeHandle h3 = mesh_.halfedge_handle(e3, 0u);
    HalfEdgeHandle h4 = mesh_.halfedge_handle(e4, 0u);

    HalfEdgeHandle h5 = mesh_.halfedge_handle(e5, 0u);
    HalfEdgeHandle h6 = mesh_.halfedge_handle(e6, 0u);
    HalfEdgeHandle h7 = mesh_.halfedge_handle(e7, 0u);
    HalfEdgeHandle h8 = mesh_.halfedge_handle(e8, 0u);

    EXPECT_EQ(v1, mesh_.halfedge(h1).from_vertex());
    EXPECT_EQ(v2, mesh_.halfedge(h1).to_vertex());
    EXPECT_EQ(v2, mesh_.halfedge(h2).from_vertex());
    EXPECT_EQ(v3, mesh_.halfedge(h2).to_vertex());
    EXPECT_EQ(v3, mesh_.halfedge(h3).from_vertex());
    EXPECT_EQ(v4, mesh_.halfedge(h3).to_vertex());
    EXPECT_EQ(v4, mesh_.halfedge(h4).from_vertex());
    EXPECT_EQ(v1, mesh_.halfedge(h4).to_vertex());

    EXPECT_EQ(v5, mesh_.halfedge(h5).from_vertex());
    EXPECT_EQ(v6, mesh_.halfedge(h5).to_vertex());
    EXPECT_EQ(v6, mesh_.halfedge(h6).from_vertex());
    EXPECT_EQ(v7, mesh_.halfedge(h6).to_vertex());
    EXPECT_EQ(v7, mesh_.halfedge(h7).from_vertex());
    EXPECT_EQ(v8, mesh_.halfedge(h7).to_vertex());
    EXPECT_EQ(v8, mesh_.halfedge(h8).from_vertex());
    EXPECT_EQ(v5, mesh_.halfedge(h8).to_vertex());

    // Check opposite halfedges
    EXPECT_EQ(v2, mesh_.opposite_halfedge(h1).from_vertex());
    EXPECT_EQ(v1, mesh_.opposite_halfedge(h1).to_vertex());
    EXPECT_EQ(v3, mesh_.opposite_halfedge(h2).from_vertex());
    EXPECT_EQ(v2, mesh_.opposite_halfedge(h2).to_vertex());
    EXPECT_EQ(v4, mesh_.opposite_halfedge(h3).from_vertex());
    EXPECT_EQ(v3, mesh_.opposite_halfedge(h3).to_vertex());
    EXPECT_EQ(v1, mesh_.opposite_halfedge(h4).from_vertex());
    EXPECT_EQ(v4, mesh_.opposite_halfedge(h4).to_vertex());

    EXPECT_EQ(v6, mesh_.opposite_halfedge(h5).from_vertex());
    EXPECT_EQ(v5, mesh_.opposite_halfedge(h5).to_vertex());
    EXPECT_EQ(v7, mesh_.opposite_halfedge(h6).from_vertex());
    EXPECT_EQ(v6, mesh_.opposite_halfedge(h6).to_vertex());
    EXPECT_EQ(v8, mesh_.opposite_halfedge(h7).from_vertex());
    EXPECT_EQ(v7, mesh_.opposite_halfedge(h7).to_vertex());
    EXPECT_EQ(v5, mesh_.opposite_halfedge(h8).from_vertex());
    EXPECT_EQ(v8, mesh_.opposite_halfedge(h8).to_vertex());

    // Add a face via vertices
    std::vector<VertexHandle> vertices;
    vertices.push_back(v2); vertices.push_back(v6);
    vertices.push_back(v7); vertices.push_back(v3);
    FaceHandle f1 = mesh_.add_face(vertices);

    EXPECT_EQ(0, f1);

    // Get halfedges of face
    std::vector<HalfEdgeHandle> halfedges = mesh_.face(f1).halfedges();

    std::vector<HalfEdgeHandle>::iterator it = halfedges.begin();

    EXPECT_EQ(8, mesh_.edge_handle(*it)); ++it;
    EXPECT_EQ(5, mesh_.edge_handle(*it)); ++it;
    EXPECT_EQ(9, mesh_.edge_handle(*it)); ++it;
    EXPECT_EQ(1, mesh_.edge_handle(*it));

    // Add invalid face
    halfedges.clear();
    halfedges.push_back(mesh_.halfedge_handle(e1, 0)); halfedges.push_back(mesh_.halfedge_handle(e2, 0));
    halfedges.push_back(mesh_.halfedge_handle(e7, 0)); halfedges.push_back(mesh_.halfedge_handle(e4, 0));

    FaceHandle fI = mesh_.add_face(halfedges);

    EXPECT_EQ(PolyhedralMesh::InvalidFaceHandle, fI);

    // Now add valid face via edges
    halfedges.clear();
    halfedges.push_back(mesh_.halfedge_handle(e1, 0)); halfedges.push_back(mesh_.halfedge_handle(e2, 0));
    halfedges.push_back(mesh_.halfedge_handle(e3, 0)); halfedges.push_back(mesh_.halfedge_handle(e4, 0));

    FaceHandle f2 = mesh_.add_face(halfedges);

    EXPECT_EQ(1, f2);

    // Get halfedges of face
    halfedges = mesh_.face(f2).halfedges();
    int handle = 0;
    for(it = halfedges.begin(); it != halfedges.end(); ++it) {
        EXPECT_EQ(handle, mesh_.edge_handle(*it)); handle++;
    }
}

TEST_F(PolyhedralMeshBase, VolumeMeshConnectivity) {

    generatePolyhedralMesh(mesh_);

    // Add invalid cell
    std::vector<HalfFaceHandle> hfaces;
    hfaces.push_back(1); hfaces.push_back(5);
    hfaces.push_back(7); hfaces.push_back(9);
    hfaces.push_back(10); hfaces.push_back(21);
    CellHandle i_cell = mesh_.add_cell(hfaces);

    EXPECT_EQ(PolyhedralMesh::InvalidCellHandle, i_cell);

    mesh_.update_adjacencies();

    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(1));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(2));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(5));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(7));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(9));
    EXPECT_EQ(CellHandle(0), mesh_.incident_cell(10));

    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(3));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(12));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(15));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(17));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(19));
    EXPECT_EQ(CellHandle(1), mesh_.incident_cell(20));

    // Test adjacency function
    HalfFaceHandle ad_hf1 = mesh_.adjacent_halfface_in_cell(1, 3);
    // Should be halfface 5
    EXPECT_EQ(HalfFaceHandle(5), ad_hf1);

    HalfFaceHandle ad_hf2 = mesh_.adjacent_halfface_in_cell(1, 7);
    // Should be halfface 7
    EXPECT_EQ(HalfFaceHandle(7), ad_hf2);

    HalfFaceHandle ad_hf3 = mesh_.adjacent_halfface_in_cell(5, 24);
    // Should be invalid
    EXPECT_EQ(PolyhedralMesh::InvalidHalfFaceHandle, ad_hf3);

    HalfFaceHandle ad_hf4 = mesh_.adjacent_halfface_in_cell(12, 24);
    // Should be invalid
    EXPECT_EQ(HalfFaceHandle(20), ad_hf4);

    HalfFaceHandle ad_hf5 = mesh_.adjacent_halfface_in_cell(0, 0);
    // Should be invalid
    EXPECT_EQ(PolyhedralMesh::InvalidHalfFaceHandle, ad_hf5);

    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(2u,  mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_faces());
}

TEST_F(PolyhedralMeshBase, VolumeMeshNormals) {

    generatePolyhedralMesh(mesh_);

    // Request normals
    mesh_.request_face_normals();

    Vec3d n_x(1.0, 0.0, 0.0);
    Vec3d n_y(0.0, 1.0, 0.0);
    Vec3d n_z(0.0, 0.0, 1.0);

    // Should be negative z-axis
    Vec3d n = mesh_.face_normal(0);
    EXPECT_DOUBLE_EQ(-n_z[0], n[0]);
    EXPECT_DOUBLE_EQ(-n_z[1], n[1]);
    EXPECT_DOUBLE_EQ(-n_z[2], n[2]);

    // Should be positive x-axis
    n = mesh_.face_normal(2);
    EXPECT_DOUBLE_EQ(n_x[0], n[0]);
    EXPECT_DOUBLE_EQ(n_x[1], n[1]);
    EXPECT_DOUBLE_EQ(n_x[2], n[2]);

    // Should be positive y-axis
    n = mesh_.face_normal(4);
    EXPECT_DOUBLE_EQ(n_y[0], n[0]);
    EXPECT_DOUBLE_EQ(n_y[1], n[1]);
    EXPECT_DOUBLE_EQ(n_y[2], n[2]);

    // Should be positive y-axis
    n = mesh_.face_normal(5);
    EXPECT_DOUBLE_EQ(n_y[0], n[0]);
    EXPECT_DOUBLE_EQ(n_y[1], n[1]);
    EXPECT_DOUBLE_EQ(n_y[2], n[2]);
}

TEST_F(PolyhedralMeshBase, PolyhedralMeshStatusTest) {

    generatePolyhedralMesh(mesh_);

    // Request status
    mesh_.request_status();

    // Select a few faces
    mesh_.status(FaceHandle(1)).set_tagged(true);
    mesh_.status(FaceHandle(4)).set_tagged(true);

    mesh_.status(HalfFaceHandle(21)).set_deleted(true);
    mesh_.status(HalfFaceHandle(0)).set_deleted(true);

    mesh_.status(VertexHandle(3)).set_selected(true);
    mesh_.status(VertexHandle(8)).set_selected(true);

    EXPECT_TRUE(mesh_.status(FaceHandle(1)).tagged());
    EXPECT_TRUE(mesh_.status(FaceHandle(4)).tagged());
    EXPECT_FALSE(mesh_.status(FaceHandle(7)).tagged());
    EXPECT_FALSE(mesh_.status(FaceHandle(2)).tagged());

    EXPECT_TRUE(mesh_.status(HalfFaceHandle(21)).deleted());
    EXPECT_TRUE(mesh_.status(HalfFaceHandle(0)).deleted());
    EXPECT_FALSE(mesh_.status(HalfFaceHandle(13)).deleted());
    EXPECT_FALSE(mesh_.status(HalfFaceHandle(20)).deleted());

    EXPECT_TRUE(mesh_.status(VertexHandle(3)).selected());
    EXPECT_TRUE(mesh_.status(VertexHandle(8)).selected());
    EXPECT_FALSE(mesh_.status(VertexHandle(1)).selected());
    EXPECT_FALSE(mesh_.status(VertexHandle(9)).selected());
}

TEST_F(PolyhedralMeshBase, PolyhedralMeshProperties) {

    generatePolyhedralMesh(mesh_);

    // Request properties
    OpenVolumeMesh::VPropHandleT<Vec3d> vp;

    mesh_.add_property(vp, "some_v_prop");
    EXPECT_TRUE(mesh_.get_property_handle(vp, "some_v_prop"));

    for(PolyhedralMesh::VertexIter v_it = mesh_.v_iter(); v_it.valid(); ++v_it) {
        mesh_.property(vp, *v_it) = Vec3d(1.0, 0.0, 0.0);
    }

    for(PolyhedralMesh::VertexIter v_it = mesh_.v_iter(); v_it.valid(); ++v_it) {
        Vec3d t;
        t = mesh_.property(vp, *v_it);
        EXPECT_DOUBLE_EQ(1.0, t[0]);
        EXPECT_DOUBLE_EQ(0.0, t[1]);
        EXPECT_DOUBLE_EQ(0.0, t[2]);
    }

    VertexHandle vh = mesh_.add_vertex(Vec3d(3.0,3.0,3.0));
    mesh_.property(vp, vh) = Vec3d(0.0);
    Vec3d p = mesh_.property(vp, vh);
    EXPECT_DOUBLE_EQ(0.0, p[0]);
    EXPECT_DOUBLE_EQ(0.0, p[1]);
    EXPECT_DOUBLE_EQ(0.0, p[2]);

    // Request edge properties
    OpenVolumeMesh::EPropHandleT<unsigned int> ep;

    mesh_.add_property(ep, "some_e_prop");
    EXPECT_TRUE(mesh_.get_property_handle(ep, "some_e_prop"));

    unsigned int i = 0;
    for(PolyhedralMesh::EdgeIter e_it = mesh_.e_iter(); e_it.valid(); ++e_it) {
        mesh_.property(ep, *e_it) = i++;
    }

    i = 0;
    for(PolyhedralMesh::EdgeIter e_it = mesh_.e_iter(); e_it.valid(); ++e_it) {
        EXPECT_EQ(i++, mesh_.property(ep, *e_it));
    }

    // Request halfface properties
    OpenVolumeMesh::HFPropHandleT<bool> hfp;

    mesh_.add_property(hfp, "some_hfp_prop");
    EXPECT_TRUE(mesh_.get_property_handle(hfp, "some_hfp_prop"));

    bool b = false;
    for(PolyhedralMesh::HalfFaceIter hf_it = mesh_.hf_iter(); hf_it.valid(); ++hf_it) {
        mesh_.property(hfp, *hf_it) = b;
        b = !b;
    }

    b = false;
    for(PolyhedralMesh::HalfFaceIter hf_it = mesh_.hf_iter(); hf_it.valid(); ++hf_it) {
        EXPECT_EQ(b, mesh_.property(hfp, *hf_it));
        b = !b;
    }

    // Request halfface properties
    OpenVolumeMesh::CPropHandleT<std::string> cp;

    mesh_.add_property(cp, "some_c_prop");
    EXPECT_TRUE(mesh_.get_property_handle(cp, "some_c_prop"));

    for(PolyhedralMesh::CellIter c_it = mesh_.c_iter(); c_it.valid(); ++c_it) {
        mesh_.property(cp, *c_it) = std::string("MyTestString");
    }

    for(PolyhedralMesh::CellIter c_it = mesh_.c_iter(); c_it.valid(); ++c_it) {
        EXPECT_EQ(std::string("MyTestString"), mesh_.property(cp, *c_it));
    }
}

TEST_F(PolyhedralMeshBase, STLCompliance) {

    generatePolyhedralMesh(mesh_);

    Print p;
    p.mute(true);
    //std::cerr << "Vertices:" << std::endl;
    std::for_each(mesh_.vertices_begin(), mesh_.vertices_end(), p);
    //std::cerr << "Edges:" << std::endl;
    std::for_each(mesh_.edges_begin(), mesh_.edges_end(), p);
    //std::cerr << "HalfEdges:" << std::endl;
    std::for_each(mesh_.halfedges_begin(), mesh_.halfedges_end(), p);
    //std::cerr << "Faces:" << std::endl;
    std::for_each(mesh_.faces_begin(), mesh_.faces_end(), p);
    //std::cerr << "HalfFaces:" << std::endl;
    std::for_each(mesh_.halffaces_begin(), mesh_.halffaces_end(), p);
    //std::cerr << "Cells:" << std::endl;
    std::for_each(mesh_.cells_begin(), mesh_.cells_end(), p);
}

/*
 * Hexahedral mesh tests
 */

TEST_F(HexahedralMeshBase, SimpleHexMeshNavigation) {

    generateHexahedralMesh(mesh_);

    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(11u, mesh_.n_faces());
    EXPECT_EQ(2u,  mesh_.n_cells());

    EXPECT_EQ(1,  mesh_.xfront_halfface(CellHandle(0)));
    EXPECT_EQ(2, mesh_.xback_halfface(CellHandle(0)));
    EXPECT_EQ(5,  mesh_.yfront_halfface(CellHandle(0)));
    EXPECT_EQ(6,  mesh_.yback_halfface(CellHandle(0)));
    EXPECT_EQ(8,  mesh_.zfront_halfface(CellHandle(0)));
    EXPECT_EQ(11,  mesh_.zback_halfface(CellHandle(0)));

    EXPECT_EQ(12, mesh_.opposite_halfface_handle_in_cell(
            HalfFaceHandle(3), CellHandle(1)));

    mesh_.update_adjacencies();

    EXPECT_EQ(HalfFaceHandle(20), mesh_.adjacent_halfface_on_sheet(
            HalfFaceHandle(9), HalfEdgeHandle(12)));
    EXPECT_EQ(HalfFaceHandle(21), mesh_.adjacent_halfface_on_sheet(
            HalfFaceHandle(8), HalfEdgeHandle(12)));

    HexahedralMesh::CellSheetCellIter csc_it = mesh_.csc_iter(0, HexahedralMesh::YF);
    EXPECT_EQ(CellHandle(1), *csc_it);

    HexahedralMesh::HalfFaceSheetHalfFaceIter hfshf_it = mesh_.hfshf_iter(5);
    EXPECT_EQ(HalfFaceHandle(15), *hfshf_it);
    hfshf_it = mesh_.hfshf_iter(6);
    EXPECT_EQ(HalfFaceHandle(16), *hfshf_it);
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest1) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.delete_vertex(VertexHandle(0));

    mesh_.garbage_collection();

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_vertices());
    EXPECT_EQ(17u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest2) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.delete_vertex(VertexHandle(0));

    mesh_.garbage_collection();

    mesh_.garbage_collection();

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_vertices());
    EXPECT_EQ(17u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest3) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.delete_edge(EdgeHandle(0));

    mesh_.garbage_collection(false);

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(19u, mesh_.n_edges());
    EXPECT_EQ(9u, mesh_.n_faces());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTest4) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.delete_edge(EdgeHandle(5));

    mesh_.garbage_collection(false);

    EXPECT_EQ(0u, mesh_.n_cells());
    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(19u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness1) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.update_adjacencies();

    mesh_.delete_vertex(VertexHandle(0));

    mesh_.garbage_collection();

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness2) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.update_adjacencies();

    mesh_.delete_edge(EdgeHandle(0));

    mesh_.garbage_collection();

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness3) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.update_adjacencies();

    mesh_.delete_face(FaceHandle(0));

    mesh_.garbage_collection();

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness4) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.update_adjacencies();

    mesh_.delete_cell(CellHandle(0));

    mesh_.garbage_collection();

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(6u, mesh_.n_faces());
    EXPECT_EQ(12u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestManifoldness5) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    mesh_.update_adjacencies();

    mesh_.delete_edge(EdgeHandle(5));

    mesh_.garbage_collection();

    EXPECT_EQ(0u, mesh_.n_cells());
    EXPECT_EQ(0u, mesh_.n_faces());
    EXPECT_EQ(0u, mesh_.n_edges());
    EXPECT_EQ(0u, mesh_.n_vertices());
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestProps1) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    // Set some properties
    OpenVolumeMesh::FPropHandleT<int> face_prop;
    mesh_.add_property(face_prop, "face_prop");

    mesh_.property(face_prop, FaceHandle(0)) = 11;
    mesh_.property(face_prop, FaceHandle(1)) = 10;
    mesh_.property(face_prop, FaceHandle(2)) = 9;
    mesh_.property(face_prop, FaceHandle(3)) = 8;
    mesh_.property(face_prop, FaceHandle(4)) = 7;
    mesh_.property(face_prop, FaceHandle(5)) = 6;
    mesh_.property(face_prop, FaceHandle(6)) = 5;
    mesh_.property(face_prop, FaceHandle(7)) = 4;
    mesh_.property(face_prop, FaceHandle(8)) = 3;
    mesh_.property(face_prop, FaceHandle(9)) = 2;
    mesh_.property(face_prop, FaceHandle(10)) = 1;

    mesh_.delete_vertex(VertexHandle(0));

    mesh_.garbage_collection();

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(11u, mesh_.n_vertices());
    EXPECT_EQ(17u, mesh_.n_edges());
    EXPECT_EQ(8u, mesh_.n_faces());

    std::set<int> fprops;
    for(HexahedralMesh::FaceIter f_it = mesh_.f_iter(); f_it.valid(); ++f_it) {
        fprops.insert(mesh_.property(face_prop, *f_it));
    }

    EXPECT_EQ(0, fprops.count(11));
    EXPECT_EQ(1, fprops.count(10));
    EXPECT_EQ(1, fprops.count(9));
    EXPECT_EQ(0, fprops.count(8));
    EXPECT_EQ(0, fprops.count(7));
    EXPECT_EQ(1, fprops.count(6));
    EXPECT_EQ(1, fprops.count(5));
    EXPECT_EQ(1, fprops.count(4));
    EXPECT_EQ(1, fprops.count(3));
    EXPECT_EQ(1, fprops.count(2));
    EXPECT_EQ(1, fprops.count(1));
}

TEST_F(HexahedralMeshBase, GarbageCollectionTestProps2) {

    generateHexahedralMesh(mesh_);

    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_cell_status();

    // Set some properties
    OpenVolumeMesh::FPropHandleT<int> face_prop;
    mesh_.add_property(face_prop, "face_prop");

    mesh_.property(face_prop, FaceHandle(0)) = 11;
    mesh_.property(face_prop, FaceHandle(1)) = 10;
    mesh_.property(face_prop, FaceHandle(2)) = 9;
    mesh_.property(face_prop, FaceHandle(3)) = 8;
    mesh_.property(face_prop, FaceHandle(4)) = 7;
    mesh_.property(face_prop, FaceHandle(5)) = 6;
    mesh_.property(face_prop, FaceHandle(6)) = 5;
    mesh_.property(face_prop, FaceHandle(7)) = 4;
    mesh_.property(face_prop, FaceHandle(8)) = 3;
    mesh_.property(face_prop, FaceHandle(9)) = 2;
    mesh_.property(face_prop, FaceHandle(10)) = 1;

    mesh_.delete_face(FaceHandle(0));

    mesh_.garbage_collection();

    EXPECT_EQ(1u, mesh_.n_cells());
    EXPECT_EQ(12u, mesh_.n_vertices());
    EXPECT_EQ(20u, mesh_.n_edges());
    EXPECT_EQ(10u, mesh_.n_faces());

    std::set<int> fprops;
    for(HexahedralMesh::FaceIter f_it = mesh_.f_iter(); f_it.valid(); ++f_it) {
        fprops.insert(mesh_.property(face_prop, *f_it));
    }

    EXPECT_EQ(0, fprops.count(11));
    EXPECT_EQ(1, fprops.count(10));
    EXPECT_EQ(1, fprops.count(9));
    EXPECT_EQ(1, fprops.count(8));
    EXPECT_EQ(1, fprops.count(7));
    EXPECT_EQ(1, fprops.count(6));
    EXPECT_EQ(1, fprops.count(5));
    EXPECT_EQ(1, fprops.count(4));
    EXPECT_EQ(1, fprops.count(3));
    EXPECT_EQ(1, fprops.count(2));
    EXPECT_EQ(1, fprops.count(1));
}

#endif // INCLUDE GUARD
