#ifndef UNITTESTS_BASICS_HH_
#define UNITTESTS_BASICS_HH_

#include <iostream>

#include <gtest/gtest.h>

#include "unittests_common.hh"

#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>

using namespace OpenVolumeMesh;

TEST_F(PolyhedralMeshBase, PropertySmartPointerTest1) {

    generatePolyhedralMesh(mesh_);

    EXPECT_EQ(12u, mesh_.n_vertices());

    VertexPropertyT<double> v_prop_d = mesh_.request_vertex_property<double>("MyVPropDbl");

    for(int i = 0; i < 1; ++i) {

        VertexPropertyT<float> v_prop = mesh_.request_vertex_property<float>("MyVProp");

        v_prop[0] = 1.4f;

        VertexPropertyT<float> v_prop2(v_prop);

        VertexPropertyT<float> v_prop3 = v_prop;

        EXPECT_EQ(12u, v_prop3->n_elements());

        EXPECT_EQ(2u, mesh_.n_vertex_props());

        VertexPropertyT<float> v_prop_duplicate = mesh_.request_vertex_property<float>("MyVProp");

        EXPECT_EQ(3u, mesh_.n_vertex_props());

        EXPECT_FLOAT_EQ(1.4f, v_prop3[0]);

        VertexPropertyT<std::string> v_prop_duplicate_2 = mesh_.request_vertex_property<std::string>("MyVProp");

        EXPECT_EQ(4u, mesh_.n_vertex_props());
    }

    v_prop_d.set_persistent();

    EXPECT_EQ(1u, mesh_.n_vertex_props());

    mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));

    EXPECT_EQ(13u, mesh_.n_vertices());

    VertexPropertyT<double> v_prop_d2 = mesh_.request_vertex_property<double>("MyVPropDbl");

    EXPECT_EQ(1u, mesh_.n_vertex_props());

    HalfEdgePropertyT<int> he_prop = mesh_.request_halfedge_property<int>("MyHEProp");

    EXPECT_EQ(40u, he_prop->n_elements());

    mesh_.add_edge(VertexHandle(0), VertexHandle(2));

    EXPECT_EQ(42u, he_prop->n_elements());
}

TEST_F(HexahedralMeshBase, PropertySmartPointerPersistencyTest1) {

    generateHexahedralMesh(mesh_);

    for(int i = 0; i < 1; ++i) {

        VertexPropertyT<float> v_prop = mesh_.request_vertex_property<float>("FloatVProp");

        v_prop[0] = 24.5f;
        v_prop[11] = 2.34f;

        v_prop.set_persistent();
    }

    VertexPropertyT<float> v_prop2 = mesh_.request_vertex_property<float>("FloatVProp");

    EXPECT_FLOAT_EQ(24.5f, v_prop2[0]);
    EXPECT_FLOAT_EQ(2.34f, v_prop2[11]);
}

TEST_F(PolyhedralMeshBase, StatusTest) {

    generatePolyhedralMesh(mesh_);

    EXPECT_EQ(12u, mesh_.n_vertices());

    StatusAttrib status(mesh_);
}

#endif /* UNITTESTS_BASICS_HH_ */
