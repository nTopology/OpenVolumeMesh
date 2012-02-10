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

#ifndef STATUSATTRIB_HH_
#define STATUSATTRIB_HH_

#include <cassert>

#include "../Core/OpenVolumeMeshHandle.hh"
#include "OpenVolumeMeshStatus.hh"
#include "../Core/PropertyDefines.hh"

namespace OpenVolumeMesh {

// Forward declaration
class TopologyKernel;

class StatusAttrib {
public:
    explicit StatusAttrib(TopologyKernel& _kernel);
    ~StatusAttrib();

    const OpenVolumeMeshStatus& operator[](const VertexHandle& _h) const {
        return v_status_[_h];
    }

    OpenVolumeMeshStatus& operator[](const VertexHandle& _h) {
        return v_status_[_h];
    }

    const OpenVolumeMeshStatus& operator[](const EdgeHandle& _h) const {
        return e_status_[_h];
    }

    OpenVolumeMeshStatus& operator[](const EdgeHandle& _h) {
        return e_status_[_h];
    }

    const OpenVolumeMeshStatus& operator[](const HalfEdgeHandle& _h) const {
        return he_status_[_h];
    }

    OpenVolumeMeshStatus& operator[](const HalfEdgeHandle& _h) {
        return he_status_[_h];
    }

    const OpenVolumeMeshStatus& operator[](const FaceHandle& _h) const {
        return f_status_[_h];
    }

    OpenVolumeMeshStatus& operator[](const FaceHandle& _h) {
        return f_status_[_h];
    }

    const OpenVolumeMeshStatus& operator[](const HalfFaceHandle& _h) const {
        return hf_status_[_h];
    }

    OpenVolumeMeshStatus& operator[](const HalfFaceHandle& _h) {
        return hf_status_[_h];
    }

    const OpenVolumeMeshStatus& operator[](const CellHandle& _h) const {
        return c_status_[_h];
    }

    OpenVolumeMeshStatus& operator[](const CellHandle& _h) {
        return c_status_[_h];
    }

    const OpenVolumeMeshStatus& mesh_status() const {
        OpenVolumeMeshHandle h(0);
        return m_status_[h];
    }

    OpenVolumeMeshStatus& mesh_status() {
        OpenVolumeMeshHandle h(0);
        return m_status_[h];
    }

    /**
     * \brief Delete all entities that have been marked as deleted
     *
     * This function deletes all entities that have been marked as deleted.
     * It proceeds bottom-up, starting with the vertices. All higher
     * dimensional entities that are incident to a deleted entity are
     * automatically marked deleted, too. Once this first pass is through,
     * one can additionally delete all resulting non-manifold configurations
     * in a second pass (triggered by the parameter of this function).
     * This step proceeds as follows: Delete all n-dimensional entities
     * (starting with n = 2), that are not incident to at least one
     * entity of dimension n + 1. Note that the second pass requires bottom-up
     * adjacencies to be available. Compute them by calling update_adjacencies().
     *
     * @param _preserveManifoldness Pass true if the mesh is required to stay three-manifold
     */
    void garbage_collection(bool _preserveManifoldness = false);

private:

    TopologyKernel& kernel_;

    VertexPropertyT<OpenVolumeMeshStatus> v_status_;
    EdgePropertyT<OpenVolumeMeshStatus> e_status_;
    HalfEdgePropertyT<OpenVolumeMeshStatus> he_status_;
    FacePropertyT<OpenVolumeMeshStatus> f_status_;
    HalfFacePropertyT<OpenVolumeMeshStatus> hf_status_;
    CellPropertyT<OpenVolumeMeshStatus> c_status_;
    MeshPropertyT<OpenVolumeMeshStatus> m_status_;
};

} // Namespace OpenVolumeMesh

#endif /* STATUSATTRIB_HH_ */
