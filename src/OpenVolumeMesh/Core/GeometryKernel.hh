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

#ifndef GEOMETRYKERNEL_HH_
#define GEOMETRYKERNEL_HH_

#include <cassert>

#include "../Geometry/VectorT.hh"
#include "TopologyKernel.hh"

namespace OpenVolumeMesh {

template <class VecT>
class GeometryKernel : public TopologyKernel {
public:

    typedef VecT PointT;

    /// Constructor
    GeometryKernel() {}

    /// Destructor
    ~GeometryKernel() {}

    /// Override of empty add_vertex function
    virtual VertexHandle add_vertex() { return add_vertex(VecT()); }

    /// Add a geometric point to the mesh
    VertexHandle add_vertex(const VecT& _p) {

        // Store vertex in list
        vertices_.push_back(_p);

        // Resize vertex props
        resize_vprops(vertices_.size());

        // Get handle of recently created vertex
        return VertexHandle(vertices_.size() - 1);
    }

    /// Set the coordinates of point _vh
    void set_vertex(const VertexHandle& _vh, const VecT& _p) {

        assert(_vh.idx() < (int)vertices_.size());

        vertices_[_vh] = _p;
    }

    /// Get point _vh's coordinates
    const VecT& vertex(const VertexHandle& _vh) const {
        return vertices_[_vh];
    }

    /// Override of n_vertices()
    virtual unsigned int n_vertices() const {
        return vertices_.size();
    }

    virtual VertexIter delete_vertex(const VertexHandle& _h) {
        assert(_h.idx() < (int)n_vertices());

        VertexIter nV = TopologyKernel::delete_vertex(_h);

        vertices_.erase(vertices_.begin() + _h.idx());

        return nV;
    }

    virtual void clear() {

        vertices_.clear();
        TopologyKernel::clear();
    }

private:

    std::vector<VecT> vertices_;
};

} // Namespace OpenVolumeMesh

#endif /* GEOMETRYKERNEL_HH_ */
