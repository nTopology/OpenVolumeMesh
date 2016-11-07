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


#ifndef ACG_UTILS_SMARTPOINTER_HH
#define ACG_UTILS_SMARTPOINTER_HH

/**********************************************
 * Warning! This header file is duplicated in *
 * OpenFlipper with the same header guard, as *
 * ACG/Utils/SmartPointer.hh.                 *
 * If you change this file, you should change *
 * that file as well.                         *
 **********************************************/

#include <memory>

// legacy code may depend on this define:
#define ACG_UNIQUE_POINTER_SUPPORTED 1

namespace ptr {
    using std::shared_ptr;
    using std::make_shared;
    using std::unique_ptr;
#if __cplusplus >= 201402L
    using std::make_unique;
#else
    template<typename T, typename... Args>
    std::unique_ptr<T>
    make_unique(Args&&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }
#endif // C++14
} // namespace ptr

#endif // ACG_UTILS_SMARTPOINTER_HH
