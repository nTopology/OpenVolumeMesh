/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
\*===========================================================================*/

#ifndef ACG_UTILS_SMARTPOINTER_HH
#define ACG_UTILS_SMARTPOINTER_HH

/*************************************************
 * Warning! This header file is duplicated in    *
 * OpenVolumeMesh with the same header guard,    *
 * as src/OpenVolumeMesh/System/MemoryInclude.hh.*
 * If you change this file, you should change    *
 * that file as well.                            *
 *************************************************/

#include <memory>

// Default for C++ 11 and higher
#if ((defined(_MSC_VER) && (_MSC_VER >= 1900)) || __cplusplus > 199711L || defined(__GXX_EXPERIMENTAL_CXX0X__))

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


#else // deprecated things

    /** This set of defines maps the pointer namespaces to the namespace ptr depending on
     *  the current architecture and compilers.
     */
    #if ( (__cplusplus >= 201103L)  || (__STDC_VERSION__ >= 201112L) )
       // C++11:
       #include <memory>
       namespace ptr = std;
       #define ACG_UNIQUE_POINTER_SUPPORTED 1
    #elif defined(__GXX_EXPERIMENTAL_CXX0X__)
       // C++11 via -std=c++0x on gcc:
       #include <memory>
       namespace ptr = std;
       #define ACG_UNIQUE_POINTER_SUPPORTED 1
    #else
       // C++98 and TR1:
       #if (_MSC_VER >= 1600)
         // VStudio 2010 supports some C++11 features
         #include <memory>
         namespace ptr = std;
         #define ACG_UNIQUE_POINTER_SUPPORTED 1
       #elif (_MSC_VER >= 1500)
         // hope for TR1 equivalents
        #if(_HAS_TR1)
         #include <memory>
         namespace ptr = std::tr1;
         #define ACG_UNIQUE_POINTER_SUPPORTED 0
        #else
         #pragma warning "TR1 not available! Please install Visual Studio Service Pack 1!"
        #endif
       #else
         // hope for TR1 equivalents
         // check for clang5 but switch to tr1 if clang uses libstdc++
         #if defined(__clang_major__) && (__clang_major__ >= 5) && !defined(__GLIBCXX__ )
          // Mavericks special treatment
          #include <memory>
          namespace ptr = std;
        #else
          #include <tr1/memory>
          namespace ptr = std::tr1;
        #endif
        #define ACG_UNIQUE_POINTER_SUPPORTED 0
       #endif
    #endif
#endif




#endif // ACG_UTILS_SMARTPOINTER_HH

