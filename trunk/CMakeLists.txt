#
# Top-level CMakeLists.txt for ZSIG library examples
#
project(zsig)

cmake_minimum_required(VERSION 2.6)
cmake_policy(VERSION 2.6)
set(CMAKE_CXX_COMPILER g++)

message(STATUS "Building in ${BUILDENV} Environment")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_BUILD_TYPE Release)

set(BASE_FLAGS "-Wno-deprecated -march=native -mfpmath=sse -fopenmp")

if( CMAKE_SIZEOF_VOID_P MATCHES 8 ) # 64-bit O.S.
  set(BASE_FLAGS "${BASE_FLAGS} -m64")
endif()

set(RELEASE_FLAGS "-O3 -DNDEBUG -ffast-math -Wall -Wno-parentheses -funroll-loops -frounding-math -fomit-frame-pointer -ftree-vectorize -fno-math-errno -funsafe-math-optimizations -Wall -Wextra -Wno-parentheses -Wno-unused -Wno-sign-compare -Wno-unknown-pragmas")
set(DEBUG_FLAGS "-O0 -ggdb")

set(CMAKE_CXX_FLAGS "${BASE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${BASE_FLAGS} ${RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${BASE_FLAGS} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${BASE_FLAGS} ${DEBUG_FLAGS} ${RELEASE_FLAGS}")

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

include_directories(${INCLUDE_DIR})

add_subdirectory(src)
