#
# Top-level CMakeLists.txt for ZSIG library examples
#
project(zsig)

cmake_minimum_required(VERSION 2.4)
cmake_policy(VERSION 2.4)
set(CMAKE_CXX_COMPILER g++)

message(STATUS "Building in ${BUILDENV} Environment")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_BUILD_TYPE Release)

if( CMAKE_SIZEOF_VOID_P MATCHES 4 ) # 32-bit O.S.
  set(BASE_FLAGS "-Wno-deprecated -march=native -mfpmath=sse -fopenmp")
else( CMAKE_SIZEOF_VOID_P MATCHES 4 ) # 64-bit O.S.
  set(BASE_FLAGS "-Wno-deprecated -m64 -march=native -mfpmath=sse -fopenmp")
endif( CMAKE_SIZEOF_VOID_P MATCHES 4 )

set(RELEASE_FLAGS "-O3 -DNDEBUG -ffast-math -Wall -Wno-parentheses -funroll-loops -frounding-math -fomit-frame-pointer -ftree-vectorize -fno-math-errno -funsafe-math-optimizations -Wall -Wextra -Wno-parentheses -Wno-unused -Wno-sign-compare -Wno-unknown-pragmas")
set(DEBUG_FLAGS "-O0 -ggdb")

set(CMAKE_CXX_FLAGS "${BASE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${BASE_FLAGS} ${RELEASE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${BASE_FLAGS} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${BASE_FLAGS} ${DEBUG_FLAGS} ${RELEASE_FLAGS}")

# Set this variable to true if you want to compile examples using OpenGL + GLUT
set(ENABLE_GL TRUE)

if(APPLE AND ENABLE_GL)
  set(CMAKE_CXX_LINK_FLAGS "-framework OpenGL -framework GLUT")
endif(APPLE AND ENABLE_GL)

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

include_directories(${INCLUDE_DIR})

add_subdirectory(src)
