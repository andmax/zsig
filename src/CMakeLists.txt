#
# Src-level CMakeLists.txt for ZSIG library examples
#
set(EX1_SRCS example_1.cc)
set(EX2_SRCS example_2.cc)

add_executable(ex1 ${EX1_SRCS})
add_executable(ex2 ${EX2_SRCS})

if( ENABLE_GL )
  link_libraries(exGL GL GLU glut)
  add_executable(exGL example_gl.cc)
endif( ENABLE_GL )
