#
# Src-level CMakeLists.txt for ZSIG library examples
#
set(EX1_SRCS example_zpol.cc)
set(EX2_SRCS example_zproj.cc)

add_executable(ex1 ${EX1_SRCS})
add_executable(ex2 ${EX2_SRCS})
