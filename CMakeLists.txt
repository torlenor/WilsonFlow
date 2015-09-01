cmake_minimum_required(VERSION 2.8)

project( FLOW_PROJECT ) # this is case sensitive 

######################################

SET(CODE_FILES  flow.cpp init.cpp options.cpp matrixhandling.cpp su3.cc)

set(CMAKE_BUILD_TYPE Release)

# SET(GCC_ADDITIONAL_COMPILE_FLAGS  "-O2 -Wall -std=c++11 -DDEBUG")
SET(GCC_ADDITIONAL_COMPILE_FLAGS  "-O3 -Wall -DNDEBUG -fopenmp")
SET(GCC_ADDITIONAL_LINK_FLAGS     "")
SET(GCC_ADDITIONAL_LIBRARIES      "")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_ADDITIONAL_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_ADDITIONAL_LINK_FLAGS}" )

add_executable(flow ${CODE_FILES})

target_link_libraries(flow ${GCC_ADDITIONAL_LIBRARIES})
