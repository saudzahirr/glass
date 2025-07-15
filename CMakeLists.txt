cmake_minimum_required(VERSION 3.22)
project(glass VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
enable_language(Fortran)
enable_language(C)

file(GLOB CMAKE_FILES "${CMAKE_SOURCE_DIR}/cmake/*.cmake")

foreach(CMAKE_FILE ${CMAKE_FILES})
    include(${CMAKE_FILE})
endforeach()

SYSTEM_INFO()

add_subdirectory(src)
add_subdirectory(cglass)
add_subdirectory(tests)
