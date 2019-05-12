cmake_minimum_required (VERSION 2.8)
cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 11)

add_definitions(--std=c++11 )

project (WaveSolver CXX)

set( LEARN_CMAKE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
set( LEARN_CMAKE_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} )

if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/data/")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data")
    message("Created directory ${CMAKE_CURRENT_BINARY_DIR}/data")
else()
    # file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/data/") #doesn't work somehow
    file(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/data/*")
    foreach(file ${files})
        file(REMOVE "${file}")
    endforeach()
    message("Cleaned previous output from directory ${CMAKE_CURRENT_BINARY_DIR}/data/*")
endif()

if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/graph/")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/graph")
    message("Created directory ${CMAKE_CURRENT_BINARY_DIR}/graph")
else()
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/graph/*")
    # file(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/graph/*") #doesn't work somehow
    file(GLOB files "${CMAKE_CURRENT_BINARY_DIR}/graph/*")
    foreach(file ${files})
        file(REMOVE "${file}")
    endforeach()
    message("Cleaned previous output from directory ${CMAKE_CURRENT_BINARY_DIR}/graph/*")
endif()

include_directories ("${PROJECT_SOURCE_DIR}/../")
include_directories ("${PROJECT_SOURCE_DIR}/Catch2")
include_directories ("${PROJECT_SOURCE_DIR}/Wave1DSolver")
include_directories ("${PROJECT_SOURCE_DIR}/Guitar")
include_directories ("${PROJECT_SOURCE_DIR}/Impulse")
include_directories ("${PROJECT_SOURCE_DIR}/Test")
include_directories ("${PROJECT_SOURCE_DIR}/TestQuadratic")
include_directories ("${PROJECT_SOURCE_DIR}/TestConstant")


add_library(Wave1DSolver  Wave1DSolver/Wave1DSolver.cpp)
add_library(Guitar        Guitar/Guitar.cpp)
add_library(Impulse       Impulse/Impulse.cpp)
add_library(Test          Test/Test.cpp)
add_library(TestQuadratic TestQuadratic/TestQuadratic.cpp)
add_library(TestConstant  TestConstant/TestConstant.cpp)

#set (LEARN_CMAKE_VERSION_MAJOR 1)
#set (LEARN_CMAKE_VERSION_MINOR 0)
#set (LEARN_CMAKE_VERSION_PATCH 0)

add_executable(WaveSolver main.cpp)
target_link_libraries (WaveSolver Wave1DSolver Guitar Test TestQuadratic TestConstant Impulse)
