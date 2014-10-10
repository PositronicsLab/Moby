cmake_minimum_required(VERSION 2.8.3)

project(test)

include (CheckIncludeFiles)
include (CheckLibraryExists)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)
set (PROJECT_SRC_DIR .)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
find_package(Boost REQUIRED)

# you may have to add additional includes in here...
add_definitions(-DBUILD_DOUBLE)
add_definitions(-DUSE_OSG)
add_definitions(-DSAFESTATIC=static)
add_definitions(-DUSE_GLPK)

include_directories(
                    .
                    ../../include
                    /usr/include
                    /usr/include/libxml2
                    /usr/local/include
                    ${Boost_INCLUDE_DIR}
                    /opt/local/include
                    /opt/X11/include
                    )

link_directories(
                  /usr/local/lib
                  /opt/X11/lib
                  )

add_library(${PROJECT_NAME}Plugin MODULE feeder.cpp)
target_link_libraries(${PROJECT_NAME}Plugin Moby Ravelin)
