project( TetraMaker )
cmake_minimum_required(VERSION 2.8)

include( Common.cmake )

#set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH} )

set( TetraMaker_HEADERS
    GetPot
   )

set( TetraMaker_INCLUDE_DIRS
        ${CMAKE_CURRENT_SOURCE_DIR}
   )

include_directories( ${TetraMaker_INCLUDE_DIRS} )

add_executable( tetramaker main.cpp ${TetraMaker_HEADERS} )
