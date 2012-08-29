################################################################
################################################################
#
# To use this file, add "include( Common.cmake )" to the top of your
# CMakeLists.txt
#
################################################################
################################################################


################## AVOID BUILDING IN-SOURCE ####################
# Custom command to cleanup cmake polution when user accidently attempts an
# insoruce build.
#ADD_CUSTOM_TARGET( cleansrc )
#ADD_CUSTOM_COMMAND(
#        TARGET  cleansrc
#        COMMAND rm -rf `find . -name CMakeFiles -or -name Makefile -or -name cmake_install.cmake -or -name CMakeCache.txt`
#        )
## Force good behavior: user shall not try an insource build.
#if( ${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR} )
#    MESSAGE( WARNING "\nERROR: You must build outside of the source tree.\n"
#                         "       Recommend running 'make cleansrc' to clean the source tree.\n" )
#endif()



################### SET BUILD TYPE OPTIONS ######################
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
IF( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
    # Verbose compile when debugging, and lots of optimization otherwise
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -g -pg")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g -pg")
ENDIF( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
 IF( "${CMAKE_BUILD_TYPE}" STREQUAL "Release" )
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -Wall -Wextra -O3 -funroll-loops -finline-functions -mmmx -msse2 ")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -O3 -funroll-loops -finline-functions -mmmx -msse2 ")
ENDIF( "${CMAKE_BUILD_TYPE}" STREQUAL "Release" )

############### Load CPATH and LIBRARY_PATH ##################
string( REPLACE ":" ";" COMMON_INCLUDE_DIRS "$ENV{CPATH}" ) 
string( REPLACE ":" ";" COMMON_LINK_DIRS "$ENV{LIBRARY_PATH}" ) 
list( REMOVE_DUPLICATES COMMON_INCLUDE_DIRS )
list( REMOVE_DUPLICATES COMMON_LINK_DIRS )

include_directories( ${COMMON_INCLUDE_DIRS} )
link_directories( ${COMMON_LINK_DIRS} )

# Hack to get eclipse working
if( APPLE )
    SET( CMAKE_ECLIPSE_EXECUTABLE /Applications/eclipse/eclipse CACHE STRING "" FORCE )
endif()

################################################################
# macros for library writers.
macro( add_to_lib_include_dirs lib dirs )
    list( APPEND ${lib}_INCLUDE_DIRS ${dirs} )
    set( ${lib}_INCLUDE_DIRS ${${lib}_INCLUDE_DIRS} CACHE INTERNAL "" FORCE )
endmacro()

macro( add_to_lib_libraries lib libs )
#    message( STATUS "DEALING with ${libs}" )
    foreach( l ${libs} )
#        message( STATUS "LIB: ${l}" )
        get_target_property( _LIBRARY ${l} LOCATION )
        if( ${_LIBRARY} STREQUAL "_LIBRARY-NOTFOUND" )
#        message( STATUS "LOCATION: ${_LIBRARY}" )
        else()
#            message( STATUS "looking for ${l} got ${_LIBRARY}" )
#            message( STATUS "******* adding '${_LIBRARY}' ")
            list( INSERT ${lib}_LIBRARIES 0 ${_LIBRARY} )
            set( ${lib}_LIBRARIES ${${lib}_LIBRARIES} CACHE INTERNAL "" FORCE )
#            message( STATUS "******* result '${${lib}_LIBRARIES}' ")
        endif()
    endforeach()
endmacro()

macro( add_to_lib_link_directories dirs )
    list( APPEND ${lib}_LINK_DIRECTORIES ${dirs} )
    set( ${lib}_LINK_DIRECTORIES ${${lib}_LINK_DIRECTORIES} CACHE INTERNAL "" FORCE )
endmacro()


ADD_CUSTOM_TARGET( WriteConfig )
ADD_CUSTOM_COMMAND(
        TARGET WriteConfig 
        COMMAND echo test >> test.cmake
        )


macro( export_library LibName )
    #add_to_lib_include_dirs( SimpleGui dirs )
    #add_to_lib_libraries( SimpleGui libs )
#    get_property( dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES )
#    message( STATUS "got ${dirs}" )

#    get_property( link_dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES )
#    message( STATUS "got ${link_dirs}" )
#    message( STATUS "Got IncDirs: ${IncDirs}" )
#    message( STATUS "Got Libs: ${Libs}" )

#    message( STATUS "Want SimpleGui_LIBRARIES: ${${LibName}_LIBRARIES} " )
#    message( STATUS "Want SimpleGui_INCLUDE_DIRS: ${${LibName}_INCLUDE_DIRS}")

    # Code to cleanup apple frameworks to "export" properly.
    # This is needed because spaces in the list of LIBRARIES in the .cmake we are
    # creating here MUST BE ESCAPPED.  Otherwise, when an external project imports
    # the list, CMAKE turns the space into a semicolon, resulting in mistakes
    # like "-lframework;-lIOKit" instead of "-framework IOKit".
    set( _${LibName}_LIBRARIES "" )
    foreach( lib ${${LibName}_LIBRARIES})
        string( REPLACE "framework " "framework\\ " lib ${lib} )
    list( APPEND _${LibName}_LIBRARIES ${lib} )
    endforeach()
    set( ${LibName}_LIBRARIES ${_${LibName}_LIBRARIES} CACHE INTERNAL "" FORCE )

    # add our library onto the front
    get_target_property( _LIBRARY ${LibName} LOCATION )
    if( NOT ${_LIBRARY} STREQUAL "_LIBRARY-NOTFOUND" )
        list( INSERT ${LibName}_LIBRARIES 0 ${_LIBRARY} )
    endif()

    file( WRITE ${LibName}Config.cmake.in "# automatically generated ${LibName}.cmake.in file\n" )
    if( NOT "${${LibName}_LIBRARIES}" STREQUAL "" )
        file( APPEND ${LibName}Config.cmake.in 
                "SET( ${LibName}_LIBRARIES  \@${LibName}_LIBRARIES\@ CACHE INTERNAL "
                "\"${LibName} libraries\" FORCE )\n"
                "mark_as_advanced( ${LibName}_LIBRARIES )\n"
            )
    endif()
    if( NOT "${${LibName}_INCLUDE_DIRS}" STREQUAL "" )
        file( APPEND ${LibName}Config.cmake.in 
            "SET( ${LibName}_INCLUDE_DIRS  \@${LibName}_INCLUDE_DIRS\@ CACHE INTERNAL "
            "\"${LibName} include directories\" FORCE )\n"
            "mark_as_advanced( ${LibName}_INCLUDE_DIRS )\n"
        )
    endif()
    if( NOT "${${LibName}_LIBRARY_DIRS}" STREQUAL "" )
        file( APPEND ${LibName}Config.cmake.in 
           "SET( ${LibName}_LIBRARY_DIRS  \@${LibName}_LIBRARY_DIRS\@ CACHE INTERNAL "
           "\"${LibName} library directories\" FORCE )\n"
           "mark_as_advanced( ${LibName}_LIBRARY_DIRS )\n"
           )
    endif()


    # Create the LibNameConfig.cmake file for the build tree.
    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/${LibName}Config.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/${LibName}Config.cmake @ONLY IMMEDIATE )

    # Add module to CMake package registery.
    export( PACKAGE ${LibName} )

endmacro()


