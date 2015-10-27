# setup some things initially
cmake_minimum_required (VERSION 2.6)
project (Moby)
include (CheckIncludeFiles)
include (CheckLibraryExists)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)

# setup include directories
include_directories ("include")

# setup library sources
set (SOURCES AABB.cpp ArticulatedBody.cpp Base.cpp BoundingSphere.cpp BoxPrimitive.cpp BV.cpp CCD.cpp CollisionDetection.cpp CollisionGeometry.cpp CompGeom.cpp ConePrimitive.cpp ConstraintSimulator.cpp ConstraintStabilization.cpp ContactParameters.cpp CP.cpp CylinderPrimitive.cpp DampingForce.cpp Dissipation.cpp DynamicBody.cpp FixedJoint.cpp Gears.cpp GJK.cpp GravityForce.cpp HeightmapPrimitive.cpp ImpactConstraintHandler.cpp ImpactConstraintHandlerNQP.cpp ImpactConstraintHandlerLCP.cpp ImpactConstraintHandlerQP.cpp IndexedTetraArray.cpp IndexedTriArray.cpp Joint.cpp LCP.cpp Log.cpp OBB.cpp OSGGroupWrapper.cpp PenaltyConstraintHandler.cpp PlanePrimitive.cpp PolyhedralPrimitive.cpp Polyhedron.cpp Primitive.cpp PrismaticJoint.cpp RCArticulatedBody.cpp RevoluteJoint.cpp RigidBody.cpp SDFReader.cpp Simulator.cpp SparseJacobian.cpp SpherePrimitive.cpp SphericalJoint.cpp SSL.cpp SSR.cpp StokesDragForce.cpp TessellatedPolyhedron.cpp Tetrahedron.cpp ThickTriangle.cpp TimeSteppingSimulator.cpp TorusPrimitive.cpp Triangle.cpp TriangleMeshPrimitive.cpp UnilateralConstraint.cpp UniversalJoint.cpp URDFReader.cpp Visualizable.cpp XMLReader.cpp XMLTree.cpp XMLWriter.cpp)
#set (SOURCES MCArticulatedBody.cpp)

# build options
option (BUILD_SHARED_LIBS "Build Moby as a shared library?" ON)
option (BUILD_TOOLS "Build the programs in the tools subdirectory?" ON)
option (USE_OSG "Build against OpenSceneGraph library?" ON)
option (USE_AP "Use the Anitescu-Potra contact model?" OFF)
option (USE_QLCPD "Build against QLCPD library?" OFF)
option (VISUALIZE_INERTIA "Visualize moments of inertia?" OFF)
option (PROFILE "Build for profiling?" OFF)

# look for QLCPD
find_library(QLCPD_FOUND qlcpd-dense /usr/local/lib /usr/lib)
if (USE_QLCPD AND QLCPD_FOUND)
  set (SOURCES ${SOURCES} QLCPD.cpp)
  set (EXTRA_LIBS qlcpd-dense)
  set_source_files_properties(src/ImpactConstraintHandler.cpp PROPERTIES COMPILE_FLAGS -DUSE_QLCPD)
  set_source_files_properties(src/ImpactConstraintHandlerQP.cpp PROPERTIES COMPILE_FLAGS -DUSE_QLCPD)
endif (USE_QLCPD AND QLCPD_FOUND)

# modify C++ flags
add_definitions (-DSAFESTATIC=static)
if (OMP)
  find_package (OpenMP REQUIRED)
  include_directories (${OPENMP_INCLUDE_DIRS})
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
endif (OMP)
if (PROFILE)
  set_source_files_properties(example/driver.cpp PROPERTIES COMPILE_FLAGS -DGOOGLE_PROFILER) 
  set (EXTRA_LIBS ${EXTRA_LIBS} profiler)
endif (PROFILE)

# fix the C++ linking error on 64-bit Linux
set (CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -ldl")

# find libraries
find_package (Ravelin REQUIRED)
find_package (LibXml2 REQUIRED)
find_package (Boost REQUIRED)
find_package (IPOPT)
get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
find_package (QHULL REQUIRED)
find_package (osg)
find_package (osgViewer)
find_package (osgDB)
find_package (osgGA)
find_package (osgFX)
find_package (osgUtil)
find_package (OpenThreads)

# see whether the Anitescu-Potra model should be activated
if (USE_AP)
  set_source_files_properties(src/ImpactConstraintHandler.cpp PROPERTIES COMPILE_FLAGS -DUSE_AP_MODEL) 
endif (USE_AP)
# see whether IPOPT was detected
if (USE_IPOPT AND IPOPT_FOUND)
  add_definitions (-DHAVE_CSTDDEF)
  add_definitions (-DHAVE_IPOPT)
  include_directories (IPOPT_INCLUDE_DIRS)
  link_directories (IPOPT_LIBRARY_DIRS)
  if (APPLE)
    find_package (BLAS REQUIRED)
  else (APPLE)
    find_package (CBLAS REQUIRED)
  endif (APPLE)
#  set (EXTRA_LIBS "${EXTRA_LIBS} ${BLAS_LIBRARIES}")
  find_library(IPOPT ipopt /usr/local/lib /usr/lib)
  find_library(IPOPT_MUMPS coinmumps /usr/local/lib /usr/lib)
  find_library(IPOPT_METIS coinmetis /usr/local/lib /usr/lib)
#  find_library(GFORTRAN gfortran /usr/local/lib /usr/lib)
#  set (EXTRA_LIBS ${EXTRA_LIBS} ${IPOPT_MUMPS} ${IPOPT_METIS} ${GFORTRAN})
  set (EXTRA_LIBS ${EXTRA_LIBS} ${IPOPT} ${IPOPT_MUMPS} ${IPOPT_METIS})
  set (SOURCES ${SOURCES} LCP_IPOPT.cpp NQP_IPOPT.cpp)
endif (USE_IPOPT AND IPOPT_FOUND)

# setup OSG includes and additional sources
if (USE_OSG AND OSG_FOUND)
  if (OSGDB_FOUND)
    if (OSGVIEWER_FOUND)
      if (OPENTHREADS_FOUND)
        include_directories (${OSG_INCLUDE_DIR})
        include_directories (${OSGDB_INCLUDE_DIR})
        include_directories (${OSGGA_INCLUDE_DIR})
        include_directories (${OSGVIEWER_INCLUDE_DIR})
        include_directories (${OPENTHREADS_INCLUDE_DIR})
        add_definitions (-DUSE_OSG)
        set (SOURCES ${SOURCES} OsgTorus.cpp)
      else (OPENTHREADS_FOUND)
        unset (OSG_FOUND)
      endif (OPENTHREADS_FOUND)
    else (OSGVIEWER_FOUND)
      unset (OSG_FOUND)
    endif (OSGVIEWER_FOUND)
  else (OSGDB_FOUND)
    unset (OSG_FOUND)
  endif (OSGDB_FOUND)
endif (USE_OSG AND OSG_FOUND)
if (NOT OSG_FOUND)
  set (USE_OSG FALSE)
endif (NOT OSG_FOUND)
if (NOT USE_OSG)
  set (VISUALIZE_INERTIA FALSE)
endif (NOT USE_OSG)
if (VISUALIZE_INERTIA)
  add_definitions (-DVISUALIZE_INERTIA)
endif (VISUALIZE_INERTIA)

# setup include directories, compiler flags, and libraries for optional pkgs
if (LIBXML2_FOUND)
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${LIBXML2_DEFINITIONS})
  include_directories (${LIBXML2_INCLUDE_DIR})
endif (LIBXML2_FOUND)

# setup include directories, compiler flags, and libraries for required pkgs
include_directories (${BOOST_INCLUDE_DIR})
include_directories (${QHULL_INCLUDE_DIR})
include_directories (${RAVELIN_INCLUDE_DIR})

# configure checks: optional libraries
CHECK_LIBRARY_EXISTS(odepack dlsode_ "" HAVE_ODEPACK)

# setup necessary library list
set (REQLIBS qhull)

# build against ODEPACK?
if (HAVE_ODEPACK)
  add_definitions (-DUSE_ODEPACK)
  set (EXTRA_LIBS ${EXTRA_LIBS} odepack)
endif (HAVE_ODEPACK)

# prepend "src/" to each source file
foreach (i ${SOURCES})
  set (LIBSOURCES ${LIBSOURCES} "${CMAKE_SOURCE_DIR}/src/${i}")
endforeach (i)

# create the library
add_library(Moby "" "" ${LIBSOURCES})
target_link_libraries (Moby ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${QHULL_LIBRARIES} ${RAVELIN_LIBRARIES} ${EXTRA_LIBS})

# link optional libraries
if (OMP)
  target_link_libraries (Moby ${OPENMP_LIBRARIES})
endif (OMP)
if (LIBXML2_FOUND)
  target_link_libraries (Moby ${LIBXML2_LIBRARIES})
endif (LIBXML2_FOUND)
if (USE_OSG AND OSG_FOUND)
  target_link_libraries (Moby ${OSG_LIBRARIES})
  target_link_libraries (Moby ${OSGVIEWER_LIBRARIES})
  target_link_libraries (Moby ${OSGDB_LIBRARIES})
  target_link_libraries (Moby ${OSGGA_LIBRARIES})
  target_link_libraries (Moby ${OPENTHREADS_LIBRARY})
  target_link_libraries (Moby ${OSGUTIL_LIBRARIES})
  target_link_libraries (Moby ${OSGMANIPULATOR_LIBRARIES})
endif (USE_OSG AND OSG_FOUND)

add_library(MobyDriver "" "" example/driver.cpp)
target_link_libraries (MobyDriver Moby)

# build the tools?
if (BUILD_TOOLS)
  add_executable(moby-driver example/main.cpp)
  if (USE_OSG AND OSG_FOUND)
    add_executable(moby-view example/view.cpp)
  endif (USE_OSG AND OSG_FOUND)

  # tools
  add_executable(moby-render example/render.cpp)
  add_executable(moby-regress example/regress.cpp)
  add_executable(moby-compare-trajs example/compare-trajs.cpp)
#  add_executable(moby-conv-decomp example/conv-decomp.cpp)
  add_executable(moby-convexify example/convexify.cpp)
  add_executable(moby-adjust-center example/adjust-center.cpp)
  add_executable(moby-center example/center.cpp)
  target_link_libraries(moby-driver MobyDriver Moby)
  if (USE_OSG AND OSG_FOUND)
    target_link_libraries(moby-render ${OSG_LIBRARIES})
    target_link_libraries(moby-render ${OSGVIEWER_LIBRARIES})
    target_link_libraries(moby-render ${OSGDB_LIBRARIES})
    target_link_libraries(moby-render ${OSGGA_LIBRARIES})
    target_link_libraries(moby-render ${OPENTHREADS_LIBRARIES})
    target_link_libraries(moby-view ${OSG_LIBRARIES})
    target_link_libraries(moby-view ${OSGVIEWER_LIBRARIES})
    target_link_libraries(moby-view ${OSGDB_LIBRARIES})
    target_link_libraries(moby-view ${OSGGA_LIBRARIES})
    target_link_libraries(moby-view ${OPENTHREADS_LIBRARIES})
  endif (USE_OSG AND OSG_FOUND)

  # tools
  target_link_libraries(moby-render Moby)
  target_link_libraries(moby-regress Moby)
  target_link_libraries(moby-compare-trajs Moby)
#  target_link_libraries(moby-conv-decomp Moby)
  target_link_libraries(moby-convexify Moby)
#  target_link_libraries(moby-output-symbolic Moby)
  target_link_libraries(moby-adjust-center Moby)
  target_link_libraries(moby-center Moby)
endif (BUILD_TOOLS)

# create environment variables file
file (WRITE ${CMAKE_BINARY_DIR}/setup.sh "#!/bin/bash\nexport MOBY_PLUGIN_PATH=${CMAKE_BINARY_DIR}\n\n")

# plugins
# sphere plugin
add_library(sphere-plugin MODULE example/contact_simple/sphere-plugin.cpp)
target_link_libraries(sphere-plugin Moby)

# box plugin
add_library(box-plugin MODULE example/contact_simple/box-plugin.cpp)
target_link_libraries(box-plugin Moby)

# rimless wheel plugins
add_library(rimless-wheel-coldet-plugin MODULE example/rimless-wheel/coldet-plugin.cpp)
add_library(rimless-wheel-init MODULE example/rimless-wheel/init.cpp)
target_link_libraries(rimless-wheel-init Moby)
target_link_libraries(rimless-wheel-coldet-plugin Moby)

# ball in parabloid plugins
add_library(ball-parabloid-coldet-plugin MODULE example/ball-in-parabloid/coldet-plugin.cpp)
add_library(ball-parabloid-init MODULE example/ball-in-parabloid/init.cpp)
target_link_libraries(ball-parabloid-init Moby)
target_link_libraries(ball-parabloid-coldet-plugin Moby)

# passive dynamic walker plugins
add_library(passive-walker-init MODULE example/passive-walker/init.cpp)
add_library(passive-walker-coldet-plugin MODULE example/passive-walker/coldet-plugin.cpp)
target_link_libraries(passive-walker-init Moby)
target_link_libraries(passive-walker-coldet-plugin Moby)

# contact constrained pendulum plugins
add_library(contact-constrained-pendulum-coldet-plugin MODULE example/contact_simple/contact-constrained-pendulum-coldet-plugin.cpp)
add_library(contact-constrained-pendulum-init MODULE example/contact_simple/contact-constrained-pendulum-init.cpp)
target_link_libraries(contact-constrained-pendulum-coldet-plugin Moby)
target_link_libraries(contact-constrained-pendulum-init Moby)

# setup source file
install (CODE "file(WRITE setup.sh \"#!/bin/bash\nexport MOBY_PLUGIN_PATH=${CMAKE_INSTALL_PREFIX}/lib/Moby/plugins\n\n\")")
install (FILES ${CMAKE_BINARY_DIR}/setup.sh DESTINATION share/Moby)

# setup install locations for plugins
install (TARGETS passive-walker-init DESTINATION lib/Moby/plugins)
install (TARGETS passive-walker-coldet-plugin DESTINATION lib/Moby/plugins)
install (TARGETS contact-constrained-pendulum-coldet-plugin DESTINATION lib/Moby/plugins)
install (TARGETS contact-constrained-pendulum-init DESTINATION lib/Moby/plugins)
install (TARGETS passive-walker-init DESTINATION lib/Moby/plugins)
install (TARGETS passive-walker-coldet-plugin DESTINATION lib/Moby/plugins)
install (TARGETS ball-parabloid-init DESTINATION lib/Moby/plugins)
install (TARGETS ball-parabloid-coldet-plugin DESTINATION lib/Moby/plugins)
install (TARGETS rimless-wheel-init DESTINATION lib/Moby/plugins)
install (TARGETS rimless-wheel-coldet-plugin DESTINATION lib/Moby/plugins)
install (TARGETS box-plugin DESTINATION lib/Moby/plugins)
install (TARGETS sphere-plugin DESTINATION lib/Moby/plugins)

# setup install locations for library
install (TARGETS Moby DESTINATION lib)
install (TARGETS MobyDriver DESTINATION lib)

# setup install locations for binaries
install (TARGETS moby-driver DESTINATION bin)
if (USE_OSG AND OSG_FOUND)
  install (TARGETS moby-view DESTINATION bin)
  install (TARGETS moby-render DESTINATION bin)
endif (USE_OSG AND OSG_FOUND)
#install (TARGETS moby-conv-decomp DESTINATION bin)
install (TARGETS moby-convexify DESTINATION bin)
install (TARGETS moby-adjust-center DESTINATION bin)
install (TARGETS moby-center DESTINATION bin)

# setup install locations for headers
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/Moby DESTINATION include)

# setup install locations for examples
install (DIRECTORY ${CMAKE_SOURCE_DIR}/example/contact_simple DESTINATION share/Moby)
install (DIRECTORY ${CMAKE_SOURCE_DIR}/example/chain-contact DESTINATION share/Moby)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

