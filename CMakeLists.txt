# setp some things initially
cmake_minimum_required (VERSION 2.6)
project (Moby)
include (CheckIncludeFiles)
include (CheckLibraryExists)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)

# setup include directories
include_directories ("include")

# setup library sources
set (SOURCES AABB.cpp AccelerationEventHandler.cpp ADF.cpp ArticulatedBody.cpp Base.cpp BoundingSphere.cpp BoxPrimitive.cpp BulirschStoerIntegrator.cpp BV.cpp CCD.cpp CollisionGeometry.cpp CompGeom.cpp ConePrimitive.cpp ContactParameters.cpp CRBAlgorithm.cpp CylinderPrimitive.cpp DampingForce.cpp DynamicBody.cpp EulerIntegrator.cpp Event.cpp EventDrivenSimulator.cpp FixedJoint.cpp FSABAlgorithm.cpp GJK.cpp GravityForce.cpp HeightmapPrimitive.cpp ImpactEventHandler.cpp ImpactEventHandlerNQP.cpp ImpactEventHandlerQP.cpp IndexedTetraArray.cpp IndexedTriArray.cpp Integrator.cpp Joint.cpp LCP.cpp Log.cpp OBB.cpp ODEPACKIntegrator.cpp OSGGroupWrapper.cpp Polyhedron.cpp Primitive.cpp PrismaticJoint.cpp RCArticulatedBody.cpp RevoluteJoint.cpp RigidBody.cpp RNEAlgorithm.cpp Rosenbrock4Integrator.cpp RungeKuttaFehlbergIntegrator.cpp RungeKuttaIntegrator.cpp RungeKuttaImplicitIntegrator.cpp Simulator.cpp SingleBody.cpp Spatial.cpp SpherePrimitive.cpp SphericalJoint.cpp SSL.cpp SSR.cpp StokesDragForce.cpp Tetrahedron.cpp ThickTriangle.cpp Triangle.cpp TriangleMeshPrimitive.cpp UniversalJoint.cpp URDFReader.cpp VariableEulerIntegrator.cpp VariableStepIntegrator.cpp Visualizable.cpp XMLReader.cpp XMLTree.cpp XMLWriter.cpp)
#set (SOURCES MCArticulatedBody.cpp)

# build options 
option (BUILD_SHARED_LIBS "Build Moby as a shared library?" ON)
option (BUILD_TOOLS "Build the programs in the tools subdirectory?" ON)
option (USE_OSG "Build against OpenSceneGraph library?" ON)
option (USE_IPOPT "Build against IPOPT library?" ON)
option (PROFILE "Build for profiling?" OFF)
option (OMP "Build with OpenMP support?" OFF)
option (THREADSAFE "Build Moby to be threadsafe? (slower)" OFF)

# check options are valid
if (THREADSAFE)
  if (OMP)
    unset (OMP)
  endif (OMP)
endif (THREADSAFE)

# modify C++ flags
if (THREADSAFE)
  add_definitions (-DSAFESTATIC=)
  add_definitions (-DTHREADSAFE)
else (THREADSAFE)
  add_definitions (-DSAFESTATIC=static)
endif (THREADSAFE)
if (OMP)
  find_package (OpenMP REQUIRED)
  include_directories (${OPENMP_INCLUDE_DIRS})
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
endif (OMP)
if (PROFILE)
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-pg -g")
endif (PROFILE)

# find libraries
find_package (LibXml2 REQUIRED)
find_package (Boost REQUIRED)
find_package (IPOPT)
get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
find_package (QHULL REQUIRED)
find_package (osg)
find_package (osgViewer)
find_package (osgDB)
find_package (osgGA)
find_package (OpenThreads)

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

# setup include directories, compiler flags, and libraries for optional pkgs
if (LIBXML2_FOUND)
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${LIBXML2_DEFINITIONS})
  include_directories (${LIBXML2_INCLUDE_DIR})
endif (LIBXML2_FOUND)

# setup include directories, compiler flags, and libraries for required pkgs
include_directories (${BOOST_INCLUDE_DIR})
include_directories (${QHULL_INCLUDE_DIR})

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
target_link_libraries (Moby ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${QHULL_LIBRARIES} Ravelin ${EXTRA_LIBS})

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
endif (USE_OSG AND OSG_FOUND)  

# build the tools?
if (BUILD_TOOLS)
  add_executable(moby-driver example/driver.cpp)
  if (USE_OSG AND OSG_FOUND)
    add_executable(moby-view example/view.cpp)
  endif (USE_OSG AND OSG_FOUND)

  # plugins
  #add_library(grasp-controller "MODULE" example/banditII/grasp-controller.cpp)
  #add_library(bandit-controller "MODULE" example/banditII/bandit-controller.cpp)
  #add_library(pinch-controller "MODULE" example/banditII/pinch-controller.cpp)
  #add_library(push-controller "MODULE" example/banditII/push-controller.cpp)
  #add_library(noball-controller "MODULE" example/banditII/noball-controller.cpp)
  #add_library(pinch-release-controller "MODULE" example/banditII/pinch-release-controller.cpp)
  #add_library(block-collider-controller "MODULE" example/banditII/block-collider-controller.cpp)

  # tools
  add_executable(moby-regress example/regress.cpp)
  add_executable(moby-compare-trajs example/compare-trajs.cpp)
#  add_executable(moby-conv-decomp example/conv-decomp.cpp)
  add_executable(moby-convexify example/convexify.cpp)
#  add_executable(moby-output-symbolic example/output-symbolic.cpp)
  add_executable(moby-adjust-center example/adjust-center.cpp)
  add_executable(moby-center example/center.cpp)
  target_link_libraries(moby-driver Moby)
  if (USE_OSG AND OSG_FOUND)
    target_link_libraries(moby-view ${OSG_LIBRARIES})
    target_link_libraries(moby-view ${OSGVIEWER_LIBRARIES})
    target_link_libraries(moby-view ${OSGDB_LIBRARIES})
    target_link_libraries(moby-view ${OSGGA_LIBRARIES})
    target_link_libraries(moby-view ${OPENTHREADS_LIBRARIES})
  endif (USE_OSG AND OSG_FOUND)

  # uncomment to build example plugins
#  target_link_libraries(controller Moby)
  #target_link_libraries(grasp-controller Moby)
  #target_link_libraries(bandit-controller Moby)
  #target_link_libraries(pinch-controller Moby)
  #target_link_libraries(push-controller Moby)
  #target_link_libraries(noball-controller Moby)
  #target_link_libraries(pinch-release-controller Moby)
  #target_link_libraries(block-collider-controller Moby)

  # tools
  target_link_libraries(moby-regress Moby)
  target_link_libraries(moby-compare-trajs Moby)
#  target_link_libraries(moby-conv-decomp Moby)
  target_link_libraries(moby-convexify Moby)
#  target_link_libraries(moby-output-symbolic Moby)
  target_link_libraries(moby-adjust-center Moby)
  target_link_libraries(moby-center Moby)
endif (BUILD_TOOLS)

# setup install locations
install (TARGETS Moby DESTINATION lib)
install (TARGETS moby-driver DESTINATION bin)
if (USE_OSG AND OSG_FOUND)
  install (TARGETS moby-view DESTINATION bin)
endif (USE_OSG AND OSG_FOUND)
#install (TARGETS moby-conv-decomp DESTINATION bin)
install (TARGETS moby-convexify DESTINATION bin)
install (TARGETS moby-adjust-center DESTINATION bin)
install (TARGETS moby-center DESTINATION bin)
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/Moby DESTINATION include)

