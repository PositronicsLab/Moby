# setup some things initially
cmake_minimum_required (VERSION 2.6)
project (Moby)
include (CheckIncludeFiles)
include (CheckLibraryExists)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)

# setup include directories
include_directories ("include")
include_directories ("/Users/drum/repo/Ravelin/include")

# setup link directories
link_directories ("/Users/drum/repo/Ravelin/debug")

# setup library sources
set (SOURCES DeformableBody.cpp GravityForce.cpp CRBAlgorithm.cpp RNEAlgorithm.cpp RCArticulatedBody.cpp RigidBody.cpp LCP.cpp OSGGroupWrapper.cpp Integrator.cpp SpherePrimitive.cpp ConePrimitive.cpp CylinderPrimitive.cpp BoxPrimitive.cpp GaussianMixture.cpp CSG.cpp FSABAlgorithm.cpp Simulator.cpp EventDrivenSimulator.cpp Spatial.cpp ArticulatedBody.cpp Primitive.cpp Visualizable.cpp AABB.cpp SSL.cpp SSR.cpp BoundingSphere.cpp CompGeom.cpp Polyhedron.cpp Integrator.cpp ODEPACKIntegrator.cpp EulerIntegrator.cpp VariableStepIntegrator.cpp Rosenbrock4Integrator.cpp RungeKuttaIntegrator.cpp RungeKuttaFehlbergIntegrator.cpp RungeKuttaImplicitIntegrator.cpp ContactParameters.cpp Log.cpp XMLWriter.cpp ThickTriangle.cpp Triangle.cpp Tetrahedron.cpp Base.cpp DynamicBody.cpp SingleBody.cpp XMLTree.cpp FixedJoint.cpp XMLReader.cpp DampingForce.cpp RevoluteJoint.cpp CollisionDetection.cpp CollisionGeometry.cpp MeshDCD.cpp PrismaticJoint.cpp SphericalJoint.cpp UniversalJoint.cpp BV.cpp C2ACCD.cpp OBB.cpp IndexedTetraArray.cpp TriangleMeshPrimitive.cpp IndexedTriArray.cpp Joint.cpp ImpactEventHandlerQP.cpp ImpactEventHandler.cpp Event.cpp StokesDragForce.cpp VariableEulerIntegrator.cpp GeneralizedCCD.cpp)
#set (SOURCES URDFReader.cpp MCArticulatedBody.cpp PSDeformableBody.cpp DeformableCCD.cpp)

# build options 
option (BUILD_SHARED_LIBS "Build Moby as a shared library?" ON)
option (BUILD_TOOLS "Build the programs in the tools subdirectory?" ON)
option (USE_OSG "Build against OpenSceneGraph library?" ON)
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
get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
find_package (QHULL REQUIRED)
find_package (osg)
find_package (osgViewer)
find_package (osgDB)
find_package (osgGA)
find_package (OpenThreads)

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
include_directories (${BOOST_INCLUDE_DIRS})
include_directories (${QHULL_INCLUDE_DIRS})

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
  add_library(grasp-controller "MODULE" example/banditII/grasp-controller.cpp)
  add_library(bandit-controller "MODULE" example/banditII/bandit-controller.cpp)
  add_library(pinch-controller "MODULE" example/banditII/pinch-controller.cpp)
  add_library(push-controller "MODULE" example/banditII/push-controller.cpp)
  add_library(noball-controller "MODULE" example/banditII/noball-controller.cpp)
  add_library(pinch-release-controller "MODULE" example/banditII/pinch-release-controller.cpp)
  add_library(block-collider-controller "MODULE" example/banditII/block-collider-controller.cpp)

  # tools
  add_executable(moby-regress example/regress.cpp)
  add_executable(moby-compare-trajs example/compare-trajs.cpp)
  add_executable(moby-conv-decomp example/conv-decomp.cpp)
  add_executable(moby-convexify example/convexify.cpp)
  add_executable(moby-output-symbolic example/output-symbolic.cpp)
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
  target_link_libraries(grasp-controller Moby)
  target_link_libraries(bandit-controller Moby)
  target_link_libraries(pinch-controller Moby)
  target_link_libraries(push-controller Moby)
  target_link_libraries(noball-controller Moby)
  target_link_libraries(pinch-release-controller Moby)
  target_link_libraries(block-collider-controller Moby)

  # tools
  target_link_libraries(moby-regress Moby)
  target_link_libraries(moby-compare-trajs Moby)
  target_link_libraries(moby-conv-decomp Moby)
  target_link_libraries(moby-convexify Moby)
  target_link_libraries(moby-output-symbolic Moby)
  target_link_libraries(moby-adjust-center Moby)
  target_link_libraries(moby-center Moby)
endif (BUILD_TOOLS)

# setup install locations
install (TARGETS Moby DESTINATION lib)
install (TARGETS moby-driver DESTINATION bin)
if (USE_OSG AND OSG_FOUND)
  install (TARGETS moby-view DESTINATION bin)
endif (USE_OSG AND OSG_FOUND)
install (TARGETS moby-conv-decomp DESTINATION bin)
install (TARGETS moby-convexify DESTINATION bin)
install (TARGETS moby-adjust-center DESTINATION bin)
install (TARGETS moby-center DESTINATION bin)
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/Moby DESTINATION include)

