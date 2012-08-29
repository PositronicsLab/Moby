import os
import string

print 'WARNING: Building Moby through SCons is now deprecated!'
print '         The CMake build process is the new way to build Moby'

##############################################################################
# WARNING: if include / library paths are not being found, do not alter this
# file; edit BUILD-OPTIONS.py instead. 
##############################################################################

# determines whether pkg-config exists
def CheckPKGConfig(context):
  context.Message('Checking for pkg-config... ')
  ret = context.TryAction('pkg-config --version')[0]
  context.Result(ret)
  return ret

# determines whether a given package exists in pkg-config
def CheckPKG(context, name):
  context.Message('Checking for (pkg-config) %s...' % name)
  ret = context.TryAction('pkg-config --exists \'%s\'' % name)[0]
  context.Result(ret)
  return ret

# determines whether a given executable is in the path
def inpath(executable):
  fout, fin = os.popen2("which " + executable)
  fout.close()
  s = fin.readline();
  if (s == ""):
    return 0
  slist = s.split();
  if (slist[0] == "no"):
    return 0
  else:
    return 1

# setup an options variable
vars = Variables("BUILD-OPTIONS.py")

# set default paths
DEFAULT_INCLUDE_PATHS = ['/usr/include/libxml2']
DEFAULT_LIB_PATHS = ['/usr/lib']
DEFAULT_INSTALL_PATH = '/usr/local'

# get options 
vars.Add('CXX', 'C++ compiler', 'g++')
vars.Add('CXXFLAGS', 'Additional C++ compiler options', '')
vars.Add(EnumVariable('REAL_TYPE', 'Type of floating number', 'double', allowed_values=('float', 'double', 'arbitrary')))
vars.Add('PRECISION', 'Arbitrary precision bits', '128')
vars.Add(BoolVariable('THREAD_SAFE', 'Set to true to build Moby so that it is thread safe (slower)', 0))
vars.Add(BoolVariable('SHARED_LIBRARY', 'Set to true to build Moby as a shared library', 0))
vars.Add(BoolVariable('USE_OSG', 'Set to true to build Moby to use OpenSceneGraph visualization', 1))
vars.Add(BoolVariable('USE_OPENMP', 'Set to true to build Moby to use OpenMP multithreading', 1))
vars.Add(BoolVariable('BUILD_EXAMPLES', 'Set to true to build example controllers and utilities in the examples directory', 1))
vars.Add(BoolVariable('DEBUG', 'Set to false to build optimized', 1))
vars.Add(BoolVariable('PROFILE', 'Set to true to build for profiling', 0))
vars.Add('INSTALL_PATH', 'Root path to which to install the Moby library and header files', DEFAULT_INSTALL_PATH)
vars.Add(BoolVariable('USE_PATH', 'Set to true to build to use the PATH solver', 0))
vars.Add('INCLUDE_PATHS', 'Additional, colon separated paths to search for include files', '')
vars.Add('LIB_PATHS', 'Additional, colon separated paths to search for libraries', '')
vars.Add('FRAMEWORK_PATHS', 'Additional, colon separated paths to search for frameworks (OS X only)', '')

# setup a construction environment
env = Environment(variables=vars)
Help(vars.GenerateHelpText(env))

# see whether to build with OSG support
USE_OSG = env['USE_OSG']

# setting some options causes others to be altered
if not USE_OSG:
  env['BUILD_EXAMPLES']=False

# do a pre-configure
preconf = Configure(env, custom_tests = { 'CheckPKGConfig' : CheckPKGConfig,
                                          'CheckPKG' : CheckPKG })

# initially indicate unable to use pkg-config for the following configurations
USE_LIBXML2_PKGCONF = False

# see whether we can use pkg-config to determine paths
if preconf.CheckPKGConfig():
  if preconf.CheckPKG('libxml-2.0'):
    USE_LIBXML2_PKGCONF = True

# attempt to determine paths via pkg-config
env = preconf.Finish()

# setup C++ compiler
__CXX = env['CXX']
__CXXFLAGS = env['CXXFLAGS']

# set link flags to empty initially
__LINKFLAGS = ""

# turn fPIC on -- necessary for dynamic libs that want to link w/Moby
__CXXFLAGS = __CXXFLAGS + ' -fPIC'

# turn all warnings on
__CXXFLAGS = __CXXFLAGS + ' -Wall'

# see whether we are installing Moby
INSTALL_PATH = env['INSTALL_PATH']

# determine whether we are using the PATH solver
USE_PATH = bool(env['USE_PATH'])
if USE_PATH:
  __CXXFLAGS = __CXXFLAGS + " -DUSE_PATH "

# set appropriate options for debug / release
if bool(env['DEBUG']):
  __CXXFLAGS = __CXXFLAGS + " -g "
else:
  __CXXFLAGS = __CXXFLAGS + " -g -O2 "

# set appropriate option for OpenMP
if bool(env['USE_OPENMP']):
  __CXXFLAGS = __CXXFLAGS + " -fopenmp "
  __LINKFLAGS = __LINKFLAGS + " -fopenmp "

# see whether to build with OSG support
if bool(USE_OSG):
  __CXXFLAGS = __CXXFLAGS + ' -DUSE_OSG'  

# see whether we are building thread safe
if bool(env['THREAD_SAFE']):
  __CXXFLAGS = __CXXFLAGS + ' -DSAFESTATIC'
else:
  __CXXFLAGS = __CXXFLAGS + ' -DSAFESTATIC=static'

# see what type to use for Real
REAL_TYPE = env['REAL_TYPE']
if REAL_TYPE == "double":
  __CXXFLAGS = __CXXFLAGS + ' -DBUILD_DOUBLE'
  FLOAT_PRECISION = 64
elif REAL_TYPE == "float":
  __CXXFLAGS = __CXXFLAGS + ' -DBUILD_SINGLE'
  FLOAT_PRECISION = 32
else:
  __CXXFLAGS = __CXXFLAGS + ' -DBUILD_ARBITRARY_PRECISION'
  if env['PRECISION']:
    FLOAT_PRECISION = int(env['PRECISION'])
    __CXXFLAGS = __CXXFLAGS + ' -DPRECISION=' + str(FLOAT_PRECISION)

# see whether to build Kinematic and examples
BUILD_EXAMPLES = bool(env['BUILD_EXAMPLES'])
SHARED_LIBRARY = bool(env['SHARED_LIBRARY'])

# set appropriate options for profiling
if bool(env['PROFILE']):
  __CXXFLAGS = __CXXFLAGS + " -g -pg -DNDEBUG -O2 "
  __LINKFLAGS = __LINKFLAGS + "-pg "

# setup constant paths
moby_include_path = "#/include"
moby_lib_path = "#/lib"

# setup the include path
inc_path = [ moby_include_path ] 
for i in env['INCLUDE_PATHS'].split(':'):
  inc_path.append(i)

# setup the library path
lib_path = [ moby_lib_path ]
for i in env['LIB_PATHS'].split(':'):
  lib_path.append(i)

# setup the frameworks path
framework_path = env['FRAMEWORK_PATHS'].split(':')

# get the platform
platform = str(Platform())

# determine set of frameworks 
if platform == 'darwin':
  __FRAMEWORKS = ['vecLib']
  __FRAMEWORKS.append(['OpenGL'])

# now, setup library paths if we are building examples, kinematic, or shared
# library...  this will differ depending on whether we are building on linux 
# or OS X we will also re-init the environment for configuration
if BUILD_EXAMPLES or SHARED_LIBRARY:
  if platform == 'darwin' or platform == 'posix':
    benv = Environment(CXX=__CXX, CXXFLAGS=__CXXFLAGS, CPPPATH=inc_path, LIBPATH=lib_path, FRAMEWORKPATH=framework_path, FRAMEWORKS=__FRAMEWORKS)

  else:
    print 'Platform ' + platform + ' is unsupported!'

else:
  # setup environment without libraries
  benv = Environment(CXX=__CXX, CXXFLAGS=__CXXFLAGS, CPPPATH=inc_path)

# determine paths setup through pkg-config
if USE_LIBXML2_PKGCONF:
  benv.ParseConfig('pkg-config --cflags --libs libxml-2.0')

# begin configuration
conf = Configure(benv)

# verify that include files can be found
if conf.CheckCHeader('glpk.h'):
  USE_GLPK = True
else:
  USE_GLPK = False
if REAL_TYPE=="arbitrary":
  if not conf.CheckCHeader('gmp.h'):
    print 'Unable to find GMP -- build without arbitrary precision or set '
    print 'INCLUDE_PATHS to point to GMP'
    Exit(1)
  if not conf.CheckCHeader('mpfr.h'):
    print 'Unable to find MPFR -- build without arbitrary precision or set '
    print 'INCLUDE_PATHS to point to MPFR'
    Exit(1)
if not conf.CheckCXXHeader('map'):
  print 'Unable to find STL -- C++ compiler non-standard?'
  Exit(1)
if not conf.CheckCHeader('libxml/parser.h'):
  print 'Unable to find libxml2 includes -- set INCLUDE_PATHS correctly!'
  Exit(1)
if not conf.CheckCXXHeader('boost/shared_ptr.hpp'):
  print 'Unable to find BOOST includes -- set INCLUDE_PATHS correctly!'
  Exit(1)
if not conf.CheckCHeader(['stdio.h', 'qhull/qhull_a.h']):
  print 'Unable to find QHULL includes -- set INCLUDE_PATHS correctly!'
  Exit(1)
if USE_OSG and not conf.CheckCXXHeader('osg/MatrixTransform'):
  print 'Unable to find OpenSceneGraph includes -- set INCLUDE_PATHS correctly!'
  Exit(1)

# verify that libraries can be found -- only necessary if building examples,
# or a shared library
if BUILD_EXAMPLES or SHARED_LIBRARY:
  if USE_GLPK and not conf.CheckLib('glpk'):
    print 'Found GLPK header but not the library!'
    USE_GLPK = False
  if not conf.CheckLib('odepack'):
    USE_ODEPACK = False
  else:
    USE_ODEPACK = True
    __CXXFLAGS = __CXXFLAGS + " -DUSE_ODEPACK"
  if not conf.CheckLib('xml2'):
    print 'Unable to find libxml2 library -- set LIB_PATHS correctly!'
    Exit(1)
  if not conf.CheckLib('qhull'):
    print 'Unable to find qhull library -- set LIB_PATHS correctly!'
    Exit(1)

  if REAL_TYPE == "arbitrary": 
    if not conf.CheckLib('libgmp'):
      print 'Unable to find GMP library -- set LIB_PATHS correctly!'
      Exit(1)
    if not conf.CheckLib('libmpfr'):
      print 'Unable to find MPFR library -- set LIB_PATHS correctly!'
      Exit(1)
  if USE_PATH and not conf.CheckLib(['path46', 'c']):
    print 'Unable to find PATH library -- set LIB_PATHS correctly!'
    Exit(1)
  if conf.CheckLib('lapack'):
    LAPACK_LIB = 'lapack'
  elif conf.CheckLib('lapack-3'):
    LAPACK_LIB = 'lapack-3'
  elif conf.CheckLib('lapack-3.1'):
    LAPACK_LIB = 'lapack-3.1'
  else:
    print 'Unable to find lapack library -- set LIB_PATHS correctly!'
    Exit(1)

# finish configuration
benv = conf.Finish()

# see whether to add USE_GLPK to the compile flags
if USE_GLPK:
  __CXXFLAGS = __CXXFLAGS + " -DUSE_GLPK"

# now, setup libraries if we are building examples...  this 
# will differ depending on whether we are building on linux or OS X
# we will also re-init the environment for configuration
if BUILD_EXAMPLES or SHARED_LIBRARY:
  if platform == 'darwin':
    libs = [ 'xml2' ]
    if USE_GLPK:
      libs.append('glpk')
    if USE_ODEPACK:
      libs.append('odepack')
    libs.extend(['pthread', 'qhull', 'm', 'dl'])
    env = Environment(CXX=__CXX, CXXFLAGS=__CXXFLAGS, CPPPATH=inc_path, 
                      FRAMEWORKS=__FRAMEWORKS, FRAMEWORKPATH=framework_path, 
                      LIBS=libs, LIBPATH=lib_path, LINKFLAGS=__LINKFLAGS, 
                      TOOLS=['default', 'qt'])
  elif platform == 'posix':
    # setup the library paths
    libs = [ 'xml2' ]
    if USE_PATH:
      libs.append('path46')
    if REAL_TYPE == "arbitrary":
      libs.extend(['mpfr', 'gmp'])
    if USE_ODEPACK:
      libs.append('odepack')
    if USE_GLPK:
      libs.append('glpk')
    libs.extend(['dl', 'pthread', LAPACK_LIB, 'cblas', 'qhull', 'm'])
    env = Environment(CXX=__CXX, CPPPATH=inc_path, LIBS=libs, LIBPATH=lib_path,
                      CXXFLAGS=__CXXFLAGS, LINKFLAGS=__LINKFLAGS, 
                      TOOLS=['default', 'qt'])
                      
  else:
    print 'Platform ' + platform + ' is unsupported!'
# just building the library itself
else:
    env = Environment(CXX=__CXX, CPPPATH=inc_path, CXXFLAGS=__CXXFLAGS) 
                      
# again determine paths setup through pkg-config
if USE_LIBXML2_PKGCONF:
  env.ParseConfig('pkg-config --cflags --libs libxml-2.0')

# setup the sources
sources = ['src/Base.cpp', 'src/ArticulatedBody.cpp',  
      'src/CollisionGeometry.cpp', 'src/Event.cpp', 
      'src/BoxPrimitive.cpp', 'src/CylinderPrimitive.cpp',
      'src/FSABAlgorithm.cpp', 'src/ThickTriangle.cpp', 'src/ConePrimitive.cpp',
      'src/GravityForce.cpp', 'src/cblas.cpp', 
      'src/LinAlg.cpp', 'src/AABB.cpp', 'src/CSG.cpp', 
      'src/BV.cpp', 'src/EventDrivenSimulator.cpp', 
      'src/Vector2.cpp', 'src/Vector3.cpp', 'src/SVector6.cpp',
      'src/VectorN.cpp', 'src/Matrix2.cpp', 'src/Matrix3.cpp', 
      'src/Matrix4.cpp', 'src/MatrixN.cpp', 'src/MatrixNN.cpp',
      'src/SMatrix6N.cpp', 'src/SMatrix6.cpp', 'src/Quat.cpp',
      'src/AAngle.cpp', 'src/ContactParameters.cpp',
      'src/BoundingSphere.cpp', 'src/ImpactEventHandler.cpp', 
      'src/Primitive.cpp', 'src/Integrator.cpp', 'src/StokesDragForce.cpp',
      'src/FixedJoint.cpp', 'src/DeformableCCD.cpp',
      'src/RCArticulatedBody.cpp', 
      'src/Joint.cpp', 'src/DampingForce.cpp',
      'src/Optimization.cpp', 'src/DynamicBody.cpp',
      'src/RigidBody.cpp',
      'src/Simulator.cpp', 
      'src/Triangle.cpp', 'src/TriangleMeshPrimitive.cpp',
      'src/Log.cpp', 
      'src/CRBAlgorithm.cpp', 'src/DeformableBody.cpp', 'src/Tetrahedron.cpp',
      'src/IndexedTetraArray.cpp', 
      'src/SpherePrimitive.cpp', 'src/RevoluteJoint.cpp',
      'src/CompGeom.cpp', 'src/Polyhedron.cpp', 'src/SparseMatrixN.cpp', 
      'src/CollisionDetection.cpp', 'src/RNEAlgorithm.cpp', 
      'src/XMLTree.cpp', 'src/MCArticulatedBody.cpp', 'src/SparseVectorN.cpp',
      'src/XMLWriter.cpp', 'src/PrismaticJoint.cpp', 
      'src/XMLReader.cpp', 'src/MeshDCD.cpp', 'src/PSDeformableBody.cpp',
      'src/RCArticulatedBodyFwdDynAlgo.cpp', 
      'src/SphericalJoint.cpp', 'src/SpatialTransform.cpp', 
      'src/UniversalJoint.cpp', 'src/SpatialRBInertia.cpp',
      'src/SpatialABInertia.cpp',
      'src/OBB.cpp', 'src/IndexedTriArray.cpp', 'src/SSL.cpp',
      'src/Visualizable.cpp',
      'src/GeneralizedCCD.cpp', 'src/SSR.cpp', 'src/C2ACCD.cpp']

# add mpreal++ only if building with arbitrary precision
if REAL_TYPE == "arbitrary":
  sources.append(['src/mpreal.cpp', 'src/blas-ap.cpp', 'src/lapack-ap.cpp', 'src/f2c-ap.cpp'])

# add the path solver interface to the list, if desired
if USE_PATH:
  sources.append('src/PathLCPSolver.cpp')

# add necessary sources for visualization
if USE_OSG:
  sources.append('src/OSGGroupWrapper.cpp')

# sort the list so we know the build order
sources.sort();

# build the Moby library
if SHARED_LIBRARY:
	libfname = env.SharedLibrary('lib/libmoby', sources)
else:
	libfname = env.Library('lib/libmoby', sources)

# build the examples and utilities
if BUILD_EXAMPLES:
	env_copy = env.Clone()
	SConscript(['example/SConscript'], exports=['env_copy', 'USE_OSG'])

# setup list of binaries
if BUILD_EXAMPLES:
	binaries = [ 	'example/center', 'example/convexify',
			'example/driver', 
			'example/extract-contacts.py' ]
else:
	binaries = []

# setup list of headers
headers = [	'include/Moby/AAngle.h', 
		'include/Moby/Base.h',
		'include/Moby/BoundingSphere.h',
		'include/Moby/BoxPrimitive.h',
		'include/Moby/BV.h',
		'include/Moby/BV.inl',
		'include/Moby/C2ACCD.h',
		'include/Moby/C2ACCD.inl',
		'include/Moby/CollisionDetection.h',
		'include/Moby/CollisionDetection.inl',
		'include/Moby/CollisionGeometry.h',
		'include/Moby/CollisionGeometry.inl',
		'include/Moby/CollisionMethod.h',
		'include/Moby/CollisionMethod.inl',
		'include/Moby/CompGeom.h',
		'include/Moby/CompGeom.inl',
		'include/Moby/ConePrimitive.h',
		'include/Moby/Constants.h',
		'include/Moby/CRBAlgorithm.h',
		'include/Moby/CylinderPrimitive.h',
		'include/Moby/DummyBV.h',
		'include/Moby/DeformableBody.h',
		'include/Moby/DeformableBody.inl',
		'include/Moby/DegenerateTriangleException.h',
		'include/Moby/DynamicBody.h',
		'include/Moby/EulerIntegrator.h',
		'include/Moby/EulerIntegrator.inl',
		'include/Moby/Event.h',
		'include/Moby/Event.inl',
		'include/Moby/FixedJoint.h',
		'include/Moby/FSABAlgorithm.h',
		'include/Moby/GeneralizedCCD.h',
		'include/Moby/GeneralizedCCD.inl',
		'include/Moby/GravityForce.h',
		'include/Moby/IndexedTetraArray.h',
		'include/Moby/IndexedTetraArray.inl',
		'include/Moby/IndexedTriArray.h',
		'include/Moby/IndexedTriArray.inl',
		'include/Moby/IndexedTri.h',
		'include/Moby/Integrator.h',
		'include/Moby/Integrator.inl',
		'include/Moby/InvalidIndexException.h',
		'include/Moby/InvalidTransformException.h',
		'include/Moby/Joint.h',
		'include/Moby/LinAlg.h',
		'include/Moby/Log.h',
		'include/Moby/mpreal.h',
		'include/Moby/Matrix2.h',
		'include/Moby/Matrix3.h',
		'include/Moby/Matrix4.h',
		'include/Moby/MatrixN.h',
		'include/Moby/MatrixN.inl',
		'include/Moby/MatrixNN.h',
		'include/Moby/MatrixNN.inl',
		'include/Moby/MCArticulatedBody.h',
		'include/Moby/MeshDCD.h',
		'include/Moby/MeshDCD.inl',
		'include/Moby/MissizeException.h',
		'include/Moby/NonconvexityException.h',
		'include/Moby/NullPointerException.h',
		'include/Moby/NumericalException.h',
		'include/Moby/OBB.h',
		'include/Moby/OBB.inl',
		'include/Moby/Optimization.h',
		'include/Moby/OSGGroupWrapper.h',
		'include/Moby/PathLCPSolver.h',
		'include/Moby/Plane.h',
		'include/Moby/PlanePrimitive.h',
		'include/Moby/Polyhedron.h',
		'include/Moby/Polyhedron.inl',
		'include/Moby/Primitive.h',
		'include/Moby/Primitive.inl',
		'include/Moby/PrismaticJoint.h',
		'include/Moby/PSDeformableBody.h',
		'include/Moby/Quat.h',
		'include/Moby/RCArticulatedBodyFwdDynAlgo.h',
		'include/Moby/RCArticulatedBody.h',
		'include/Moby/RCArticulatedBodyDampening.h',
		'include/Moby/RCArticulatedBodyInvDynAlgo.h',
		'include/Moby/RecurrentForce.h',
		'include/Moby/RevoluteJoint.h',
		'include/Moby/RigidBody.h',
		'include/Moby/RigidBody.inl',
		'include/Moby/RNEAlgorithm.h',
		'include/Moby/RungeKuttaFehlbergIntegrator.h',
		'include/Moby/RungeKuttaFehlbergIntegrator.inl',
		'include/Moby/RungeKuttaImplicitIntegrator.h',
		'include/Moby/RungeKuttaImplicitIntegrator.inl',
		'include/Moby/RungeKuttaIntegrator.h',
		'include/Moby/RungeKuttaIntegrator.inl',
		'include/Moby/Simulator.h',
		'include/Moby/Simulator.inl',
		'include/Moby/SingularException.h',
		'include/Moby/SMatrix6.h',
		'include/Moby/SMatrix6N.h',
		'include/Moby/SpherePrimitive.h',
		'include/Moby/SphericalJoint.h',
		'include/Moby/SSL.h',
		'include/Moby/SSL.inl',
		'include/Moby/SSR.h',
		'include/Moby/SSR.inl',
		'include/Moby/StokesDragForce.h',
		'include/Moby/SVector6.h',
		'include/Moby/Tetrahedron.h',
		'include/Moby/TetraMeshPrimitive.h',
		'include/Moby/ThickTriangle.h',
		'include/Moby/Triangle.h',
		'include/Moby/TriangleMeshPrimitive.h',
		'include/Moby/Types.h',
		'include/Moby/UniversalJoint.h',
		'include/Moby/VariableEulerIntegrator.h',
		'include/Moby/VariableEulerIntegrator.inl',
		'include/Moby/Vector2.h',
		'include/Moby/Vector3.h',
		'include/Moby/VectorN.h',
		'include/Moby/VectorN.inl',
		'include/Moby/Visualizable.h',
		'include/Moby/XMLReader.h',
		'include/Moby/XMLTree.h',
		'include/Moby/XMLWriter.h']

# write out the package config file
#f = open('moby.pc', 'w')
#f.write('Name: Moby\n')
#f.write('Description: Moby Library\n')
#f.write('Version: 2.x\n')
#f.write('Libs: -Lmoby\n')
#f.write('Libs.private: ')
##for i in lib_path:
#  f.write('-L')
#  f.write(i)
#  f.write(' ')
#for i in libs:
#  f.write('-l')
#  f.write(i)
#  f.write(' ')
#f.write('\n') 
#f.write('Cflags: ')
#f.write(__CXXFLAGS)
#for i in inc_path:
#  f.write('-I')
#  f.write(i)
#  f.write(' ')
#f.write('\n') 
#f.close()

# Installation stuff
Alias("install", env.Install(INSTALL_PATH+'/include/Moby', headers))
Alias("install", env.Install(INSTALL_PATH+'/lib', 'lib/libmoby.so'))
Alias("install", env.Install(INSTALL_PATH+'/bin', binaries))
Alias("install", env.Install(INSTALL_PATH+'/lib/pkgconfig', 'moby.pc'))

