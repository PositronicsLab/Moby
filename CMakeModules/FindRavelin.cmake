# Find Ravelin header and library.
#

# This module defines the following uncached variables:
#  RAVELIN_FOUND, if false, do not try to use Ravelin.
#  RAVELIN_INCLUDE_DIRS, where to find Ravelin/Ravelin_a.h.
#  RAVELIN_LIBRARIES, the libraries to link against to use the Ravelin library
#  RAVELIN_LIBRARY_DIRS, the directory where the Ravelin library is found.

find_path(
  RAVELIN_INCLUDE_DIR
  Ravelin/VectorNd.h 
  PATHS /usr/local/include /usr/include
)

if( RAVELIN_INCLUDE_DIR )
  find_library(
    RAVELIN_LIBRARY
    NAMES libRavelin Ravelin 
    PATHS /usr/local/lib /usr/lib
  )
  if( RAVELIN_LIBRARY )
    set(RAVELIN_LIBRARY_DIR "")
    get_filename_component(RAVELIN_LIBRARY_DIRS ${RAVELIN_LIBRARY} PATH)
    # Set uncached variables as per standard.
    set(RAVELIN_FOUND ON)
    set(RAVELIN_INCLUDE_DIRS ${RAVELIN_INCLUDE_DIR})
    set(RAVELIN_LIBRARIES ${RAVELIN_LIBRARY})
  endif(RAVELIN_LIBRARY)
else(RAVELIN_INCLUDE_DIR)
  message(FATAL_ERROR "FindRavelin: Could not find VectorNd.h")
endif(RAVELIN_INCLUDE_DIR)
	    
if(RAVELIN_FOUND)
  if(NOT RAVELIN_FIND_QUIETLY)
    message(STATUS "FindRavelin: Found both VectorNd.h and libRavelin.a")
  endif(NOT RAVELIN_FIND_QUIETLY)
else(RAVELIN_FOUND)
  if(RAVELIN_FIND_REQUIRED)
    message(FATAL_ERROR "FindRavelin: Could not find VectorNd.h and/or libRavelin.a")
  endif(RAVELIN_FIND_REQUIRED)
endif(RAVELIN_FOUND)
