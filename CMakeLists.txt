#=============================================================================
# AdH - Adaptive Hydraulics
# General executable
#=============================================================================

##############################################################################
# Define the source files needed to build Pre-AdH
#list(APPEND adh_sources
#    grid_interpolator.c
#    read_adh_grid.c
#    )

# Build executable
add_executable(grid_interpolator ${adh_sources})

# Cygwin
  if(CYGWIN AND CMAKE_COMPILER_IS_GNUCC)
    # If compiling on Windows in cygwin v3.x then we need the -mno-cygwin flag
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
      OUTPUT_VARIABLE GCC_VERSION)
    if(GCC_VERSION VERSION_LESS 4.0) # no longer supported in GCC v4
      set(compiler_flags "${compiler_flags} -mno-cygwin")
    endif()
  endif(CYGWIN AND CMAKE_COMPILER_IS_GNUCC)

  target_link_libraries(grid_interpolator
    #    adh_fe
    # adh_initio
    # adh_structs
    #adh_columns
    #adh_initio
    #adh_solver
    #adh_structs
    #adh_fe
    #adh_meteor
    #adh_grid
    #adh_node
    #adh_elem
    #adh_friction
    #adh_tools
    #adh_structs
    #adh_debug
    #adh_messg
    #adh_meteor 
      )

  #add_definitions(-D_DEBUG_KITCHEN_SINK -D_DEBUG_TEMPORARY -D_DEBUG -g)

##############################################################################
# Check build type (prevent misspelling)
if(NOT CMAKE_BUILD_TYPE MATCHES "^(|Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(FATAL_ERROR "CMAKE_BUILD_TYPE parameter should be left empty or set to Debug (-D_DEBUG -g), Release, or RelWithDebInfo (-g).")
endif(NOT CMAKE_BUILD_TYPE MATCHES "^(|Debug|Release|RelWithDebInfo|MinSizeRel)$")
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-D_DEBUG -g)
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
  add_definitions(-g)
endif(CMAKE_BUILD_TYPE MATCHES "Debug")

# Build options for
# Determine what level of debugging we are using
set(BUILD_DEBUG_LEVEL "0" CACHE STRING
  "0: None  1: -g  2: -g + _DEBUG  3: Level 2 + Temporary  4: Level 2 + Kitchen Sink  5: All Debugging on")
# Determine if we are using some debugging so that debuggers can show the code
if(BUILD_DEBUG_LEVEL LESS 0 OR BUILD_DEBUG_LEVEL GREATER 5)
  message(FATAL_ERROR "BUILD_DEBUG_LEVEL parameter must be an integer value of the range 0 to 5, inclusive.")
elseif(BUILD_DEBUG_LEVEL EQUAL 1)
  add_definitions(-g)
# Determine if we are using the memory checking debugging
elseif(BUILD_DEBUG_LEVEL EQUAL 2)
  add_definitions(-D_DEBUG -g)
# Determine if we are using temporary user-defined debug functions
elseif(BUILD_DEBUG_LEVEL EQUAL 3)
  add_definitions(-D_DEBUG_TEMPORARY -D_DEBUG -g)
# Determine if we are using full-check-everything debugging
elseif(BUILD_DEBUG_LEVEL EQUAL 4)
  add_definitions(-D_DEBUG_KITCHEN_SINK -D_DEBUG -g)
# Determine if we are using all debugging potential
elseif(BUILD_DEBUG_LEVEL EQUAL 5)
  add_definitions(-D_DEBUG_KITCHEN_SINK -D_DEBUG_TEMPORARY -D_DEBUG -g)
endif(BUILD_DEBUG_LEVEL LESS 0 OR BUILD_DEBUG_LEVEL GREATER 5)
# Check whether the build type and debug level are consistent
if((BUILD_DEBUG_LEVEL EQUAL 0 AND CMAKE_BUILD_TYPE MATCHES "^(Debug|RelWithDebInfo)$")
OR (BUILD_DEBUG_LEVEL GREATER 0 AND CMAKE_BUILD_TYPE MATCHES "^(Release|MinSize)$"))
  message(FATAL_ERROR "BUILD_DEBUG_LEVEL is inconsistent with CMAKE_BUILD_TYPE; both parameters should include or exclude debug definitions.")
endif((BUILD_DEBUG_LEVEL EQUAL 0 AND CMAKE_BUILD_TYPE MATCHES "^(Debug|RelWithDebInfo)$") OR
   (BUILD_DEBUG_LEVEL GREATER 0 AND CMAKE_BUILD_TYPE MATCHES "^(Release|MinSize)$"))

  # XDMF (Kitware/Paraview Binary)
    target_link_libraries(grid_interpolator
        ${BLAS_LIBRARIES}
        ${XML_LIBRARY}
        ${HDF5_LIBRARY}
        ${SZIP_LIBRARY}
        ${ZLIB_LIBRARY}
      )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -ldl -llzma")
      message("HDF5 INCLUDE DIR: ${HDF5_INCLUDE_DIR}")
      message("HDF5 libraries: ${HDF5_LIBRARY}")
      message("XML INCLUDE DIR: ${XML_INCLUDE_DIR}")
      message("XML libraries: ${XML_LIBRARY}")

      #  if(PARMETIS_FOUND)
      #target_link_libraries(grid_interpolator ${PARMETIS_LIBRARIES})
      #endif(PARMETIS_FOUND)

##############################################################################
#include_directories("${CMAKE_SOURCE_DIR}/tools/STAND_ALONE_UTILITIES/grid_interpolator")
#include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/debug/include")
