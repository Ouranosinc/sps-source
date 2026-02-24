cmake_minimum_required(VERSION 3.16)

message(STATUS "(EC) Generating spsdyn libraries Makefile")

#----- Append EC specific module path
list(APPEND ${CMAKE_SOURCE_DIR}/cmake_rpn/modules CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH})

ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}" LANGUAGES Fortran)
set(PROJECT_VERSION ${VERSION}${STATE})
set(spsdyn_VERSION ${PROJECT_VERSION} CACHE INTERNAL "spsdyn version" FORCE) # Needed for cascaded version identification
message(STATUS "(EC) ${PROJECT_NAME} version = ${PROJECT_VERSION}")

ec_build_info()            # Generate build information
include(ec_compiler_presets)

# To build without OpenMP, you have to add -DWITH_OPENMP=FALSE to the cmake command line
set(WITH_OPENMP TRUE CACHE BOOL "Control whether to use OpenMP")
include(ec_openmp)

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILER_FLAGS}")
endif()

file(GLOB F_FILES src/base/*.F90)

# Adding specific flags for sps
if (("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM") AND NOT ("${CMAKE_SYSTEM_NAME}" STREQUAL "CrayLinuxEnvironment"))
  # Flag for MKL is different depending on Intel compiler version
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 2021)
    set(MKL "-qmkl")
  else()
    set(MKL "-mkl")
  endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL} -static-intel")
set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -${MKL}")
  # flags in addition to extra checks defined in cmake_rpn compiler presets
  # https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/alphabetical-option-list.html
  if (EXTRA_CHECKS)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 2025)
      # maximum is 3 with Intel 2025
      set(OPT-REPORT "3")
    else()
      set(OPT-REPORT "5")
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -qopt-report=${OPT-REPORT} -C -init=snan,arrays -warn nointerfaces -check noarg_temp_created")
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  find_package(LAPACK)
  find_package(BLAS)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )

add_library(spsdyn STATIC ${F_FILES})

target_include_directories(spsdyn PRIVATE
   include
)
target_link_libraries(spsdyn PRIVATE
   modelutils
   modelutils_tmg_stubs
   modelutils_ov_ifort
   ${phy_LIBRARIES}
)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  if (EXTRA_CHECKS AND "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    target_compile_options(spsdyn PRIVATE "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:Fortran>>:-fcheck=bounds>" )
    message(STATUS "(EC) Extra debug flags for spsdyn: -fcheck=bounds")
  endif()
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/sps_version.in ${CMAKE_CURRENT_SOURCE_DIR}/include/sps_version.inc @ONLY)

message(STATUS "(EC) Generating mainspsdm Makefile")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/main/mainspsdm.in ${CMAKE_CURRENT_SOURCE_DIR}/src/main/mainspsdm.F90)

add_executable(mainspsdm src/main/mainspsdm.F90)

target_link_libraries(mainspsdm
   spsdyn
   vgrid::vgrid
   tdpack::tdpack
   rpncomm::rpncomm
   rmn::rmn-ompi
   OpenMP::OpenMP_Fortran
   MPI::MPI_Fortran
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)

target_link_libraries(modelutils PRIVATE modelutils_tmg_stubs)
target_link_libraries(rpnphy PRIVATE modelutils_ov_ifort)

install(TARGETS mainspsdm RUNTIME DESTINATION bin)
