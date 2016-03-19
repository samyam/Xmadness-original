# NOTE: There is not a consistent ABI for Fortran's LOGICAL. Though it is
#       standard for it to be represented in the same manner as INTEGER,
#       the value for .true. is often +1 but sometimes -1, whereas .false.
#       is almost always 0. We will thus standardize on .true.=+1 and .false=0,
#       though this is known to occasionally break. Unfortunately it is
#       nontrivial to determine this value since the same Fortran compiler
#       which compiled ScaLAPACK would need to be known in order to create
#       a custom routine which sets the input to .true. or .false.
if(NOT ELEM_FORT_LOGICAL)
  set(ELEM_FORT_LOGICAL int)
endif()
if(NOT ELEM_FORT_TRUE)
  set(ELEM_FORT_TRUE 1)
endif()
if(NOT ELEM_FORT_FALSE)
  set(ELEM_FORT_FALSE 0)
endif()

# Go ahead and check for Fortran, but keep in mind that CMake's 'OPTIONAL' 
# argument for enable_language is still completely broken as of 2.8.8
workaround_9220(Fortran FORTRAN_WORKS)
if(FORTRAN_WORKS)
  enable_language(Fortran OPTIONAL)
  set(ELEM_HAVE_F90_INTERFACE FALSE)
  if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    include(FortranCInterface)
    FortranCInterface_VERIFY(CXX)
    if(FortranCInterface_VERIFIED_CXX)
      set(ELEM_HAVE_F90_INTERFACE TRUE)
      FortranCInterface_HEADER(
        ${CMAKE_CURRENT_BINARY_DIR}/include/elemental/FCMangle.h 
        MACRO_NAMESPACE "FC_")
      install(FILES ${PROJECT_BINARY_DIR}/include/elemental/FCMangle.h
              DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/madness/external/elemental/"
              COMPONENT elemental)
    endif()
  else()
    message(STATUS "${CMAKE_Fortran_COMPILER} does not appear to support F90")
  endif()
else()
  message(STATUS "Could not find working Fortran compiler")
endif()
