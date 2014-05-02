cmake_minimum_required(VERSION 2.8)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

project(ImageGraphCutSegmentation)
enable_testing()

add_subdirectory(CMakeHelpers)

# Compile with c++11 functionality.
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++11") # gcc >= 4.7
endif(UNIX)

# ITK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                          ITKImageIntensity ITKImageFeature
                          ITKMathematicalMorphology
                          ITKBinaryMathematicalMorphology ITKDistanceMap
                          ITKTestKernel)
INCLUDE( ${USE_ITK_FILE} )

if( "${ITK_VERSION_MAJOR}" LESS 4 )
  MESSAGE(FATAL_ERROR "You must build ImageGraphCutSegmentation with ITK >= 4.0!")
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

# Boost
set(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE(Boost 1.42 COMPONENTS required)

INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS})

# Submodules
UseSubmodule(Mask ImageGraphCutSegmentation)

# Build the graph cut library
add_subdirectory(Kolmogorov)
set(ImageGraphCutSegmentation_libraries ${ImageGraphCutSegmentation_libraries} KolmogorovMaxFlow)

# Allow this project to be detected and used as a submodule
CreateSubmodule(ImageGraphCutSegmentation)

# Give the compiler all of the required include directories
include_directories(${ImageGraphCutSegmentation_include_dirs})

# Make the h/hpp files appear in a QtCreator project
add_custom_target(ImageGraphCut SOURCES
ImageGraphCut.h ImageGraphCut.hpp PixelDifference.h)

option(ImageGraphCutSegmentation_BuildExamples "Build ImageGraphCutSegmentation examples?")
if(ImageGraphCutSegmentation_BuildExamples)
  add_subdirectory(Examples)
endif()

option(ImageGraphCutSegmentation_BuildTests "Build ImageGraphCutSegmentation tests?")
if(ImageGraphCutSegmentation_BuildTests)
  add_subdirectory(Tests)
  add_subdirectory(Kolmogorov/Tests)
endif()
