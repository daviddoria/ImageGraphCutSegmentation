cmake_minimum_required(VERSION 2.8)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

Project(ImageGraphCutSegmentation)

# Compile with c++11 functionality.
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++11") # gcc >= 4.7
endif(UNIX)

# ITK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                          ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                          ITKBinaryMathematicalMorphology ITKDistanceMap)
INCLUDE( ${USE_ITK_FILE} )

if( "${ITK_VERSION_MAJOR}" LESS 4 )
  MESSAGE(FATAL_ERROR "You must build ImageGraphCutSegmentation with ITK >= 4.0!")
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

add_subdirectory(Kolmogorov)

# Mask submodule
get_property(FoundMask GLOBAL PROPERTY MaskIncludeDirs SET)
if(NOT FoundMask)
  add_subdirectory(Mask)
endif()

get_property(MaskIncludeDirs GLOBAL PROPERTY MaskIncludeDirs)
set(ImageGraphCutSegmentation_include_dirs ${ImageGraphCutSegmentation_include_dirs} ${MaskIncludeDirs})
get_property(MaskLibs GLOBAL PROPERTY MaskLibs)
set(ImageGraphCutSegmentation_libraries ${ImageGraphCutSegmentation_libraries} ${MaskLibs})

# Create graph library
add_library(ImageGraphCut ImageGraphCut.cxx)
target_link_libraries(ImageGraphCut KolmogorovMaxFlow)
set(ImageGraphCutSegmentation_libraries ${ImageGraphCutSegmentation_libraries} ImageGraphCut) # Add the library to the list of libraries used

# Allow this project to be detected and used as a submodule
set(ImageGraphCutSegmentation_include_dirs ${ImageGraphCutSegmentation_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(GLOBAL PROPERTY ImageGraphCutSegmentationIncludeDirs ${ImageGraphCutSegmentation_include_dirs})
set_property(GLOBAL PROPERTY ImageGraphCutSegmentationLibs ${ImageGraphCutSegmentation_libraries})

# Give the compiler all of the required include directories
include_directories(${ImageGraphCutSegmentation_include_dirs})

option(ImageGraphCutSegmentation_BuildExamples "Build ImageGraphCutSegmentation examples?")
if(ImageGraphCutSegmentation_BuildExamples)
  add_subdirectory(Examples)
endif()

option(ImageGraphCutSegmentation_BuildTests "Build ImageGraphCutSegmentation tests?")
if(ImageGraphCutSegmentation_BuildTests)
  add_subdirectory(Tests)
  add_subdirectory(Kolmogorov/Tests)
endif()
