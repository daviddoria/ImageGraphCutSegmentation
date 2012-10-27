enable_testing()

ADD_EXECUTABLE(ImageGraphCutSegmentationTest ImageGraphCutSegmentationTest.cpp)
TARGET_LINK_LIBRARIES(ImageGraphCutSegmentationTest
${ITK_LIBRARIES}
${ImageGraphCutSegmentation_libraries})

add_test(NAME ImageGraphCutSegmentationTest COMMAND ImageGraphCutSegmentationTest
${CMAKE_SOURCE_DIR}/data/soldier.png ${CMAKE_SOURCE_DIR}/data/foreground.png ${CMAKE_SOURCE_DIR}/data/background.png output.png)
