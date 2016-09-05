cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(oilprj)

#include_directories( ../../include )

# set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# set the build type
set(CMAKE_BUILD_TYPE Release)
# confirmation messages
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")

find_package(PkgConfig)
pkg_check_modules(OPENCV REQUIRED opencv)

add_executable(oiledge-detector ./src/oiledge-detector.cpp ./src/OillineDetector.cpp ./src/tool.cpp)
target_link_libraries(oiledge-detector ssm /usr/local/lib/kaw_lib.a ${OPENCV_LIBRARIES})

add_executable(oiledge-viewer ./src/oiledge-viewer.cpp) 
target_link_libraries(oiledge-viewer ssm)

add_executable(oiledge-sample ./src/oiledge-sample.cpp ./src/OillineDetector.cpp ./src/OillineDetector2.cpp ./src/tool.cpp)
target_link_libraries(oiledge-sample ssm /usr/local/lib/kaw_lib.a ${OPENCV_LIBRARIES})


