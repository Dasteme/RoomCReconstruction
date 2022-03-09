cmake_minimum_required(VERSION 3.20)
project(RoomCReconstruction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)

set(HEADER_FILES
        external/tinyply/tinyply.hpp
        external/nanoflann/nanoflann.hpp
        external/ts/pc/parallel.hpp
        external/ts/pc/pc_io.hpp
        external/ts/pc/pc_tools.hpp
        external/ts/pc/pc_base.hpp
        external/ts/pc/pc_score.hpp
        src/wallextender.hpp
        src/TriangleNode3D.hpp
        src/rcr_io.hpp
        src/Helper.hpp
        )

set(SOURCE_FILES
        external/tinyply/tinyply.cpp
        external/ts/pc/pc_io.cpp
        external/ts/pc/pc_tools.cpp
        external/ts/pc/pc_score.cpp
        src/main.cpp
        src/wallextender.cpp
        src/TriangleNode3D.cpp
        src/rcr_io.cpp
        src/Helper.cpp
        )



# Create target for compiling and running
set(EXECUTABLE_NAME RoomCReconstruction)
add_executable(${EXECUTABLE_NAME} ${HEARDER_FILES} ${SOURCE_FILES})




target_include_directories(
        ${EXECUTABLE_NAME}
        PRIVATE include
        PRIVATE external/ts
        PRIVATE external/parallel
        PRIVATE SYSTEM external)


target_link_libraries(
        ${EXECUTABLE_NAME}
        PRIVATE Eigen3::Eigen
        PRIVATE Threads::Threads)