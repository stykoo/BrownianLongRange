cmake_minimum_required(VERSION 2.8)

PROJECT(BrownianLongRange)

set(EXECUTABLE_NAME "BrownianLongRange")

# Build in release mode by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "-W -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "-Ofast -g -pg") # Profiling

# This is not a strict requirement
set(CMAKE_CXX_STANDARD 14)

# Packages needed
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Threads REQUIRED)
find_package(SFML REQUIRED COMPONENTS system window graphics)
include_directories(${HDF5_INCLUDE_DIRS})

# Executable
file(
	GLOB
	source_files
	src/*
)

add_executable(
	${EXECUTABLE_NAME}
	${source_files}
)

# Linking
target_link_libraries(
	${EXECUTABLE_NAME}
	${HDF5_LIBRARIES}
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	${SFML_LIBRARIES}
	${VTK_LIBRARIES}
)
