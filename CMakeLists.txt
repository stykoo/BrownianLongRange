cmake_minimum_required(VERSION 2.8)

PROJECT(BrownianLongRange)

# Build in release mode by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "-W -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -pg")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	set(CMAKE_CXX_FLAGS "-W -Wall -axSSE4.1,AVX,CORE-AVX2,CORE-AVX512")
endif()

# This is not a strict requirement
set(CMAKE_CXX_STANDARD 14)

# Packages needed
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules" ${CMAKE_MODULE_PATH})
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Threads REQUIRED)
find_package(SFML REQUIRED COMPONENTS system window graphics)
find_package(MKL)
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${MKL_INCLUDE_DIRS})

if (MKL_FOUND)
	set(EXECUTABLE_NAME "BrownianLongRange_MKL")
else()
	set(EXECUTABLE_NAME "BrownianLongRange")
endif()


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

if (MKL_FOUND)
	target_compile_definitions(${EXECUTABLE_NAME} PRIVATE USE_MKL)
endif()

# Linking
target_link_libraries(
	${EXECUTABLE_NAME}
	${HDF5_LIBRARIES}
	${GSL_LIBRARIES}
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	${SFML_LIBRARIES}
)

target_link_libraries(${EXECUTABLE_NAME} -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group pthread dl)
