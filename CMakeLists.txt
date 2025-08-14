# CMake file for parallel MFEM projects
# Set cmake minimum version

cmake_minimum_required(VERSION 3.12.0...4.0.0)
# Project name
project(mfem_projects LANGUAGES CXX)


# Set and check compiler
# set(CMAKE_CXX_COMPILER mpic++)
# message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")

# Set C++ standard.
set(CMAKE_CXX_STANDARD 17)  # MFEM needs 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find mpi. mpic++ is set automatically when CMake finds mpi.
find_package(MPI REQUIRED)

# Tell CMake how to use OpenMP on macOS with AppleClang
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/include")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "/opt/homebrew/Cellar/libomp/20.1.8/lib/libomp.dylib")

# Find openmp.
find_package(OpenMP REQUIRED)

# Only set default if not provided
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose Debug or Release" FORCE)
endif()


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# Path to MFEM build
set(MFEM "${CMAKE_SOURCE_DIR}/../../software/mfem")

# Path to HYPRE build
set(HYPRE "${CMAKE_SOURCE_DIR}/../../software/hypre")

# Path to METIES build
set(METIS "${CMAKE_SOURCE_DIR}/../../software/metis-5.1.0")

# Path to MFEMPLUS
set(MFEMPLUS "${CMAKE_SOURCE_DIR}/../../software/mfemplus")

# Path nlohmann-json
set(JSON "/opt/homebrew/Cellar/nlohmann-json/3.12.0/include")

# Path to OpenMP
set(OPENMP "/opt/homebrew/Cellar/libomp/20.1.8/include/")

# Get all cpp files in drivers directory
file(GLOB DRV "${CMAKE_SOURCE_DIR}/source/*.cpp")

# Loop over all .cpp files
foreach(drv ${DRV})

	# Get filename without extension
	 get_filename_component(build_name ${drv} NAME_WE)


	# Add source for the build
	add_executable(${build_name} ${drv} 
				     ${MFEMPLUS}/custombilininteg.cpp
				     ${MFEMPLUS}/customfunctions.cpp 
				     ${MFEMPLUS}/customoperators.cpp
				     ${MFEMPLUS}/computequantities.cpp)
	
	# Include headers
	target_include_directories(${build_name} PRIVATE ${MFEM}
							# ${MFEM}/include
						    ${HYPRE}/src/hypre/include
					      ${HYPRE}/src/utilities
					      ${METIS}/include
					      ${MFEMPLUS}
				        ${JSON}
				        ${MPI_CXX_INCLUDE_DIRS}
				        ${OpenMP_CXX_INCLUDE_DIRS}
								${OPENMP}
						    include		   					
								)

	# Link directory to search for libmfem.a, libhypre.a, libmetis.a
	target_link_directories(${build_name} PRIVATE ${MFEM}
					   	    ${HYPRE}/src/lib
						    	${METIS}/lib/libmetis 
						    	)

	# Link libraries to the executable
	target_link_libraries(${build_name} PRIVATE mfem
						    pthread
						    dl
						    HYPRE	
						    metis
						    ${MPI_CXX_LIBRARIES}
    						OpenMP::OpenMP_CXX
						    )

	message(STATUS "Configured target: ${build_name}")

endforeach()
