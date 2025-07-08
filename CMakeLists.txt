# CMake file for parallel MFEM projects
# Set cmake minimum version

cmake_minimum_required(VERSION 3.12.0...4.0.0)

# Set and check compiler
set(CMAKE_CXX_COMPILER mpic++)
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")

# Set C++ standard.
set(CMAKE_CXX_STANDARD 17)  # MFEM needs 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find mpi. mpic++ is set automatically when CMake finds mpi.
# find_package(MPI REQUIRED)

# Project name
project(mfem_projects)

# Path to MFEM build
set(MFEM "${CMAKE_SOURCE_DIR}/../../software/mfem")

# Path to HYPRE build
set(HYPRE "${CMAKE_SOURCE_DIR}/../../software/hypre")

# Path to METIES build
set(METIS "${CMAKE_SOURCE_DIR}/../../software/metis-5.1.0")

# Path to MFEMplus
set(MFEMPLUS "${CMAKE_SOURCE_DIR}/../../mfemplus")

# Path nlohmann-json
set(JSON "/opt/homebrew/Cellar/nlohmann-json/3.12.0/include")

# Get all cpp files in src directory
file(GLOB SRC "${CMAKE_SOURCE_DIR}/source/*.cpp")

# Loop over all .cpp files
foreach(src ${SRC})

	# Get filename without extension
	 get_filename_component(build_name ${src} NAME_WE)


	# Add source for the build
	add_executable(${build_name} ${src} 
				     ${MFEMPLUS}/custombilininteg.cpp
				     ${MFEMPLUS}/customfunctions.cpp )
	
	# Include headers
	target_include_directories(${build_name} PRIVATE ${MFEM}
							# ${MFEM}/include
						     ${HYPRE}/src/hypre/include
						     ${HYPRE}/src/utilities
						     ${METIS}/include
						     ${MFEMPLUS}
					         ${JSON}
					         # ${MPI_CXX_INCLUDE_DIRS}
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
						    # ${MPI_CXX_LIBRARIES}
						    )

	message(STATUS "Configured target: ${build_name}")

endforeach()
