# CMake file for MFEM projects.

cmake_minimum_required(VERSION 3.12.0...4.0.0)

# Use C++17, as required by MFEM
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Project name
project(isotropic_elasticity)

# Path to MFEM serial build
set(MFEM_serial "${CMAKE_SOURCE_DIR}/../../serial/mfem")

# Path to mfemplus
set(CUSTOM_DIR "${CMAKE_SOURCE_DIR}/../../mfemplus")

# Get all .cpp files in drivers
file(GLOB DRIVERS "${CMAKE_SOURCE_DIR}/drivers/*.cpp")

# Loop over all .cpp files

foreach(driver ${DRIVERS})

	# Get filename without extension
	get_filename_component(build_name ${driver} NAME_WE)

	# Add the source file for the build
	add_executable(${build_name}	${driver}
								 	${CUSTOM_DIR}/custombilininteg.cpp
								 	${CUSTOM_DIR}/customfunctions.cpp )

	# Include headers
	target_include_directories(${build_name} PRIVATE
		${MFEM_serial}
	    ${CUSTOM_DIR}
	    /opt/homebrew/Cellar/nlohmann-json/3.12.0/include
	)

	# Link directory to search for libfem.a
	target_link_directories(${build_name} PRIVATE ${MFEM_serial})

	# Link MFEM to executable
	target_link_libraries(${build_name} 	PRIVATE mfem
											pthread
											dl)

	message(STATUS "Configured target: ${build_name}")

endforeach()