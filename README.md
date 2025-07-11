# Intro to MFEM

This is a guide to smoothly install and run finite element codes in C++ using the **MFEM** library. **MFEM** is a powerful and vast library. **Do not be intimidated by it!!** With a few weeks, you can get comfortable running your own finite element codes. If you don't know what a function or a class does, use command + left click to enter the header file where it is defined. See what it does, what constructors it has, what arguments it takes, and so on. You may have to go through several layers of files to understand. Your last resort is to just try it and output the results. 

## Coding in C++

C++ is an object oriented language. You create 'objects', which are instances of classes, that have 'methods', which are functions within a class. Every relevant piece of a finite element code is usually an object in **MFEM**, for example the finite element space, the type of element, the mesh, and so on. Therefore, it is necessary to have a basic acquaintance at least with classes and functions in C++, and as you code you will learn more.
This was the resource I used: 
> https://www.w3schools.com/cpp/default.asp

## Git

Git is the version control that you should use for your projects. It is not strictly necessary to code, but it is good practice.
This was the resource I used:
> https://www.w3schools.com/git/default.asp

Git add, commit, push, fetch, merge, and status are some common commands.

## Setting up MFEM

All of the following can be done on terminal. Set up your directories in the same way mentioned here. It will make it easier to integrate MFEM and MFEMPLUS (custom codes written for AML) into your codes and to share between AML. Create a clean new directory named "MFEM". Everything will be done installed there. Create the following subdirectories. 

- MFEM
	- software
	- projects
	- supportinglibraries

Enter the projects directory and clone this GitHub repo.

### Building MFEM

MFEM can be built both in serial and parallel. Serial MFEM codes will only run on one processor. However, for any practical large scale simulations, we need to run parallelized codes, therefore we will build MFEM parallel. However, the parallel build can also compile serial codes.
Enter the *software* directory. It should be empty. The following assumes that you have git installed on your system. Open your terminal or your favorite terminal emulator. Wezterm and iTerm work well. Run the following commands.

#### Step 1: Install C++ compiler and CMake
Check if you have a C++ compiler installed

	gcc --version
	make --version

If you get an error, run 

	xcode-select --install

A dialog box will appear, click on Install. This will install necessary compilers and build tools.

Install CMake and Make. CMake is a build system generator. It produces the makefiles, which can be run by make. Make is the actual build system, which produces the executable that you can run. Make sure Make and CMake are installed on your system. Run the following commands.

	brew install make
	make --version
	brew install cmake
	cmake --version


#### Step 2: Download and build Metis

Go back to the software directory.

	curl -sL https://github.com/mfem/tpls/raw/refs/heads/gh-pages/metis-5.1.0.tar.gz -o metis-5.1.0.tar.gz
	tar -xzf metis-5.1.0.tar.gz
	cd metis-5.1.0

Metis uses CMake to build. However, in the donwloaded Metis folder, the minimum CMake version specified is older than what MFEM uses, so we need to change that. In the current folder, open *CMakeLists.txt* in Sublime Text. In the first line, change 
> "cmake_minimum_required(VERSION 2.8)" 

to 

> "cmake_minimum_required(VERSION 3.12.0...4.0.0)"

Now we can compile the Metis library. 

	make BUILDDIR=lib config
	make BUILDDIR=lib -j8

#### Step 3: Download and build Hypre

Hypre is a library of parallel sparse arrays and linear solvers.
Go back to the software directory and do the following.

	curl -sL https://github.com/hypre-space/hypre/archive/refs/tags/v2.26.0.tar.gz -o hypre-2.26.0.tar.gz
	tar -xzf hypre-2.26.0.tar.gz
	cd hypre-2.26.0/src
	./configure --disable-fortran
	make -j8
	cd ../..
	ln -s hypre-2.26.0 hypre

#### Step 4: Clone and build MFEM

We built Hypre and Metis first because the parallel build of MFEM requires these libraries.  
Clone MFEM.

	git clone https://github.com/mfem/mfem.git
	cd mfem
	git checkout v4.8

This clones the git repository onto your local machine and switches to v4.8 branch. 
Now build it.

	make parallel -j8 MFEM_USE_MPI=YES MFEM_USE_METIS_5=YES METIS_DIR=@MFEM_DIR@/../metis-5.1.0

The build is complete. The mfem folder must contain libmfem.a, which is a compiled library of the MFEM code that our C++ codes will use.

### Install Paraview

Run 

	brew install paraview

Paraview will be used to visualize the results.

### Clone MFEMPLUS

Clone MFEMPLUS from GitHub. Go back to the software directory under the top level directory and run:

	git clone https://github.com/akulakar/mfemplus.git

It should create a new directory titled mfemplus, which should contain hpp and cpp files. You can link and use these in your code.

### Install supporting libraries

Run

	brew install nlohmann-json

This installs a library to support json objects in C++. json scripts will be used to store simulation parameters.
Then enter the supporting libraries directory and run

	git clone https://github.com/ArashPartow/exprtk.git

This library gives support to evaluate certain expressions. 

### Working on an MFEM project and executing files

Enter the projects directory and the MFEMSetUpAML directory that was cloned from GitHub.

#### 1. Subdirectories

MFEMSetUpAML has the following sub directories.

- source
- scripts
- mesh
- .vscode
- .git

It also contains CMakeLists.txt and .gitignore.

Create three more subdirectories:

	mkdir drivers results build


*source* contains "isotropic_elasticity.cpp" and "isotropic_elasticity_prl.cpp". These are the two cpp files that we will want to build and execute.  *drivers* will contain supporting cpp files. For this project it is empty.  *scripts* contains a json file to store simulation parameters and mathematica file to generate meshes. You can also have mathematica scripts to visualize certain results or shell scripts to execute certain commands, and so on. *build* will contain the executables. The executables are the compiled code that will actually be run by the computer. results will contain results such as dispalecements, stressesgit w, etc., and mesh will contain meshes.

#### 2. .vscode settings

The c_cpp_properties.json is important to get the relevant mfem and mfemplus libraries included. This is only for .vscode to know where these functions and classes are for the sake of making the .vscode experience smoother, not for compilation. tasks.json builds the code and launch.json exectutes it. However, we will build with cmake and make, and execute with terminal. So we will not rely on vscode for the build. Check the c_cpp_properties.json in the .vscode folder to see how to include paths, it is straightforward.

#### 3. CMakeLists.txt

There is a CMakeLists.txt file in the directory. It should work as is, and for any new project in the projects directory. The easiest way I found to learn CMake and how to use it was to ask ChatGPT. The details are not provided here. The general idea is that CMakeLists.txt will specify the compiler (mpicc for us), the relevant libraries (metis, hypre, and mfem) and cpp codes to be compiled. The important command to run is the following. Enter the build folder and run

	cmake ..

This command generates the makefiles in the build folder. Now the make files are generated to actually build your executables from the cpp codes. 

#### 4. Build with make and execute

Let us first build the serial code. Run:

	make isotropic_elasticity

Now the code is compiled with mpicc. This code is meant to be run on one processor. If you run it on more than one, it will just run the same code one each processor, it will not run a single parallelized executable. Execute the code with either of the following:

	mpirun -np 1 isotropic_elasticity

or

	./isotropic_elasticity

Next, let us build the parallel code. The parallel code can be run on one processor or on as many as your machine has. Generally, for a simple computation, the time taken will not greatly reduce when more than one core is used. In fact, it may increase. However, with larger computations, parallelized codes make a big difference.

	make isotropic_elasticity_prl

Execute with:

	mpirun -np <num-of-cores> isotropic_elasticity_prl


Now you should have serial and parralel mfem codes compiling without issues. The best way to learn what everything in this directory does, from the headers that are included, to the commands being used, are to first try everything out yourself. Change certain things and see what happens. Change the boundary conditions, the mesh, the element order, etc. If that doesn't help, ask ChatGPT what something does. 

## Next steps

This was not meant to be an exhaustive tutorial for MFEM by any means. In fact, it can feel quite overwhelming to get a hang of it. However, I PROMISE YOU it will not take longer than 2 weeks and an attempt at coding up something yourself to get acquainted with it. After that, it becomes fun!

