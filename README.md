# Intro to MFEM

This is a guide to smoothly install and run finite element codes in C++ using the **MFEM** library. **MFEM** is a powerful and vast library. **Do not be intimidated by it!!** With a few weeks, you can get comfortable running your own finite element codes. If you don't know what a function or a class does, see what it does, what constructors it has, what arguments it takes, and so on. You may have to go through several layers of files to understand. Your last resort is to try it and see what happens. 

## Coding in C++

C++ is an object oriented language. You create 'objects' that have 'methods'. Every relevant piece of a finite element code is usually an object in **MFEM**, for example the finite element space, the type of element, the mesh, and so on. Therefore, it is necessary to have a basic acquaintance at least with classes and functions in C++, and as you code you will learn more.
This was the resource I used: 
> https://www.w3schools.com/cpp/default.asp

## Git

Git is the version control that you should use for your projects. It is not strictly necessary to code, but it is good practice.
This was the resource I used:
> https://www.w3schools.com/git/default.asp

You can learn enough git to get started in about 1 - 2 hrs. 

## Setting up MFEM serial

The following instructions are for you to create a new directory and do the following. Not to work on the cloned directory itself. However, the cloned directory will match everything stated here.
Set up your directories in the following way. It will make it easier to integrate MFEM and MFEMPLUS (custom codes written for AML) into your codes. Create a clean new directory named "MFEM". Everything will be done there. All of the following can be done on terminal. 

- MFEM
	- projects
	- serial
	- parallel
	- supportinglibraries

### Building MFEM serial

MFEM serial is the serial version of MFEM, so designed to only run on one processor. It is easier to set up and get started with. but to run heavy computations you need to use the parallel version. However, this is a good first start.

Enter the *serial* directory under the *MFEM* directory. It should be empty. The following assumes that you have git installed on your system. Open your terminal or your favorite terminal emulator. Wezterm and iTerm are popular. Run the following commands.

#### Step 1: Install C++ compiler
Check if you have a C++ compiler installed

	gcc --version
	make --version

If you get an error, run 

	xcode-select --install

A dialog box will appear, click on Install. This will install necessary compilers and build tools.

#### Step 2: Download MFEM

	git clone https://github.com/mfem/mfem.git
	cd mfem
	git checkout v4.8

This copies the git repository onto your local machine and switches to v4.8 branch.

#### Step 3: Build MFEM serial

	make serial -j8

The build is complete. The mfem folder must contain libmfem.a, which is a compiled library of the MFEM code that our C++ codes will use.

### Download MFEMPLUS

Clone MFEMPLUS from GitHub. Exit the serial directory. You must be in your overall MFEM directory. Run the following:

	git clone https://github.com/akulakar/mfemplus.git

It should create a new directory titled mfemplus, which should contain hpp and cpp files. You can link and use these in your code.


### Install supporting libraries

Run

	brew install nlohmann-json

This installs a library to support json objects in C++. json scripts will be used to store simulation parameters.
Then enter the supporting libraries directory and run

	git clone https://github.com/ArashPartow/exprtk.git

This library gives support to evaluate certain expressions. 

### Setting up the projects directory

Enter the projects directory. Create a new project. Name it for example, "BeamBending", or anything else you want.

#### Step 1: Create directories

Enter BeamBending. Create the following subdirectories:

- drivers
- source
- scripts
- build
- results
- mesh

drivers will contain supporting cpp files. source will contain the main cpp file you will build and execute. scripts can contain, for example, json files to store simulation parameters, mathematica files to visualize certain results or create meshes, and so on. build will contain the executable. This is the compiled code that will actually be run by your computer. results will contain results, and mesh will contain meshes you want to use.

#### Step 2: .vscode settings

The c_cpp_properties.json is important to get the relevant mfem and mfemplus libraries included. This is only for .vscode to know where these functions and classes are for the sake of making the .vscode experience better, not for compilation. tasks.json builds the code and launch.json exectutes it. However, we will build with cmake and make, and execute with terminal. So we will not use these. Check the c_cpp_properties.json to see how to include paths, it is straightforward.

#### Step 3: CMakeLists.txt

CMake is a build system generator. It produces the makefiles, which can be run by make. make is the actual build system, which produces the executable that you can run. Make sure make and cmake are installed on your system. Run the following commands.

	brew install make
	gmake --version
	brew install cmake
	cmake --version

There is a CMakeLists.txt file in the directory. It should work as it is in your project too. The easiest way to learn what it does it to ask ChatGPT what it does. The important command to run is the following. Enter build and run

	cmake ..

Note the space. Now the make files are generated to actually build your cpp files. 

#### Step 4: Build

To build the cpp file into an executable that can be run by your computer, run

	make <name-of-cpp-file-in-source-folder>

If it is successful, you can execute the code with 

	./<name-of-cpp-file-in-source-folder>

Now you should have mfem serial codes compiling without issues. The best way to learn what everything in this directory does, from the headers that are included, to the commands being used, are to first try everything out yourself. Change certain things and see what happens. If that doesn't help, ask ChatGPT what something does. 


## Next steps

Once you have some practice with MFEM serial, build and use MFEM parallel.