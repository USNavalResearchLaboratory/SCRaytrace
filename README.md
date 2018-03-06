# SCRaytrace
Solar Corona Raytracing tools

## Project tree structure
- src: C++ source code.
- data: data files used by the C++ code.
- python: python wrappers. Empty for now...
- idl: IDL wrappers. Empty for now. We will see if we transfer the SolarSoft sources in here.
- docbook: Docbook source files for the user manual.

## Compiling the sources
SCraytrace now uses cmake. To compile the project once cloned from github, run the following:  
cmake -H. -Bbuild  
cmake --build build -- -j3

The compiled code will be generated in the build folder.  

For additional testing, you can run the following:  
cd build  
make test

## Dependencies
- boost C++ library
- cppunit

## Generating the code documentation with Doxygen
Use doxywizard to generate the code documentation. The Doxyfile file is the configuration file.  

## Gererating the user manual
The user manual is writen using the Docbook.  

## To Do
- Implement ray-tracing of data cube from predictive science
- Implement SoloHI and WISPR projections
- Implement a Python wrapper
- Documentation for the F-Corona/VSF models
- Generate pakage with cmake

