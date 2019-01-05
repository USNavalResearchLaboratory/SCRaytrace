# SCRaytrace
Solar Corona Raytracing tools

## Project tree structure
- src: C++ source code.
- data: data files used by the C++ code.
- python: python wrappers. Empty for now...
- idl: IDL wrappers. Empty for now. We will see if we transfer the SolarSoft IDL sources in here maybe, which would allow collaborative development, things that was not possible so far.
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
The user manual is writen using Docbook. More details will follow here on how to edit and generate the manual.  



