# SCRaytrace
Solar Corona Raytracing tools

## Project tree structure
- src: C++ source code.
- data: data files used by the C++ code.
- python: python wrappers.
- idl: IDL wrappers. Empty for now. We will see if we transfer the SolarSoft IDL sources in here maybe, which would allow collaborative development, things that was not possible so far.
- docbook: Docbook source files for the user manual. [Deprecated. Use Sphinx now. Migration ongoing.]
- sphinx: sphinx documentation
- doxygen: doxygen code documentation

## Cloning the repository
`git clone https://github.com/DrRaytrace/SCRaytrace.git`

## Compiling the sources
SCraytrace now uses cmake. To compile the project once cloned from github, run the following:  
`cmake -H. -Bbuild`  
`cmake --build build -- -j3`

The compiled code will be generated in the build folder.  

For additional testing, you can run the following:  
`cd build`  
`make test`

Individual test:  
`cd build/src`  
`./testboosttest --log_level=message`, for example.

## Clean up the build
`rm -rf build` 

or maybe try  
`make clean`

## Mac OSX
On a Mac installation, you may have to change a couple relative paths to be absolute (when executing programs in an iPython or IDL environment):  
`install_name_tool -change libboost_thread.dylib /full/path/to/libboost_thread.dylib libraytracethread.dylib`  
`install_name_tool -change libboost_system.dylib /full/path/to/libboost_system.dylib libboost_thread.dylib`

## Dependencies
### C++ only
- boost (C++ library)
    - thread
    - unit_test_framework
    - bind.hpp
    
### Documentation
- doxygen
- dot/graphviz
- sphinx


## Generating the code documentation with Doxygen
Run  
`doxygen Doxyfile`

The doc can be viewed in your browser here:  
`doxygen/html/index.html`

## Gererating the user manual
The first user manual was writen in Docbook but we are working on rewriting it using Sphinx and ReStructuredText.

To generate the Sphinx documentation:  
`cd sphinx`  
`make html`

The doc can be viewed in your browser here:  
`sphinx/build/html/index.html`
