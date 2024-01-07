# SCRaytrace

Solar Corona Raytracing tools

## Project tree structure

- src: C++ source code.
- data: data files used by the C++ code.
- python: python wrappers.
- sphinx: sphinx documentation.
- documentation: contains the user manual, in MS Word.

## Cloning the repository

`git clone https://github.com/DrRaytrace/SCRaytrace.git`

## Compiling the sources

SCraytrace now uses cmake. To compile the project once cloned from github, run the following:  
`cmake -H. -Bbuild`  
`cmake --build build -- -j3`

The compiled code will be generated in the `build` folder.  

For additional testing, you can run the following:  
`cd build`  
`make test`

Individual test:  
`cd build/src`  
`./testboosttest --log_level=message`  

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
- ntirpc

### Tip for compiling and installing Boost

Using Boost 1.72.0:  
`./bootstrap.sh --prefix=/Users/username/local --with-libraries=thread --with-libraries=test --with-libraries=date_time --with-libraries=headers`  
`./b2`  
`./b2 headers`  
`./b2 install`  

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

## Python wrapper

The Python module is in the python folder. It's name is scraytrace. It can be imported like so:  
`import scraytrace`  

The scraytrace Python wrapper can be tested by running the `test.py` in the `tests` sub-folder of `python/scraytrace`.

## News

### Version 3.1.0

- New dependency on ntirpc. This is to read xdr files. glibc does not provide this library anymore. Compiled using ntirpc: <https://github.com/linuxbox2/ntirpc.git>  

- Now uses libwcs to compute projections. This allows using pretty much all projection types implemented in libwcs.

- Use of libwcs implemented in the scraytrace python wrapper.
