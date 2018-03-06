# SCRaytrace
Solar Corona Raytracing tools

## Compiling the source from github
SCraytrace now uses cmake. To compile the project once cloned from github, run the following:  
cmake -H. -Bbuild  
cmake --build build -- -j3

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

