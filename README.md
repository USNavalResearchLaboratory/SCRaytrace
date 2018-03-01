# SCRaytrace
Solar Corona Raytracing tools

SCraytrace now uses cmake. To compile the project once cloned from github, run the following:
cmake -H. -Bbuild
cmake --build build -- -j3

For additional testing, you can run the following:
cd build
make test

Dependencies
- boost C++ library
- cppunit

