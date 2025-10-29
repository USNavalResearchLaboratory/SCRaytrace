
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 


set(Boost_USE_STATIC_LIBS OFF)
set(BUILD_SHARED_LIBS ON)
find_package(Boost 1.67 CONFIG REQUIRED COMPONENTS thread unit_test_framework)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    message("Boost include directory: ${Boost_INCLUDE_DIRS}")
endif()

find_library(BoostThread libboost_thread.dylib libboost_thread.so HINTS ~/local/lib ~/soft/lib)
find_library(BoostTest libboost_unit_test_framework.dylib libboost_unit_test_framework.so HINTS ~/local/lib ~/soft/lib)

# ---- SCRaytrace requires boost/bind.hpp
find_file(BoostBind NAMES bind.hpp PATHS
      ${Boost_INCLUDE_DIRS}/boost)

if(BoostBind)
    message("Found boost/bind.hpp")
else()
    message(FATAL_ERROR "${PROJECT_NAME} requires boost/bind.hpp. Could not find it in ${Boost_INCLUDE_DIRS}")
endif()


# ---- SCRaytrace requires boost/test/unit_test.hpp
find_file(BoostUnitTest NAMES unit_test.hpp PATHS
      ${Boost_INCLUDE_DIRS}/boost/test)

if(BoostUnitTest)
    message("Found boost/test/unit_test.hpp")
else()
    message(FATAL_ERROR "${PROJECT_NAME} requires boost/test/unit_test.hpp")
endif()

# Add tirpc include directory and library
#include_directories(/usr/include/tirpc)
#link_libraries(tirpc)
find_library(TIRPC_LIBRARY NAMES tirpc)
find_path(TIRPC_INCLUDE_DIR NAMES tirpc/rpc/rpc.h)


set(SCRAYTRACE_DATA_DIR "${CMAKE_SOURCE_DIR}/data")

# ---- define configure header file
configure_file(
    "${PROJECT_SOURCE_DIR}/src/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )
include_directories(${PROJECT_BINARY_DIR})


# Check for GCC or Clang
message("Compiler is ${CMAKE_CXX_COMPILER_ID}")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    # Common flags for GCC and Clang
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

    # Specific flags for GCC
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstack-protector-strong")
    endif()

    # Specific flags for Clang
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics -std=c++20")
    endif()

# Check for MSVC
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Flags for MSVC
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /WX")

# Handle other compilers or provide default flags
else()
    message(STATUS "Unknown C++ compiler: ${CMAKE_CXX_COMPILER_ID}. Using default flags.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2") # Example default flag
endif()

# You can also set flags based on build type (Debug/Release)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DDEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
endif()

message("CMAKE_CXX_FLAGS is set to ${CMAKE_CXX_FLAGS}")


add_library(raytrace SHARED Cmat.cpp    Cvec.cpp    Cbasis.cpp  Cbasis.h    Clos.cpp    Clos.h 
    CModelBase.cpp  CModelBase.h    raytrace.cpp    Cvec.h  Cmat.h  constant.h 
    raytrace.h  CControlPoint.h CControlPoint.cpp   
    Cnurbs.h    Cnurbs.cpp  CModel46.cpp    CModel46.h  CModel48.cpp CModel48.h 
    models01to10.cpp    models01to10.h  models11to20.cpp    models11to20.h  
    models21to30.cpp    models21to30.h  models31to40.cpp    models31to40.h  
    models41to50.cpp    models41to50.h  models51to60.cpp    models51to60.h  
    models61to70.cpp    models61to70.h  models71to80.cpp    models71to80.h
    models81to90.cpp    models81to90.h
    Cdvoxel.cpp Cdvoxel.h Cvoxvert.cpp Cvoxvert.h 
    rtcloud.cpp rtcloud.h 
    rtmiscfunc.h rtmiscfunc.cpp voxinteg.h ModelPosition.h)

add_library(raytracethread SHARED scene.h scene.cpp rtthread.h rtthread.cpp
    camera.h camera.cpp physicsbase.h physicsbase.cpp 
    physicsthomson.h physicsthomson.cpp physicsuv.h physicsuv.cpp 
    physicsisotropic.h physicsisotropic.cpp physicsmie.h physicsmie.cpp 
    miehps.h miehps.cpp physicsvsf.cpp physicsvsf.h 
    physicsvsfvarydist.cpp physicsvsfvarydist.h physicsvariablevsf.cpp physicsvariablevsf.h)

target_link_libraries(raytracethread raytrace  Boost::thread)

# Link with tirpc only if it was found
if(TIRPC_LIBRARY)
    target_link_libraries(raytracethread raytrace ${TIRPC_LIBRARY})
    # If you also checked for TIRPC_INCLUDE_DIR:
    # target_include_directories(myapp PRIVATE ${TIRPC_INCLUDE_DIR})
endif()


add_executable(raytracerun main.cpp)
target_link_libraries(raytracerun raytrace)

add_executable(rtthread rtthreadmain.cpp)
target_link_libraries(rtthread raytracethread)

add_executable(hpsmie oonewcallhpsmie.cpp)
target_link_libraries(hpsmie raytrace raytracethread Boost::thread)

# ---- tests
add_executable(testmodel26 model26test.cpp)
target_link_libraries(testmodel26 raytracethread Boost::unit_test_framework)

add_executable(testmodel67 model67test.cpp )
target_link_libraries(testmodel67 raytracethread Boost::unit_test_framework)

add_executable(testmodel77 model77test.cpp )
target_link_libraries(testmodel77 raytracethread Boost::unit_test_framework)

add_executable(testrtmiscfunc rtmiscfunctest.cpp  rtmiscfunc.h)
target_link_libraries(testrtmiscfunc raytracethread Boost::unit_test_framework)

add_executable(testcvec cvectest.cpp)
target_link_libraries(testcvec raytracethread Boost::unit_test_framework)

add_executable(testcamera cameratest.cpp)
target_link_libraries(testcamera raytracethread Boost::unit_test_framework)

add_executable(testCbasis Cbasistest.cpp)
target_link_libraries(testCbasis raytracethread Boost::unit_test_framework)

add_executable(testscene scenetest.cpp)
target_link_libraries(testscene raytracethread Boost::unit_test_framework)

add_executable(testphysics physicstest.cpp)
target_link_libraries(testphysics raytracethread Boost::unit_test_framework)


add_executable(testboosttest boostTestTest.cpp)
target_link_libraries(testboosttest Boost::unit_test_framework) # ${BoostTest}



add_test(NAME testmodel26 COMMAND testmodel26)
add_test(NAME testmodel67 COMMAND testmodel67)
add_test(NAME testmodel77 COMMAND testmodel77)
add_test(NAME testrtmiscfunc COMMAND testrtmiscfunc)
add_test(NAME testcvec COMMAND testcvec)
add_test(NAME testcamera COMMAND testcamera)
add_test(NAME testCbasis COMMAND testCbasis)
add_test(NAME testscene COMMAND testscene)
add_test(NAME testphysics COMMAND testphysics)
add_test(NAME boostTestTest COMMAND testboosttest)

