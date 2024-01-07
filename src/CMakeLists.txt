
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 


set(Boost_USE_STATIC_LIBS OFF)
set(BUILD_SHARED_LIBS ON)

find_package(Boost 1.67.0 REQUIRED COMPONENTS thread unit_test_framework)

include_directories(~/Documents/cpp/ntirpc-master/ntirpc) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    message("Boost include directory: ${Boost_INCLUDE_DIRS}")
endif()

find_library(BoostThread libboost_thread.dylib libboost_thread.so HINTS ~/local/lib ~/.local/lib ~/soft/lib )
find_library(BoostTest libboost_unit_test_framework.dylib libboost_unit_test_framework.so HINTS ~/local/lib ~/.local/lib ~/soft/lib)
find_library(NTIRPC_LIBRARY libntirpc.so libntirpc.dylib HINTS ~/.local/lib REQUIRED)
find_library(WCSLIB_LIBRARY libwcs.so)


if(NTIRPC_LIBRARY)
    message(STATUS "Found NTIRPC: ${NTIRPC_LIBRARY}")
else()
    message("Did not find NTIRPC")
endif() 

message(STATUS "BoostThread: ${BoostThread}")

# ---- SCRaytrace requires boost/bind.hpp
find_file(BoostBind NAMES bind.hpp PATHS
      ${Boost_INCLUDE_DIRS}/boost)

if(BoostBind)
    message(STATUS "Found boost/bind.hpp")
else()
    message(FATAL_ERROR "${PROJECT_NAME} requires boost/bind.hpp. Could not find it in ${Boost_INCLUDE_DIRS}")
endif()


# ---- SCRaytrace requires boost/test/unit_test.hpp
find_file(BoostUnitTest NAMES unit_test.hpp PATHS
      ${Boost_INCLUDE_DIRS}/boost/test)

if(BoostUnitTest)
    message(STATUS "Found boost/test/unit_test.hpp")
else()
    message(FATAL_ERROR "${PROJECT_NAME} requires boost/test/unit_test.hpp")
endif()


set(SCRAYTRACE_DATA_DIR "${CMAKE_SOURCE_DIR}/data")

# ---- define configure header file
configure_file(
    "${PROJECT_SOURCE_DIR}/src/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )
include_directories(${PROJECT_BINARY_DIR})



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
target_link_libraries(raytracethread raytrace  Boost::thread ${NTIRPC_LIBRARY} ${WCSLIB_LIBRARY})



add_executable(raytracerun main.cpp)
target_link_libraries(raytracerun raytrace ${WCSLIB_LIBRARY})

add_executable(rtthread rtthreadmain.cpp)
target_link_libraries(rtthread raytracethread ${WCSLIB_LIBRARY})

add_executable(hpsmie oonewcallhpsmie.cpp)
target_link_libraries(hpsmie raytrace raytracethread Boost::thread ${WCSLIB_LIBRARY})

add_executable(testwcs testwcs.cpp)
target_link_libraries(testwcs raytracethread ${WCSLIB_LIBRARY})

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

