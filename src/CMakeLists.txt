


set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(Boost 1.49.0 COMPONENTS thread)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

find_library(BoostThread libboost_thread.so HINTS ~/soft/lib)


# ---- find cppunit
find_path(CPPUNIT_INCLUDE_DIR cppunit/TestCase.h
    /usr/local/include
    /usr/include
    ~/soft/include)

find_library(CPPUNIT_LIBRARY cppunit
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)

find_library(CPPUNIT_DEBUG_LIBRARY cppunit
               ${CPPUNIT_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)

if(CPPUNIT_INCLUDE_DIR)
  if(CPPUNIT_LIBRARY)
    message("Found cppunit")
    set(CPPUNIT_FOUND "YES")
    set(CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY} ${CMAKE_DL_LIBS})
    set(CPPUNIT_DEBUG_LIBRARIES ${CPPUNIT_DEBUG_LIBRARY}
${CMAKE_DL_LIBS})
  else(CPPUNIT_LIBRARY)
    message("cppunit not found")
  endif(CPPUNIT_LIBRARY)
endif(CPPUNIT_INCLUDE_DIR)

set(SCRAYTRACE_DATA_DIR "${CMAKE_SOURCE_DIR}/data")

# ---- define configure header file
configure_file(
    "${PROJECT_SOURCE_DIR}/src/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )
include_directories(${PROJECT_BINARY_DIR})


# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wl,-rpath,$ORIGIN")


add_library(raytrace SHARED Cmat.cpp    Cvec.cpp    Cbasis.cpp  Cbasis.h    Clos.cpp    Clos.h 
    CModelBase.cpp  CModelBase.h    raytrace.cpp    Cvec.h  Cmat.h  constant.h 
    raytrace.h  CControlPoint.h CControlPoint.cpp   
    Cnurbs.h    Cnurbs.cpp  CModel46.cpp    CModel46.h  CModel48.cpp CModel48.h 
    models01to10.cpp    models01to10.h  models11to20.cpp    models11to20.h  
    models21to30.cpp    models21to30.h  models31to40.cpp    models31to40.h  
    models41to50.cpp    models41to50.h  models51to60.cpp    models51to60.h  
    models61to70.cpp    models61to70.h  models71to80.cpp    models71to80.h 
    Cdvoxel.cpp Cdvoxel.h Cvoxvert.cpp Cvoxvert.h 
    rtcloud.cpp rtcloud.h cuvemission.cpp cuvemission.h 
    rtmiscfunc.h rtmiscfunc.cpp voxinteg.h ModelPosition.h)

# add_library(raytraceStatic STATIC Cmat.cpp    Cvec.cpp    Cbasis.cpp  Cbasis.h    Clos.cpp    Clos.h 
#     CModelBase.cpp  CModelBase.h    raytrace.cpp    Cvec.h  Cmat.h  constant.h 
#     raytrace.h  CControlPoint.h CControlPoint.cpp   
#     Cnurbs.h    Cnurbs.cpp  CModel46.cpp    CModel46.h  CModel48.cpp CModel48.h 
#     models01to10.cpp    models01to10.h  models11to20.cpp    models11to20.h  
#     models21to30.cpp    models21to30.h  models31to40.cpp    models31to40.h  
#     models41to50.cpp    models41to50.h  models51to60.cpp    models51to60.h  
#     models61to70.cpp    models61to70.h  models71to80.cpp    models71to80.h 
#     Cdvoxel.cpp Cdvoxel.h Cvoxvert.cpp Cvoxvert.h 
#     rtcloud.cpp rtcloud.h cuvemission.cpp cuvemission.h 
#     rtmiscfunc.h rtmiscfunc.cpp voxinteg.h ModelPosition.h)
# set_target_properties(raytraceStatic PROPERTIES OUTPUT_NAME raytrace_static)


add_library(raytracethread SHARED scene.h scene.cpp rtthread.h rtthread.cpp 
    sun.h sun.cpp camera.h camera.cpp physicsbase.h physicsbase.cpp 
    physicsthomson.h physicsthomson.cpp physicsuv.h physicsuv.cpp 
    physicsisotropic.h physicsisotropic.cpp physicsmie.h physicsmie.cpp 
    miehps.h miehps.cpp physicsvsf.cpp physicsvsf.h 
    physicsvsfvarydist.cpp physicsvsfvarydist.h)
target_link_libraries(raytracethread raytrace ${BoostThread})
# target_link_libraries(raytracethread raytrace "-Wl,--whole-archive" ${Boost_THREAD_LIBRARY} "-Wl,--no-whole-archive")
# target_link_libraries(raytracethread raytrace  ${Boost_THREAD_LIBRARY} )

# link flag -Wl,-rpath,$ORIGIN



add_executable(raytracerun main.cpp)
target_link_libraries(raytracerun raytrace)

add_executable(rtthread rtthreadmain.cpp)
target_link_libraries(rtthread raytracethread)
# target_link_libraries(rtthread ${Boost_THREAD_LIBRARY})
# target_link_libraries(rtthread "-Wl,--whole-archive" ${Boost_THREAD_LIBRARY} "-Wl,--no-whole-archive")

# add_executable(hpsmie oonewcallhpsmie.cpp)
# target_link_libraries(hpsmie raytrace raytracethread "-Wl,--whole-archive" ${Boost_THREAD_LIBRARY} "-Wl,--no-whole-archive")
# target_link_libraries(hpsmie raytrace raytracethread "-Wl,--whole-archive" ${Boost_THREAD_LIBRARY} "-Wl,--no-whole-archive")

# ---- tests
add_executable(testcuvemission cuvemissiontest.cpp cuvemissiontest.h testsuite.cpp)
target_link_libraries(testcuvemission raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testmodel26 model26test.cpp model26test.h testsuite.cpp)
target_link_libraries(testmodel26 raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testmodel67 model67test.cpp model67test.h testsuite.cpp)
target_link_libraries(testmodel67 raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testmodel77 model77test.cpp model77test.h testsuite.cpp)
target_link_libraries(testmodel77 raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testrtmiscfunc rtmiscfunctest.cpp rtmiscfunctest.h rtmiscfunc.h testsuite.cpp)
target_link_libraries(testrtmiscfunc raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testcvec cvectest.h cvectest.cpp testsuite.cpp)
target_link_libraries(testcvec raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testsun suntest.h suntest.cpp testsuite.cpp)
target_link_libraries(testsun raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testcamera cameratest.h cameratest.cpp testsuite.cpp)
target_link_libraries(testcamera raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testCbasis Cbasistest.h Cbasistest.cpp testsuite.cpp)
target_link_libraries(testCbasis raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testscene scenetest.h scenetest.cpp testsuite.cpp)
target_link_libraries(testscene raytrace raytracethread ${CPPUNIT_LIBRARIES})

add_executable(testphysics physicstest.h physicstest.cpp  testsuite.cpp)
target_link_libraries(testphysics raytrace raytracethread ${CPPUNIT_LIBRARIES})


add_test(NAME testcuvemission COMMAND testcuvemission)
add_test(NAME testmodel26 COMMAND testmodel26)
add_test(NAME testmodel67 COMMAND testmodel67)
add_test(NAME testmodel77 COMMAND testmodel77)
add_test(NAME testrtmiscfunc COMMAND testrtmiscfunc)
add_test(NAME testcvec COMMAND testcvec)
add_test(NAME testsun COMMAND testsun)
add_test(NAME testcamera COMMAND testcamera)
add_test(NAME testCbasis COMMAND testCbasis)
add_test(NAME testscene COMMAND testscene)
add_test(NAME testphysics COMMAND testphysics)



