cmake_minimum_required(VERSION 2.8.12)

# ---- Avoids warning with new policy when compiling on Mac OSX
set(CMAKE_MACOSX_RPATH 1)

project(scraytrace)
enable_testing()


add_subdirectory(src)
add_subdirectory(data)
#configure_file(data ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)


set (scraytrace_VERSION_MAJOR 3)
set (scraytrace_VERSION_MINOR 1)
set (scraytrace_VERSION_PATCH 0)




# ---- packaging
# include (InstallRequiredSystemLibraries)

# set (CPACK_RESOURCE_FILE_LICENSE  
#      "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")

# set(CPACK_PACKAGE_VERSION_MAJOR "${scraytrace_VERSION_MAJOR}")
# set(CPACK_PACKAGE_VERSION_MINOR "${scraytrace_VERSION_MINOR}")
# set(CPACK_PACKAGE_VERSION_PATCH "${scraytrace_VERSION_PATCH}")
# set(CPACK_SOURCE_GENERATOR "TBZ2")


# set(CPACK_SOURCE_PACKAGE_FILE_NAME
#   "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
# set(CPACK_SOURCE_IGNORE_FILES
#   "/build/;/.bzr/;~$;${CPACK_SOURCE_IGNORE_FILES}")

# include(CPack)

# Run this way to create a package
# For binnary:
#  cpack -C CPackConfig.cmake
# For source:
#  cpack -C CPackSourceConfig.cmake

