cmake_minimum_required(VERSION 3.1)
project(isONclust2)

set(isONclust2_VERSION_MAJOR 2)
set(isONclust2_VERSION_MINOR 3)

set(CMAKE_MACOSX_RPATH 0)
#set(CMAKE_SYSTEM_PROCESSOR x86-64)

option(BUILD_SHARED_LIBS OFF)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

set (CPACK_PACKAGE_VERSION "\"v${isONclust2_VERSION_MAJOR}.${isONclust2_VERSION_MINOR}-${GIT_COMMIT_HASH}\"")

configure_file(
    "${PROJECT_SOURCE_DIR}/src/isONclust2_config.h.in"
        "${PROJECT_BINARY_DIR}/isONclust2_config.h"
    )

include_directories("${PROJECT_BINARY_DIR}")
find_package(Threads REQUIRED)
set(HAVE_POSIX_MEMALIGN 1)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static -lpthread -lm")
    add_compile_options(
        -msse3
        # All warnings
        -Wall
        #-fopt-info-vec-all
        -Wno-unused-local-typedefs
        -Wno-unused-variable
        -Wno-strict-overflow
        -Wno-error=maybe-uninitialized
    )
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR CMAKE_COMPILER_IS_CLANGXX)
    set(CMAKE_EXE_LINKER_FLAGS "-lc++ -lc++abi -lm -lc")    
    ##--------
    #   - Added TBB_USE_GLIBCXX_VERSION macro to specify the version of GNU
    #     libstdc++ when it cannot be properly recognized, e.g. when used
    #     with Clang on Linux* OS.
    string(REPLACE "." "0" TBB_USE_GLIBCXX_VERSION ${CMAKE_CXX_COMPILER_VERSION})
    add_definitions(-DTBB_USE_GLIBCXX_VERSION=${TBB_USE_GLIBCXX_VERSION})
    ##-------
        add_compile_options(
    -msse
    -msse2
    -msse3
    -msse4.1
    -mavx
    -mavx2
        # All warnings, error-on-warn
        -Wall 
        -Wno-unknown-pragmas
        -Wno-unused-variable
    -Wno-tautological-pointer-compare
    )
    if(APPLE)
        add_compile_options(
        )
    endif()
endif()

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/vendor/bioparser/vendor/zlib)
include_directories(${PROJECT_SOURCE_DIR}/vendor/bioparser/include)
include_directories(${PROJECT_SOURCE_DIR}/vendor/tbb/include)
include_directories(${PROJECT_SOURCE_DIR}/vendor/parasail)
include_directories(${PROJECT_SOURCE_DIR}/vendor/spoa)
include_directories(${PROJECT_SOURCE_DIR}/vendor/spoa/vendor/cereal/include)

include_directories(${PROJECT_BINARY_DIR}/vendor/bioparser/vendor/zlib)
add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/tbb EXCLUDE_FROM_ALL)
add_subdirectory(vendor/parasail EXCLUDE_FROM_ALL)
add_subdirectory(vendor/spoa EXCLUDE_FROM_ALL)


message(STATUS "CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")

set(SOURCE 
    src/main.cpp
    src/args.cpp
    src/seq.cpp
    src/qualscore.cpp
    src/output.cpp
    src/hpc.cpp
    src/minimizer.cpp
    src/kmer_index.cpp
    src/p_emp_prob.cpp
    src/util.cpp
    src/pbar.cpp
    src/cluster.cpp
    src/serialize.cpp
    src/consensus.cpp
    )

add_library(lisONclust2
    src/args.cpp
    src/seq.cpp
    src/qualscore.cpp
    src/output.cpp
    src/hpc.cpp
    src/minimizer.cpp
    src/kmer_index.cpp
    src/p_emp_prob.cpp
    src/util.cpp
    src/pbar.cpp
    src/cluster.cpp
    src/serialize.cpp
    src/consensus.cpp
)
target_link_libraries(lisONclust2 tbb_static Threads::Threads parasail spoa)

add_executable(isONclust2 ${SOURCE})
set_target_properties(isONclust2 PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(isONclust2 PROPERTIES LINK_SEARCH_END_STATIC 1)
target_link_libraries(isONclust2 tbb_static bioparser parasail cereal Threads::Threads spoa)


add_executable(test_isONclust2
        test/isONclust2_test.cpp)
set_target_properties(test_isONclust2 PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(test_isONclust2 PROPERTIES LINK_SEARCH_END_STATIC 1)

add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)
target_link_libraries(test_isONclust2 bioparser parasail tbb_static lisONclust2 gtest_main Threads::Threads spoa)
