cmake_minimum_required(VERSION 3.50)
project(RangeTree LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

option(USE_MULTI_THREADED_DLL "Use Multi-threaded DLL" ON)
option(RANGETREE_BUILD_TESTS "Build the tests" OFF)

set(CMAKE_CXX_STANDARD 11)

# Enable output of compile commands during generation.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set MSVC_RUNTIME_LIBRARY etc..
IF(MSVC)
    add_compile_options(/MP)
    if(USE_MULTI_THREADED_DLL)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
    else()
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
endif()

set(CMAKE_DEBUG_POSTFIX "d") 

if(RANGETREE_BUILD_TESTS)
  add_subdirectory(tests)
endif()

#install(FILES ${HEADERS} DESTINATION include)