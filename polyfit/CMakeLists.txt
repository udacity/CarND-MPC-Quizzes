project(q1)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp)

include_directories(/usr/local/include)
include_directories(src/Eigen-3.3)
link_directories(/usr/local/lib)

add_executable(polyfit ${sources})

target_link_libraries(polyfit)
