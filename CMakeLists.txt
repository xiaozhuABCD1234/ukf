cmake_minimum_required(VERSION 3.8)
project(main)
set(CMAKE_CXX_STANDARD 17)
add_executable(main main.cpp ukf.hpp ukf.cpp)
target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(main PRIVATE -Wall -Wextra -Wpedantic)
find_package(Eigen3 REQUIRED)
target_link_libraries(main PRIVATE Eigen3::Eigen)

