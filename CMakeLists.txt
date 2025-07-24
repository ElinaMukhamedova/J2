cmake_minimum_required(VERSION 3.19.1)
project(J2)
set(CMAKE_CXX_STANDARD 20)

add_subdirectory(src)
add_subdirectory(third_party)

enable_testing()
add_subdirectory(tests)

add_executable(main main.cpp)
target_link_libraries(main J2)