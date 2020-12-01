cmake_minimum_required(VERSION 3.8)

project (example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(SOURCES example.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
