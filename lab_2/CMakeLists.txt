cmake_minimum_required(VERSION 3.16)
project("Calculate integral")

# set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(MYCOMPILE_FLAGS "-pthread")

SET(CMAKE_C_FLAGS    "${CMAKE_C_FLAGS} ${MYCOMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MYCOMPILE_FLAGS}")

file(GLOB SRC "src/*")

add_executable(lab_2 ${SRC})
target_include_directories(lab_2 PUBLIC "include")