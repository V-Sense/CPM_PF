cmake_minimum_required(VERSION 2.8)
set(MODULE_NAME PFilter)
project(${MODULE_NAME})

## Source files -----------------------------------------------------------------------------------

file(GLOB_RECURSE sources     "*.c" "*.cpp")
file(GLOB_RECURSE headers_hpp "*.hpp")
file(GLOB_RECURSE headers_h   "*.h")
list(APPEND headers ${headers_h} ${headers_hpp})


## Project setup ----------------------------------------------------------------------------------

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${OpenCV_INCLUDE_DIRS}
)

add_library(${MODULE_NAME} STATIC
        ${sources}
        ${headers}
)

SET_TARGET_PROPERTIES(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

TARGET_LINK_LIBRARIES(${MODULE_NAME}
    ${OpenCV_LIBS}
    )
