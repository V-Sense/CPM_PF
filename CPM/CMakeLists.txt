CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(MODULE_NAME CPM)
PROJECT (${MODULE_NAME})

file(GLOB_RECURSE sources     "*.c" "*.cpp")
file(GLOB_RECURSE headers_hpp "*.hpp")
file(GLOB_RECURSE headers_h   "*.h")
list(APPEND headers ${headers_h} ${headers_hpp})

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

ADD_LIBRARY(${MODULE_NAME} STATIC
        ${sources}
        ${headers}
)


TARGET_LINK_LIBRARIES(${MODULE_NAME}
    jpeg
    png
    ${OpenCV_LIBS}
    )

