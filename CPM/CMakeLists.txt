CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(MODULE_NAME CPM)
PROJECT (${MODULE_NAME})

#ADD_DEFINITIONS(-DWITH_SSE)

#find_package(OpenCV REQUIRED)
#find_package(LAPACK REQUIRED)

file(GLOB_RECURSE sources     "*.c" "*.cpp")
file(GLOB_RECURSE headers_hpp "*.hpp")
file(GLOB_RECURSE headers_h   "*.h")
list(APPEND headers ${headers_h} ${headers_hpp})

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
#C:/OpenCV/include
#C:/OpenCV/include/opencv
#C:/OpenCV/include/opencv2
)

LINK_DIRECTORIES(
#C:/OpenCV/lib/x86/Debug
#C:/OpenCV/lib/x86/Release
)

#ADD_EXECUTABLE(${MODULE_NAME} main.cpp CPM.cpp)


ADD_LIBRARY(${MODULE_NAME} STATIC
        ${sources}
        ${headers}
)


TARGET_LINK_LIBRARIES(${MODULE_NAME}
    jpeg
    png
    ${OpenCV_LIBS}
    )


#debug opencv_core310d.lib
#debug opencv_highgui310d.lib
#debug opencv_imgcodecs310d.lib
#debug opencv_imgproc310d.lib
#debug opencv_xfeatures2d310d.lib
#optimized opencv_core310.lib
#optimized opencv_highgui310.lib
#optimized opencv_imgcodecs310.lib
#optimized opencv_imgproc310.lib
#optimized opencv_xfeatures2d310.lib)
