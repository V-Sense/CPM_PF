PROJECT(CPMPF_FLOW)

set(${PROJECT_NAME}_PROJECT_SRCS 
	main.cpp
	
	cpmpf_parameters.h
	cpmpf_parameters.cpp
	flowIO.h
	flowIO.cpp
	ImageIOpfm.h
	ImageIOpfm.cpp
	utils.h
	utils.cpp
)

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CPM_SOURCE_DIR}
    ${PFilter_SOURCE_DIR}
    ${Variational_refinement_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

ADD_EXECUTABLE(${PROJECT_NAME} ${${PROJECT_NAME}_PROJECT_SRCS})

TARGET_LINK_LIBRARIES(${PROJECT_NAME}
    ${OpenMP_CXX_LIBRARIES}
    ${OpenCV_LIBS}
    CPM
    PFilter
    Variational_refinement
)

