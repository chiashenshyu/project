cmake_minimum_required(VERSION 2.8.12.2)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

##############################################################
project(PFLocalization)

set(MY_SOURCE_FILES
	${CMAKE_CURRENT_SOURCE_DIR}/common/src/dubins.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/common/src/dynamics.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/common/src/kdTreeNode.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/common/src/visualizer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PFLocalization/src/ParticleFilter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Controller/src/ppc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PathPlanning/Global_Planner/src/planner.cpp
   	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
   	)
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/common/include
	${CMAKE_CURRENT_SOURCE_DIR}/Controller/include
	${CMAKE_CURRENT_SOURCE_DIR}/PathPlanning/Global_Planner/include
	${CMAKE_CURRENT_SOURCE_DIR}/PFLocalization/include
	)

add_executable(${PROJECT_NAME} ${MY_SOURCE_FILES})

find_package(PythonLibs 2.7)
target_include_directories(${PROJECT_NAME} PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})

###############################################################
