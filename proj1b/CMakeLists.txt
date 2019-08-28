cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project1B)
SET(VTK_DIR /Users/tcolby/Documents/CIS441/VTK-8.1.2)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(project1B MACOSX_BUNDLE project1B)

if(VTK_LIBRARIES)
  target_link_libraries(project1B ${VTK_LIBRARIES})
else()
  target_link_libraries(project1B vtkHybrid)
endif()

