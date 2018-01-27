include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(subtract main.cpp)

find_package(OpenCV REQUIRED)

target_link_libraries(subtract ${OpenCV_LIBS})
