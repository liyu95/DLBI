
MESSAGE(STATUS "building main")

add_executable(lapls main.cpp)

find_package(OpenCV REQUIRED)

target_link_libraries(lapls ${OpenCV_LIBS}) #simu  #opencv_core opencv_imgproc opencv_calib3d opencv_video opencv_features2d 
#	opencv_ml opencv_highgui opencv_objdetect opencv_contrib opencv_legacy)
