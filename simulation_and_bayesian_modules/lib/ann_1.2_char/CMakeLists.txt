
set(ANN_INCLUDE_DIR ./include)
set(ANN_SRC_DIR ./src)

aux_source_directory(./src ANN_SRCS)

MESSAGE(STATUS "building ann_1.2_char")

include_directories(${ANN_INCLUDE_DIR} ${ANN_SRC_DIR})

add_library (ann_1.2_char STATIC ${ANN_SRCS})

# target_link_libraries(ann_1.1_char)

set_target_properties(ann_1.2_char PROPERTIES OUTPUT_NAME "ann_1.1_char")
set_target_properties(ann_1.2_char PROPERTIES CLEAN_DIRECT_OUTPUT 1)