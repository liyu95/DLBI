project(sr-simulator)

cmake_minimum_required(VERSION 2.8)

#the store of data in clapack or cblas is column main sequence

set(ANN_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/ann_1.1_char/include)
set(CLAPACK_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/clapack_3.1/INCLUDE)
set(MRC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/mrc)
set(BLASWRAP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/clapack_3.1/BLAS/WRAP)
set(CBLAS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/cblas/INCLUDE)
set(SBA_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/sba_1.6)
set(MATRIX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/matrix)
set(OPENCV_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/opencv_2.3.1 
    ${PROJECT_SOURCE_DIR}/lib/opencv_2.3.1/opencv ${PROJECT_SOURCE_DIR}/lib/opencv_2.3.1/opencv2 )
set(CMINPACK_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/cminpack)



include_directories(${PROJECT_SOURCE_DIR}/lib ${ANN_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src 
    ${CLAPACK_INCLUDE_DIR} ${BLASWRAP_INCLUDE_DIR} ${CBLAS_INCLUDE_DIR} ${MRC_INCLUDE_DIR} ${SBA_INCLUDE_DIR}
    ${MATRIX_INCLUDE_DIR} ${OPENCV_INCLUDE_DIR} ${CMINPACK_INCLUDE_DIR})
link_directories (${PROJECT_BINARY_DIR}/lib)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(lib)
add_subdirectory(src)
