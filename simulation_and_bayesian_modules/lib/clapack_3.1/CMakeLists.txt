
set(F2C_INCLUDE_DIR ./F2CLIBS)
set(F2C_SRC_DIR ./F2CLIBS/libf2c)

aux_source_directory(${F2C_SRC_DIR} F2C_SRC)

MESSAGE(STATUS "building f2clibs")

include_directories(${F2C_INCLUDE_DIR} ${F2C_SRC_DIR})

add_library (f2c STATIC ${F2C_INCLUDE_DIR}/f2c.h ${F2C_SRC})

target_link_libraries(f2c m)

ADD_DEFINITIONS(
    -DINTEGER_STAR_8 
    #-DKR_headers
)

set_target_properties(f2c PROPERTIES OUTPUT_NAME "f2c")
set_target_properties(f2c PROPERTIES CLEAN_DIRECT_OUTPUT 1)
############################################################################################

set(BLASWRAP_INCLUDE_DIR ./BLAS/WRAP)
set(BLASWRAP_SRC_DIR ./BLAS/SRC)

aux_source_directory(${BLASWRAP_SRC_DIR} BLASWRAP_SRCS)

MESSAGE(STATUS "building blaswrap")

include_directories(${BLASWRAP_SRC_DIR} ${BLASWRAP_INCLUDE_DIR})

add_library (blaswrap STATIC ${BLASWRAP_INCLUDE_DIR}/blaswrap.h ${BLASWRAP_SRCS})

target_link_libraries(blaswrap f2c m)

set_target_properties(blaswrap PROPERTIES OUTPUT_NAME "blaswrap")
set_target_properties(blaswrap PROPERTIES CLEAN_DIRECT_OUTPUT 1)
############################################################################################

set(CLAPACK_INCLUDE_DIR ./INCLUDE)
set(CLAPACK_SRC_DIR ./SRC)

aux_source_directory(${CLAPACK_SRC_DIR} CLAPACK_SRC)
aux_source_directory(./INSTALL CLAPACK_SRC2)

MESSAGE(STATUS "building clapack_3.1")

include_directories(${CLAPACK_INCLUDE_DIR} ${CLAPACK_SRC_DIR})

add_library (clapack STATIC ${CLAPACK_INCLUDE_DIR}/clapack.h ${CLAPACK_SRC} ${CLAPACK_SRC2})

target_link_libraries(clapack blaswrap m)

set_target_properties(clapack PROPERTIES OUTPUT_NAME "clapack")
set_target_properties(clapack PROPERTIES CLEAN_DIRECT_OUTPUT 1)