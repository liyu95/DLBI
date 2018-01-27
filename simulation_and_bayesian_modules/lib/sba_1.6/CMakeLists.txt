

include_directories(.)
# sba library source files
add_library(sba STATIC
  sba_levmar.c sba_levmar_wrap.c sba_lapack.c sba_crsm.c sba_chkjac.c
  sba.h sba_chkjac.h compiler.h
)

target_link_libraries(sba f2c clapack) #cblas

set_target_properties(sba PROPERTIES OUTPUT_NAME "sba")
set_target_properties(sba PROPERTIES CLEAN_DIRECT_OUTPUT 1)