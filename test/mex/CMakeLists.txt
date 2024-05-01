if (CMAKE_HOST_WIN32)
    message(STATUS "TestMexMmcMpc configured for Windows")
    set(TEST_MEX_MMC_MPC_MATLAB_FLAGS -nodesktop -nosplash -wait)
else()
    set(TEST_MEX_MMC_MPC_MATLAB_FLAGS -nodisplay)
endif()

add_test(NAME TestMexMmcMpc_Horizons_1_2_5_10 COMMAND matlab ${TEST_MEX_MMC_MPC_MATLAB_FLAGS} -r 
    "addpath('${CMAKE_BINARY_DIR}/mex'), \
     addpath('${CMAKE_CURRENT_SOURCE_DIR}'), \
     TestMexMmcMpc('${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc', 1, 10000), \
     clear mex, \
     TestMexMmcMpc('${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc', 2, 10000), \
     clear mex, \
     TestMexMmcMpc('${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc', 5, 10000), \
     clear mex, \
     TestMexMmcMpc('${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc', 10, 10000), \
     exit")

set(TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT 8e-3)

add_test(NAME TestMexMmcMpc_Compare_xoutN1_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/xoutN1.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/xoutN1.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexMmcMpc_Compare_uoutN1_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/uoutN1.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/uoutN1.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexMmcMpc_Compare_xoutN2_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/xoutN2.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/xoutN2.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexMmcMpc_Compare_uoutN2_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/uoutN2.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/uoutN2.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexMmcMpc_Compare_xoutN5_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/xoutN5.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/xoutN5.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexMmcMpc_Compare_uoutN5_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/uoutN5.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/uoutN5.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexMmcMpc_Compare_xoutN10_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/xoutN10.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/xoutN10.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexMmcMpc_Compare_uoutN10_With_TestMmcMpc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexMmcMpc/uoutN10.csv ${CMAKE_SOURCE_DIR}/test/data/TestMmcMpc/uoutN10.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
