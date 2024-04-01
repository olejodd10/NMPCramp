if (CMAKE_HOST_WIN32)
    message(STATUS "TestMexSimulateMmc configured for Windows")
    set(TEST_MEX_SIMULATE_MMC_MATLAB_FLAGS -nodesktop -nosplash -wait)
else()
    set(TEST_MEX_SIMULATE_MMC_MATLAB_FLAGS -nodisplay)
endif()

add_test(NAME TestMexSimulateMmc_Horizons_1_2_5_10 COMMAND matlab ${TEST_MEX_SIMULATE_MMC_MATLAB_FLAGS} -r 
    "addpath('${CMAKE_BINARY_DIR}/mex'), \
     addpath('${CMAKE_CURRENT_SOURCE_DIR}'), \
     TestMexSimulateMmc('${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc', 1, 10000), \
     clear mex, \
     TestMexSimulateMmc('${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc', 2, 10000), \
     clear mex, \
     TestMexSimulateMmc('${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc', 5, 10000), \
     clear mex, \
     TestMexSimulateMmc('${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc', 10, 10000), \
     exit")

set(TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT 8e-3)

add_test(NAME TestMexSimulateMmc_Compare_xoutN1_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc/xoutN1.csv ${CMAKE_SOURCE_DIR}/test/data/TestSimulateMmc/xoutN1.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexSimulateMmc_Compare_uoutN1_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc/uoutN1.csv ${CMAKE_SOURCE_DIR}/test/data/TestSimulateMmc/uoutN1.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexSimulateMmc_Compare_xoutN2_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc/xoutN2.csv ${CMAKE_SOURCE_DIR}/test/data/TestSimulateMmc/xoutN2.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexSimulateMmc_Compare_uoutN2_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc/uoutN2.csv ${CMAKE_SOURCE_DIR}/test/data/TestSimulateMmc/uoutN2.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexSimulateMmc_Compare_xoutN5_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc/xoutN5.csv ${CMAKE_SOURCE_DIR}/test/data/TestSimulateMmc/xoutN5.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexSimulateMmc_Compare_uoutN5_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc/uoutN5.csv ${CMAKE_SOURCE_DIR}/test/data/TestSimulateMmc/uoutN5.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexSimulateMmc_Compare_xoutN10_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc/xoutN10.csv ${CMAKE_SOURCE_DIR}/test/data/TestSimulateMmc/xoutN10.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexSimulateMmc_Compare_uoutN10_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/TestMexSimulateMmc/uoutN10.csv ${CMAKE_SOURCE_DIR}/test/data/TestSimulateMmc/uoutN10.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
