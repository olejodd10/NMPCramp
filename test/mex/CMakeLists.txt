add_test(NAME TestMexSimulateMmc_Horizons_1_2_5_10 COMMAND matlab -nodisplay -r 
    "addpath(\"${CMAKE_BINARY_DIR}/mex\"), \
     addpath(\"${CMAKE_CURRENT_SOURCE_DIR}\"), \
     TestMexSimulateMmc(\"${CMAKE_SOURCE_DIR}/test/data/MMCMEX\", 1, 100), \
     TestMexSimulateMmc(\"${CMAKE_SOURCE_DIR}/test/data/MMCMEX\", 2, 100), \
     TestMexSimulateMmc(\"${CMAKE_SOURCE_DIR}/test/data/MMCMEX\", 5, 100), \
     TestMexSimulateMmc(\"${CMAKE_SOURCE_DIR}/test/data/MMCMEX\", 10, 100), exit")

set(TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT 6e-3)

add_test(NAME TestMexSimulateMmc_Compare_xoutN1_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/MMCMEX/xoutN1.csv ${CMAKE_SOURCE_DIR}/test/data/MMC/xoutN1.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexSimulateMmc_Compare_uoutN1_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/MMCMEX/uoutN1.csv ${CMAKE_SOURCE_DIR}/test/data/MMC/uoutN1.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexSimulateMmc_Compare_xoutN2_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/MMCMEX/xoutN2.csv ${CMAKE_SOURCE_DIR}/test/data/MMC/xoutN2.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexSimulateMmc_Compare_uoutN2_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/MMCMEX/uoutN2.csv ${CMAKE_SOURCE_DIR}/test/data/MMC/uoutN2.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexSimulateMmc_Compare_xoutN5_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/MMCMEX/xoutN5.csv ${CMAKE_SOURCE_DIR}/test/data/MMC/xoutN5.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexSimulateMmc_Compare_uoutN5_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/MMCMEX/uoutN5.csv ${CMAKE_SOURCE_DIR}/test/data/MMC/uoutN5.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})

add_test(NAME TestMexSimulateMmc_Compare_xoutN10_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/MMCMEX/xoutN10.csv ${CMAKE_SOURCE_DIR}/test/data/MMC/xoutN10.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
add_test(NAME TestMexSimulateMmc_Compare_uoutN10_With_TestSimulateMmc COMMAND python ${CMAKE_SOURCE_DIR}/test/scripts/compare_csvs.py 
    ${CMAKE_SOURCE_DIR}/test/data/MMCMEX/uoutN10.csv ${CMAKE_SOURCE_DIR}/test/data/MMC/uoutN10.csv ${TEST_MEX_SIMULATE_MMC_COMPARISON_LIMIT})
